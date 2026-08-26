using Plots

# Lamp-post line profile via the file-based kerrz C API
# (`krz_tool_Lineprofile_run`). Loads a precomputed emissivity FITS file,
# builds transfer functions, and returns (g, flux) in memory.
#
# Emissivity files: kerrz-modified/emissivity/emis_a{spin}_h{height}.fits
# Missing files are generated with the kerrz CLI (Γ=2.0 is baked into the FITS).
# Do not use a=0 — kerrz panics in the elliptic integrals; use a small positive spin.

const kerrz_root = joinpath(@__DIR__, "kerrz-modified")
const libkerrz = joinpath(kerrz_root, "zig-out", "lib", "libkerrz.dylib")
const kerrz_cli = joinpath(kerrz_root, "zig-out", "bin", "kerrz")
const emissivity_dir = joinpath(kerrz_root, "emissivity")

# ---------------------------------------------------------------------
# Structs mirroring kerrz-modified/wrappers/kerrz.h
# ---------------------------------------------------------------------
struct KrzKerrMetric
    M::Cdouble
    a::Cdouble
    horizon_radius::Cdouble
    horizon_radius_negative::Cdouble
    isco::Cdouble
end

struct KrzFourVector
    t::Cdouble
    r::Cdouble
    th::Cdouble
    ph::Cdouble
end

struct KrzThreadPool
    pool::Ptr{Cvoid}
    num_threads::Csize_t
end

struct KrzToolLineprofile
    metric::KrzKerrMetric
    x_obs::KrzFourVector
    emissivity_path::Ptr{Cchar}
    r_in::Cdouble
    r_out::Cdouble
    nradii::Csize_t
    nangles::Csize_t
    ng::Csize_t
    ngstar::Csize_t
    nrsteps::Csize_t
    normalise::Cint
    _pad::Cint
end

struct KrzLineprofile
    g_grid::Ptr{Cdouble}
    flux::Ptr{Cdouble}
    num_g::Csize_t
end

@assert sizeof(KrzKerrMetric) == 40
@assert sizeof(KrzFourVector) == 32
@assert sizeof(KrzThreadPool) == 16
@assert sizeof(KrzToolLineprofile) == 144
@assert sizeof(KrzLineprofile) == 24

const KRZ_RETCODE_NAMES = ("SUCCESS", "ALLOCATION_FAILED", "THREAD_ERROR",
                           "CONVERGENCE_FAILED", "INVALID_ARGUMENT")

function _check_retcode(ret::Integer, fname::AbstractString)
    ret == 0 && return
    name = 0 <= ret < length(KRZ_RETCODE_NAMES) ? KRZ_RETCODE_NAMES[ret + 1] : "UNKNOWN"
    error("$fname failed with code $ret ($name)")
end

# ---------------------------------------------------------------------
# Path helper and on-demand emissivity generation
# ---------------------------------------------------------------------
function _height_str(h)
    isinteger(h) ? string(Int(h)) : string(h)
end

function emissivity_path(a, h)
    joinpath(emissivity_dir, "emis_a$(a)_h$(_height_str(h)).fits")
end

"""
    ensure_emissivity(a, h) -> path

Return the FITS path for `(a, h)`, running `kerrz emissivity` if the file
is missing (same flags as `scripts/precompute_emissivity.sh`).
"""
function ensure_emissivity(a, h; nphotons::Integer = 3000,
                           photon_index::Float64 = 2.0)
    path = emissivity_path(a, h)
    isfile(path) && return path
    isfile(kerrz_cli) || error(
        "kerrz CLI not found at $kerrz_cli\nBuild with: (cd kerrz-modified && zig build)")
    mkpath(emissivity_dir)
    println("Generating emissivity a=$a h=$h -> $path")
    run(`$kerrz_cli emissivity --spin $a --lamppost h:$(_height_str(h)),vr:0 --photon-index $photon_index --nphotons $nphotons -o $path`)
    isfile(path) || error("kerrz emissivity did not write $path")
    return path
end

# ---------------------------------------------------------------------
# ccall wrappers
# ---------------------------------------------------------------------
function make_metric(M::Float64, a::Float64)
    ccall((:krz_kerrMetric, libkerrz), KrzKerrMetric, (Cdouble, Cdouble), M, a)
end

function threadpool_init!(pool::Ref{KrzThreadPool}, n_threads::Integer)
    ret = ccall((:krz_ThreadPool_init, libkerrz), Cint,
                (Ptr{KrzThreadPool}, Csize_t), pool, n_threads)
    _check_retcode(ret, "krz_ThreadPool_init")
    return pool
end

threadpool_deinit!(pool::Ref{KrzThreadPool}) =
    ccall((:krz_ThreadPool_deinit, libkerrz), Cvoid, (Ptr{KrzThreadPool},), pool)

lineprofile_deinit!(lp::Ref{KrzLineprofile}) =
    ccall((:krz_Lineprofile_deinit, libkerrz), Cvoid, (Ptr{KrzLineprofile},), lp)

function _copy_c_doubles(ptr::Ptr{Cdouble}, n::Integer)
    (ptr == C_NULL || n == 0) && return Float64[]
    copy(unsafe_wrap(Vector{Float64}, ptr, n; own = false))
end

"""
    build_lineprofile(fits_path; a, θ, kwargs...)

Load `fits_path` and run `krz_tool_Lineprofile_run`. Returns `(; g, flux)`.
`θ` is observer inclination in degrees. `r_in=0` means ISCO.
`n_threads=0` uses every available core.
"""
function build_lineprofile(fits_path::AbstractString;
                           a::Float64 = 0.998,
                           θ::Float64 = 30.0,
                           M::Float64 = 1.0,
                           observer_r::Float64 = 1e7,
                           r_in::Float64 = 0.0,
                           r_out::Float64 = 400.0,
                           nradii::Integer = 100,
                           nangles::Integer = 200,
                           ng::Integer = 1000,
                           ngstar::Integer = 2800,
                           nrsteps::Integer = 3000,
                           normalise::Bool = true,
                           n_threads::Integer = 0)
    tool = KrzToolLineprofile(
        make_metric(M, a),
        KrzFourVector(0.0, observer_r, deg2rad(θ), 0.0),
        C_NULL,
        r_in,
        r_out,
        Csize_t(nradii),
        Csize_t(nangles),
        Csize_t(ng),
        Csize_t(ngstar),
        Csize_t(nrsteps),
        Cint(normalise),
        Cint(0),
    )

    pool = Ref(KrzThreadPool(C_NULL, 0))
    threadpool_init!(pool, n_threads)
    result = Ref(KrzLineprofile(C_NULL, C_NULL, 0))
    try
        GC.@preserve fits_path begin
            tool = KrzToolLineprofile(
                tool.metric, tool.x_obs,
                Base.unsafe_convert(Ptr{Cchar}, fits_path),
                tool.r_in, tool.r_out,
                tool.nradii, tool.nangles, tool.ng, tool.ngstar, tool.nrsteps,
                tool.normalise, tool._pad,
            )
            ret = ccall((:krz_tool_Lineprofile_run, libkerrz), Cint,
                        (Ptr{KrzThreadPool}, KrzToolLineprofile, Ptr{KrzLineprofile}),
                        pool, tool, result)
            _check_retcode(ret, "krz_tool_Lineprofile_run")
        end
        lp = result[]
        return (
            g = _copy_c_doubles(lp.g_grid, lp.num_g),
            flux = _copy_c_doubles(lp.flux, lp.num_g),
        )
    finally
        lineprofile_deinit!(result)
        threadpool_deinit!(pool)
    end
end

# ---------------------------------------------------------------------
# Grid of lamp-post line profiles
# ---------------------------------------------------------------------
spins = [0.01, 0.7, 0.998]
heights = [3, 8, 20]
thetas = [30, 50, 75]

@time begin
    for a in spins, h in heights
        ensure_emissivity(a, h)
    end

    panels = Plots.Plot[]
    for (i, h) in enumerate(heights)
        for (j, a) in enumerate(spins)
            p = plot(title = "a = $a, h = $h";
                     legend = (i == 1 && j == 1) ? :topright : false,
                     xlims = (0, 2))
            for θ in thetas
                fits = emissivity_path(a, h)
                profile = build_lineprofile(fits; a = Float64(a), θ = Float64(θ))
                plot!(p, profile.g, profile.flux; label = "θ = $(θ)°")
            end
            i == length(heights) && xlabel!(p, "energy shift g")
            j == 1 && ylabel!(p, "flux")
            push!(panels, p)
        end
    end

    plot(panels...;
         layout = (length(heights), length(spins)),
         size = (1400, 1000),
         plot_title = "Lamp-post line profiles")
end
