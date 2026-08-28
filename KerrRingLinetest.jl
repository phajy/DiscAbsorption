using Plots

# Ring-corona line profile via the file-based kerrz C API
# (`krz_tool_Lineprofile_run`). Loads a precomputed emissivity FITS file,
# builds transfer functions, and returns (g, flux) in memory.
#
# Adapted from the lamppost version of this script — the ccall wrappers and
# krz_tool_Lineprofile_run call are corona-agnostic (they just read whichever
# FITS file `emissivity_path` points to), so only the emissivity-generation
# step needed to change.
#
# Emissivity files: kerrz-modified/emissivity/ring/emis_ring_a{spin}_h{height}_r{radius}[...].fits
# Missing files are generated with the kerrz CLI and cached indefinitely —
# rerunning with the same (spin, height, radius, photon_index, velocity,
# nphotons) reuses the existing file instead of regenerating it, so the
# library of precomputed profiles grows over time rather than being rebuilt
# on every run.
# Do not use a=0 — kerrz panics in the elliptic integrals; use a small positive spin.
#
# VERIFY BEFORE RUNNING: the `--ring h:...,r:...` sub-key syntax and the
# top-level `--velocity` flag below are inferred by analogy with the
# lamppost command's `--lamppost h:...,vr:0` pattern — I don't have the CLI
# arg-parser source for the `emissivity` command to confirm the exact ring
# key names. Run `kerrz emissivity --help` once and adjust the single `run(...)`
# line in `ensure_emissivity_ring` below if the real syntax differs.

const kerrz_root = joinpath(@__DIR__, "kerrz-modified")
const libkerrz = joinpath(kerrz_root, "zig-out", "lib", "libkerrz.dylib")
const kerrz_cli = joinpath(kerrz_root, "zig-out", "bin", "kerrz")
const emissivity_dir = joinpath(kerrz_root, "emissivity")
const ring_emissivity_dir = joinpath(emissivity_dir, "ring")

# ---------------------------------------------------------------------
# Structs mirroring kerrz-modified/wrappers/kerrz.h
# (unchanged from the lamppost script — krz_tool_Lineprofile_run only ever
# sees a FITS path, it doesn't know or care which corona model produced it)
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
# Path helper and on-demand emissivity generation — ring corona
# ---------------------------------------------------------------------
_num_str(x) = isinteger(x) ? string(Int(x)) : string(x)

"""
    ring_emissivity_path(a, h, r; photon_index=2.0, velocity=:co_rotate, nphotons=3000)

Cache path for a ring corona's emissivity FITS file. Every parameter that
actually changes the FITS content is baked into the filename — spin, height,
and radius always; photon_index/velocity/nphotons only when they differ from
the defaults, so the common case keeps a short, readable filename while
still avoiding collisions if you later sweep those too.
"""
function ring_emissivity_path(a, h, r; photon_index::Float64 = 2.0,
                              velocity::Symbol = :co_rotate,
                              nphotons::Integer = 3000)
    gtag = photon_index == 2.0 ? "" : "_g$(_num_str(photon_index))"
    vtag = velocity === :co_rotate ? "" : "_v$(velocity)"
    ntag = nphotons == 3000 ? "" : "_n$(nphotons)"
    joinpath(ring_emissivity_dir,
             "emis_ring_a$(a)_h$(_num_str(h))_r$(_num_str(r))$(gtag)$(vtag)$(ntag).fits")
end

"""
    ensure_emissivity_ring(a, h, r; kwargs...) -> path

Return the FITS path for a ring corona at spin `a`, height `h`, radius `r`,
running `kerrz emissivity --ring ...` if the file is missing. Reuses the
existing file on every subsequent call with the same parameters — this is
what builds the on-disk library up over time instead of regenerating on
every run.
"""
function ensure_emissivity_ring(a, h, r; photon_index::Float64 = 2.0,
                                velocity::Symbol = :co_rotate,
                                nphotons::Integer = 3000)
    path = ring_emissivity_path(a, h, r; photon_index, velocity, nphotons)
    isfile(path) && return path
    isfile(kerrz_cli) || error(
        "kerrz CLI not found at $kerrz_cli\nBuild with: (cd kerrz-modified && zig build)")
    mkpath(ring_emissivity_dir)
    println("Generating ring emissivity a=$a h=$h r=$r v=$velocity -> $path")

    velocity_str = velocity === :co_rotate ? "corotate" : string(velocity)

    # SEE FILE-TOP NOTE: --ring h:,r: sub-keys and the top-level --velocity
    # flag are inferred, not confirmed — check `kerrz emissivity --help`.
    run(`$kerrz_cli emissivity --spin $a --ring-like h:$(_num_str(h)),x:$(_num_str(r))
         --velocity $velocity_str --photon-index $photon_index --nphotons $nphotons
         -o $path`)

    isfile(path) || error("kerrz emissivity did not write $path")
    return path
end

# ---------------------------------------------------------------------
# ccall wrappers (unchanged from the lamppost script — corona-agnostic)
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
# Grid of ring-corona line profiles
# ---------------------------------------------------------------------
# A ring corona is a genuine 2-parameter shape (height AND radius), unlike
# lamppost's single height — so the "size" axis below pairs the two rather
# than sweeping either alone. Adjust these pairings to whatever geometries
# you actually care about.
spins = [0.01, 0.7, 0.998]
ring_geometries = [(height = 3.0, radius = 1.5),
                    (height = 8.0, radius = 3.0),
                    (height = 20.0, radius = 5.0)]
thetas = [5, 30, 60, 85]

for a in spins, geo in ring_geometries
    ensure_emissivity_ring(a, geo.height, geo.radius)
end

panels = Plots.Plot[]
for (i, geo) in enumerate(ring_geometries)
    for (j, a) in enumerate(spins)
        p = plot(title = "a = $a, h = $(geo.height), r = $(geo.radius)";
                 legend = (i == 1 && j == 1) ? :topright : false,
                 xlims = (0, 2))
        for θ in thetas
            fits = ring_emissivity_path(a, geo.height, geo.radius)
            profile = @time build_lineprofile(fits; a = Float64(a), θ = Float64(θ))
            plot!(p, profile.g, profile.flux; label = "θ = $(θ)°")
        end
        i == length(ring_geometries) && xlabel!(p, "energy shift g")
        j == 1 && ylabel!(p, "flux")
        push!(panels, p)
    end
end

plot(panels...;
     layout = (length(ring_geometries), length(spins)),
     size = (1400, 1000),
     plot_title = "Ring-corona line profiles")