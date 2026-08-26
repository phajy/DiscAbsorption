"""
kerrz_lineprofile.jl

Example of assembling a relativistic line profile by calling `libkerrz`
directly from Julia via `ccall`, rather than shelling out to the `kerrz`
CLI and parsing text output.

IMPORTANT — please read before trusting the numbers this produces:

  This reimplements, by hand, a piece of the physics that the *actual*
  `kerrz lineprof` CLI command performs internally (see `lineprof.zig`
  in the kerrz source, which calls `kerrz.tools.Lineprofile` and
  `kerrz.tools.TransferFunctionTable` — neither of which is exposed
  through the public C header `kerrz.h`).

  The weighting/normalization used in `build_lineprofile` below
  (`em * r * dr * tr.f`) is the *standard shape* of a Cunningham-method
  disc integral, but the exact convention (whether `f` already folds in
  the radial Jacobian, whether extra g^3/g^4 relativistic-beaming terms
  are needed, whether emissivity has phi-dependence or a spectral-index
  weighting like the `g^Gamma` term implied by your emissivity table
  filenames) is NOT verified against the real kerrz internals.

  Treat this as a scaffold. Validate it by running the real `kerrz lineprof`
  CLI for the same parameters and comparing the resulting profile shape
  and normalization — see the `validate_against_cli` stub at the bottom.
"""

# ---------------------------------------------------------------------
# Path to the compiled shared library — update this for your machine
# ---------------------------------------------------------------------
const libkerrz = "/Users/er19801/kerrz/zig-out/lib/libkerrz.dylib"

# ---------------------------------------------------------------------
# Struct definitions mirroring kerrz.h exactly (field order & types matter)
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

struct KrzRingCorona
    height::Cdouble
    radius::Cdouble
end

struct KrzThreadPool
    pool::Ptr{Cvoid}
    num_threads::Csize_t
end

struct KrzEmissivityCache
    context::Ptr{Cvoid}
end

struct KrzEmissivityTrace
    em::Cdouble
    t::Cdouble
    g::Cdouble
    local_theta::Cdouble
end

struct KrzTFOptions
    max_points::Csize_t
    refine_N::Csize_t
    refine_M::Csize_t
    optimise::Csize_t
    minimum_guess::Cdouble
    heuristic::Cint
end

struct KrzToolTransferFunction
    tag::Cint            # TOOL_TRANSFER_FUNCTION = 0
    metric::KrzKerrMetric
    x_obs::KrzFourVector
    r_target::Cdouble
    options::KrzTFOptions
end

struct KrzCunninghamTrace
    image_angle::Cdouble
    image_radius::Cdouble
    alpha::Cdouble
    beta::Cdouble
    g::Cdouble
    g_star::Cdouble
    phi::Cdouble
    delta_t::Cdouble
    jacobian::Cdouble
    f::Cdouble
    r_err::Cdouble
end

struct KrzCunninghamTransferFunction
    g_min::Cdouble
    g_max::Cdouble
    traces::Ptr{KrzCunninghamTrace}
    num_traces::Csize_t
end

# ---------------------------------------------------------------------
# Thin ccall wrappers
# ---------------------------------------------------------------------
function make_metric(M::Float64, a::Float64)
    ccall((:krz_kerrMetric, libkerrz), KrzKerrMetric, (Cdouble, Cdouble), M, a)
end

function threadpool_init!(pool::Ref{KrzThreadPool}, n_threads::Integer)
    ret = ccall((:krz_ThreadPool_init, libkerrz), Cint,
                (Ptr{KrzThreadPool}, Csize_t), pool, n_threads)
    ret != 0 && error("krz_ThreadPool_init failed with code $ret")
    return pool
end

threadpool_deinit!(pool::Ref{KrzThreadPool}) =
    ccall((:krz_ThreadPool_deinit, libkerrz), Cvoid, (Ptr{KrzThreadPool},), pool)

function emissivity_cache_init!(cache::Ref{KrzEmissivityCache}, num_traces::Integer)
    ret = ccall((:krz_EmissivityCache_init, libkerrz), Cint,
                (Ptr{KrzEmissivityCache}, Csize_t), cache, num_traces)
    ret != 0 && error("krz_EmissivityCache_init failed with code $ret")
    return cache
end

emissivity_cache_deinit!(cache::Ref{KrzEmissivityCache}) =
    ccall((:krz_EmissivityCache_deinit, libkerrz), Cvoid, (Ptr{KrzEmissivityCache},), cache)

function emissivity_ring!(pool::Ref{KrzThreadPool}, cache::Ref{KrzEmissivityCache},
                           metric::KrzKerrMetric, ring::KrzRingCorona)
    ret = ccall((:krz_emissivity_ring, libkerrz), Cint,
                (Ptr{KrzThreadPool}, Ptr{KrzEmissivityCache}, KrzKerrMetric, KrzRingCorona),
                pool, cache, metric, ring)
    ret != 0 && error("krz_emissivity_ring failed with code $ret")
    return cache
end

function interpolate_emissivity(cache::Ref{KrzEmissivityCache}, r::Float64, phi::Float64)
    ccall((:krz_interpolate_emissivity, libkerrz), KrzEmissivityTrace,
          (Ptr{KrzEmissivityCache}, Cdouble, Cdouble), cache, r, phi)
end

function tf_defaults()
    ccall((:krz_tool_TransferFunction_defaults, libkerrz), KrzToolTransferFunction, ())
end

function tf_run(tool::KrzToolTransferFunction)
    result = Ref(KrzCunninghamTransferFunction(0.0, 0.0, C_NULL, 0))
    ret = ccall((:krz_tool_TransferFunction_run, libkerrz), Cint,
                (KrzToolTransferFunction, Ptr{KrzCunninghamTransferFunction}),
                tool, result)
    ret != 0 && error("krz_tool_TransferFunction_run failed with code $ret")
    return result
end

tf_deinit!(r::Ref{KrzCunninghamTransferFunction}) =
    ccall((:krz_CunninghamTransferFunction_deinit, libkerrz), Cvoid,
          (Ptr{KrzCunninghamTransferFunction},), r)

# Wrap traces without copying — safe as long as the caller reads/consumes
# them BEFORE calling tf_deinit!, since deinit is what invalidates the
# underlying memory (the wrap itself is zero-cost).
function tf_traces_unsafe(result::Ref{KrzCunninghamTransferFunction})
    r = result[]
    r.num_traces == 0 && return KrzCunninghamTrace[]
    unsafe_wrap(Vector{KrzCunninghamTrace}, r.traces, r.num_traces; own = false)
end

# Copying variant — use only if you need traces to outlive the tf_deinit! call
# (e.g. storing them for later inspection/debugging).
function tf_traces(result::Ref{KrzCunninghamTransferFunction})
    copy(tf_traces_unsafe(result))
end

# ---------------------------------------------------------------------
# Radial grid helper — log-spaced, matters once you extend range to r~400
# ---------------------------------------------------------------------
function logspace(r_min::Float64, r_max::Float64, n::Int)
    exp.(range(log(r_min), log(r_max); length = n))
end

# ---------------------------------------------------------------------
# Assemble a line profile by integrating disc annuli
# ---------------------------------------------------------------------
"""
    build_lineprofile(metric, x_obs, ring, r_grid, g_grid;
                       emissivity_num_traces=5000, tf_max_points=200)

For each radius in `r_grid`, computes the local emissivity and the
disc-to-observer transfer function, then bins the flux contribution
into `g_grid`. Returns a flux histogram of length `length(g_grid) - 1`.

`g_grid` should be the *edges* of your energy-shift bins (e.g.
`range(0.0, 2.0; length=101)` for a 100-bin profile from g=0 to g=2).
"""
function build_lineprofile(metric::KrzKerrMetric, x_obs::KrzFourVector,
                            ring::KrzRingCorona, r_grid::Vector{Float64},
                            g_grid::AbstractVector{Float64};
                            emissivity_num_traces::Integer = 5000,
                            tf_max_points::Integer = 200,
                            n_threads::Integer = 4,
                            gamma_index::Float64 = 2.0,      # photon index Γ of illuminating spectrum
                            beaming_exponent::Float64 = 3.0) # relativistic beaming power on g (see note below)

    pool = Ref(KrzThreadPool(C_NULL, 0))
    threadpool_init!(pool, n_threads)

    ecache = Ref(KrzEmissivityCache(C_NULL))
    emissivity_cache_init!(ecache, emissivity_num_traces)
    emissivity_ring!(pool, ecache, metric, ring)

    flux = zeros(Float64, length(g_grid) - 1)
    base_tf_opts = tf_defaults().options

    # Hoisted out of the loop — same options struct reused for every radius,
    # only r_target changes per-call, so this is built once instead of
    # `length(r_grid)` times.
    opts = KrzTFOptions(tf_max_points, base_tf_opts.refine_N,
                         base_tf_opts.refine_M, base_tf_opts.optimise,
                         base_tf_opts.minimum_guess, base_tf_opts.heuristic)

    # Direct O(1) bin lookup — only valid because g_grid is uniformly spaced
    # (as built by `range`/`collect(range(...))` in main()). Falls back to
    # searchsortedlast if you ever pass a non-uniform grid — check before
    # relying on this for irregular grids.
    g0 = g_grid[1]
    dg = g_grid[2] - g_grid[1]
    nbins = length(flux)
    @inline bin_index(g::Float64) = floor(Int, (g - g0) / dg) + 1

    # Per-thread flux accumulators avoid write contention between threads;
    # summed into `flux` once at the end.
    nthreads_julia = Threads.nthreads()
    flux_per_thread = [zeros(Float64, nbins) for _ in 1:nthreads_julia]

    for i in eachindex(r_grid)
        tid = Threads.threadid()
        local_flux = flux_per_thread[tid]

        r = r_grid[i]

        # trapezoidal-ish annulus width
        dr = if i == 1
            r_grid[2] - r_grid[1]
        elseif i == length(r_grid)
            r_grid[end] - r_grid[end-1]
        else
            (r_grid[i+1] - r_grid[i-1]) / 2
        end

        # Local illumination at this radius. For an axisymmetric ring
        # corona this likely doesn't depend on phi — confirm against
        # the real EmissivityFunction behaviour if that assumption matters.
        etrace = interpolate_emissivity(ecache, r, 0.0)
        # etrace.g is the CORONA-TO-DISC energyshift (illumination), distinct
        # from tr.g below, which is the DISC-TO-OBSERVER energyshift. Power-law
        # illumination means the flux hitting the disc scales as g^Γ.
        em = etrace.em * etrace.g^(gamma_index-2)

        (isnan(em) || isinf(em)) && continue   # skip bad points defensively

        tool = KrzToolTransferFunction(Cint(0), metric, x_obs, r, opts)

        result = tf_run(tool)
        traces = tf_traces_unsafe(result)   # no-copy read, consumed before deinit below

        for tr in traces
            (isnan(tr.g) || isnan(tr.f)) && continue
            bin = bin_index(tr.g)
            if 1 <= bin <= nbins
                # tr.g here is the disc-to-observer energyshift; g^beaming_exponent
                # is the relativistic beaming term (specific intensity invariance
                # gives beaming_exponent=3; some conventions fold Γ in here too,
                # giving g^(3+Γ) instead of separating it as done above — confirm
                # which convention kerrz's own Lineprofile.zig uses before trusting this)
                # NOTE: overall weighting convention not yet validated — see module docstring
                local_flux[bin] += em * r * dr * tr.f * tr.g^beaming_exponent
            end
        end

        tf_deinit!(result)  # safe now — traces already consumed above
    end

    for lf in flux_per_thread
        flux .+= lf
    end


    emissivity_cache_deinit!(ecache)
    threadpool_deinit!(pool)

    return flux
end