! vim: cc=80 wrap tw=80

! This module defines the types that the kerrz Fortran library has access to.
! It is recommended not to directly use this module but instead the
! `kerrz_interface` module, which exposes all of the same types.
module kerrz_derived_types
    use iso_c_binding, only: c_double, c_ptr, c_null_ptr, c_size_t, c_int
    implicit none

    ! Error codes that may be returned by some function calls.
    integer, parameter :: KRZ_RET_SUCCESS = 0
    integer, parameter :: KRZ_RET_ALLOCATION_FAILED = 1
    integer, parameter :: KRZ_RET_THREAD_ERROR = 2
    integer, parameter :: KRZ_RET_CONVERGENCE_FAILED = 3

    ! These are the status codes that indicate the fate of a null-geodesic.
    integer, parameter :: KRZ_STATUS_NONE = 0
    integer, parameter :: KRZ_STATUS_EVENT_HORIZON = 1
    integer, parameter :: KRZ_STATUS_INTERSECTED_DISC = 2
    integer, parameter :: KRZ_STATUS_AT_INFINITY = 3

    ! The kerrz thread pool. Construct with `krz_ThreadPool_init`, and free
    ! associated resources with `krz_ThreadPool_deinit`.
    type, bind(C) :: krz_ThreadPool
        type(c_ptr) :: pool
        integer(c_size_t) :: num_threads
    end type krz_ThreadPool

    ! The Kerr metric. Construct with `krz_KerrMetric_init`. This does not need
    ! to be freed, and lives entirely on the stack.
    type, bind(C) :: krz_KerrMetric
        real(c_double) :: M
        real(c_double) :: a
        real(c_double) :: horizon_radius, horizon_radius_negative
        real(c_double) :: isco
    end type krz_KerrMetric

    ! A ring corona abstract with some height and radius in the Boyer-Lindquist
    ! coordinates. It may be constructed directly
    ! ```fortran
    ! type(krz_RingCorona) :: ring
    ! ring = krz_RingCorona(4.0d.0, 2.0d0)
    ! ```
    type, bind(C) :: krz_RingCorona
        real(c_double) :: height
        real(c_double) :: radius
    end type krz_RingCorona

    ! The cache used to evaluate emissivity profiles, allowing resources to be
    ! reused between successive calls. Use `krz_emissivity_*` functions to
    ! populate the cache, and `krz_interpolate_emissivity` to obtain the
    ! emissivity as a function of radius.
    !
    ! Construct with `krz_EmissivityCache_init` and free with
    ! `krz_EmissivityCache_deinit`.
    type, bind(C) :: krz_EmissivityCache
        type(c_ptr) :: context
        integer(c_size_t) :: num_bins
        type(c_ptr) :: radii, emissivity
    end type krz_EmissivityCache

    ! A four-vector in the Boyer--Lindquist coordinates.
    type, bind(C) :: krz_FourVector
        real(c_double) :: t
        real(c_double) :: r
        real(c_double) :: th
        real(c_double) :: ph
    end type krz_FourVector

    ! Initial conditions for a null-geodesic.
    type, bind(C) :: krz_InitialConditions
        ! The geodesic energy.
        real(c_double) :: E
        ! The geodesic angular momentum.
        real(c_double) :: L
        ! The geodesic Carter constant.
        real(c_double) :: Q
        ! The reduced geodesic angular momentum L / E.
        real(c_double) :: lambda
        ! The reduced geodesic Carter constant Q / E^2.
        real(c_double) :: eta
        ! The initial four position.
        type(krz_FourVector) :: x_init
        ! The signs of the momenta and the (maximum) number of windings the
        ! geodesic should follow.
        real(c_double) :: theta_sign, radial_sign, windings
    end type krz_InitialConditions

    ! The result from a trace.
    type, bind(C) :: krz_TraceResult
        ! The status of the geodesic. Use `KRZ_STATUS_*` in comparisons.
        integer :: status
        ! The mino time of the geodesic.
        real(c_double) :: mino_time
        ! The end four-position of the geodesic.
        type(krz_FourVector) :: x_final
        real(c_double) :: winding
    end type krz_TraceResult

    ! Used to store a single evaluation of a (Cunningham) transfer function
    ! table.
    type, bind(C) :: krz_TransferFunctionCache
        type(c_ptr) :: cache
    end type krz_TransferFunctionCache

    ! A continuum transfer function.
    type, bind(C) :: krz_ContinuumRing
        type(c_ptr) :: cache
    end type krz_ContinuumRing

    ! A point on the observer plane that is interpolated from a (Cunningham)
    ! transfer function.
    type, bind(C) :: krz_InterpolatedPoint
        real(c_double) :: alpha, beta
        real(c_double) :: delta_t
    end type krz_InterpolatedPoint

    ! A point along the continuum transfer function for a ring-like corona.
    type, bind(C) :: krz_ContinuumRingPoint
        real(c_double) :: delta_t
        real(c_double) :: energyshift
        real(c_double) :: dcosd_dcosth
    end type krz_ContinuumRingPoint

    ! An orthonormal frame at position `x`. This is returned from one of the
    ! frame constructors, e.g. `krz_stationaryFrame`.
    type, bind(C) :: krz_OrthonormalFrame
        type(krz_FourVector) :: x
        real(c_double) :: metric_components(5)
        real(c_double) :: matrix(16)
    end type krz_OrthonormalFrame

    ! The trace solution for a continuum spectrum, in this case from the corona
    ! to the observer.
    type, bind(C) :: krz_ContinuumLamppost
        ! The corona-to-observer photon.
        type(krz_TraceResult) :: res
        ! The local sky angle, defined off of the zenith (i.e. the spin axis).
        real(c_double) :: angle_delta
        ! The Jacobian term.
        real(c_double) :: dcosd_dcosth
        ! The impact parameters.
        real(c_double) :: alpha, beta
    end type krz_ContinuumLamppost

    ! The result from interpolating an emissivity profile on the surface of the
    ! accretion disc.
    type, bind(C) :: krz_EmissivityTrace
        ! The emissivity itself.
        real(c_double) :: em
        ! The corona-to-disc light travel time.
        real(c_double) :: t
        ! The corona-to-disc energyshift.
        real(c_double) :: g
        ! The incident angle between the photon and the disc surface.
        ! A value of `0` is directly parallel with the disc normal at that
        ! radius.
        real(c_double) :: local_theta
    end type krz_EmissivityTrace

end

! This module contains all definitions, bindings, and wrapper of the kerrz
! library for Fortran 90.
!
! It is recommended to define your own `kerrz` module which uses
! `kerrz_interface`, and stores any global state and writes application-specific
! wrapper routines. That way, in your code you can make use of the semantically
! helpful `use kerrz` elsewhere and the intention and meanings be obvious.
module kerrz_interface
    use kerrz_derived_types
    implicit none

    interface

    ! Initialise a thread pool. Returns a krz_RETCODE.
    integer(c_int) function krz_ThreadPool_init(pool, num_threads)             &
        bind(C, name="krz_ThreadPool_init") result(status)
        use kerrz_derived_types, only: krz_ThreadPool, c_size_t, c_int
        type(krz_ThreadPool), intent(inout) :: pool
        integer(c_size_t), value, intent(in) :: num_threads
    end function krz_ThreadPool_init

    ! Free the thread pool.
    subroutine krz_ThreadPool_deinit(pool) bind(C, name="krz_ThreadPool_deinit")
        use kerrz_derived_types, only: krz_ThreadPool
        type(krz_ThreadPool), intent(inout) :: pool
    end subroutine krz_ThreadPool_deinit

    ! Initialise the metric and calculate fixed quantities of the spacetime.
    ! This has a different name from the C wrapper as Fortran is case
    ! insensitive.
    type(krz_KerrMetric) function krz_KerrMetric_init(M, a)                    &
        bind(C, name="krz_kerrMetric") result(met)
        use kerrz_derived_types, only: krz_KerrMetric, c_double
        real(c_double), value, intent(in) :: M, a
    end function krz_KerrMetric_init

    ! Intiialise the emissivity cache. Returns a krz_RETCODE.
    integer(c_int) function krz_EmissivityCache_init(cache, num_traces)        &
        bind(C, name="krz_EmissivityCache_init") result(status)
        use kerrz_derived_types, only: krz_EmissivityCache, c_int, c_size_t
        type(krz_EmissivityCache), intent(inout) :: cache
        integer(c_size_t), value, intent(in) :: num_traces
    end function krz_EmissivityCache_init

    ! Free the emissivity cache.
    subroutine krz_EmissivityCache_deinit(cache)                               &
        bind(C, name="krz_EmissivityCache_deinit")
        use kerrz_derived_types, only: krz_EmissivityCache
        type(krz_EmissivityCache), intent(inout) :: cache
    end subroutine krz_EmissivityCache_deinit

    ! Calculate an emissivity profile for the ring corona, making use of the
    ! thread pool and pre-allocated cache.
    integer(c_int) function krz_emissivity_ring(pool, cache, metric, ring)     &
        bind(C, name="krz_emissivity_ring") result(status)
        use kerrz_derived_types, only: c_int, krz_ThreadPool,                  &
            krz_EmissivityCache, krz_RingCorona, krz_KerrMetric
        type(krz_ThreadPool), intent(inout) :: pool
        type(krz_EmissivityCache), intent(inout) :: cache
        type(krz_KerrMetric), value, intent(in) :: metric
        type(krz_RingCorona), value, intent(in) :: ring
    end function krz_emissivity_ring

    ! Intepolate the (time-averaged) emissivity at a partciular radius and
    ! azimuth on the accretion disc from a cache.
    type(krz_EmissivityTrace) function krz_interpolate_emissivity(cache, r,    &
        phi) bind(C, name="krz_interpolate_emissivity")
        use kerrz_derived_types, only: c_double, krz_EmissivityCache,          &
            krz_EmissivityTrace
        type(krz_EmissivityCache), intent(inout) :: cache
        real(c_double), value, intent(in) :: r, phi
    end function krz_interpolate_emissivity

    ! Intepolate the continuum transfer function at a particular azimuthal point
    ! along the ring corona.
    type(krz_ContinuumRingPoint) function krz_interpolate_continuum(cache,     &
        phi) bind(C, name="krz_interpolate_continuum")
        use kerrz_derived_types, only: c_double, krz_ContinuumRing,            &
            krz_ContinuumRingPoint
        type(krz_ContinuumRing), intent(inout) :: cache
        real(c_double), value, intent(in) :: phi
    end function krz_interpolate_continuum

    ! Calculate (Cunningham) transfer functions for a particular metric and
    ! observer position between some inner and outer radius on the accretion
    ! disc.  Calculates `num_radii` worth of transfer functions. By default,
    ! these are linearly spaced between `r_min` and `r_max`.
    integer(c_int) function krz_transfer_functions(cache, metric, observer,    &
        r_min, r_max, num_radii)                                               &
        bind(C, name="krz_transfer_functions") result(status)
        use kerrz_derived_types, only: c_int, c_double,                        &
            krz_TransferFunctionCache,                                         &
            krz_KerrMetric, krz_TransferFunctionCache, krz_FourVector,         &
            c_size_t
        type(krz_TransferFunctionCache), intent(inout) :: cache
        type(krz_KerrMetric), value, intent(in) :: metric
        type(krz_FourVector), value, intent(in) :: observer
        real(c_double), value, intent(in) :: r_min, r_max
        integer(c_size_t), value, intent(in) :: num_radii
    end function krz_transfer_functions

    ! Interpolate a point on the image plane given a set of coordinates on the
    ! accretion disc.
    type(krz_InterpolatedPoint) function krz_interpolate_disc_coordinates(     &
            cache, r, phi)                                                     &
        bind(C, name="krz_interpolate_disc_coordinates") result(pt)
        use kerrz_derived_types, only: krz_InterpolatedPoint,                  &
            krz_TransferFunctionCache, c_double
        type(krz_TransferFunctionCache), intent(inout) :: cache
        real(c_double), value, intent(in) :: r, phi
    end function krz_interpolate_disc_coordinates

    ! Interpolate a point on the image plane given a set of coordinates on the
    ! accretion disc.
    integer(c_int) function krz_shadow(metric, obs_incl, alpha, beta, num_pts) &
        bind(C, name="krz_shadow") result(status)
        use kerrz_derived_types, only: c_int, krz_KerrMetric,                  &
            c_double, c_size_t
        type(krz_KerrMetric), value, intent(in) :: metric
        real(c_double), value, intent(in) :: obs_incl
        integer(c_size_t), value, intent(in) :: num_pts
        real(c_double), intent(inout) :: alpha(num_pts), beta(num_pts)
    end function krz_shadow

    ! Calculate the initial conditions for a null-geodesic from impact
    ! parameters defined on the image plane of a (distant) observer.
    !
    ! Note: this currently uses formulae for mapping the impact parameters to
    ! the initial momentum that are onlid valid in the asymptotically flat
    ! limit, and so will be erroneous if `x_init%r` is small (i.e. <100 rg).
    type(krz_InitialConditions) function krz_fromImpactParameters(             &
        metric, x_init, alpha, beta) bind(C, name="krz_fromImpactParameters")  &
        result (cond)
        use kerrz_derived_types, only:  krz_KerrMetric, c_double,              &
            krz_FourVector, krz_InitialConditions
        type(krz_KerrMetric), value, intent(in) :: metric
        type(krz_FourVector), value, intent(in) :: x_init
        real(c_double), value, intent(in) :: alpha, beta
    end function krz_fromImpactParameters

    ! Trace the geodesic to a particular polar angle.
    type(krz_TraceResult) function krz_traceToAngle(metric, ic, angle)         &
        bind(C, name="krz_traceToAngle") result (res)
        use kerrz_derived_types, only:  krz_KerrMetric, krz_InitialConditions, &
            krz_TraceResult, c_double
        type(krz_KerrMetric), value, intent(in) :: metric
        type(krz_InitialConditions), value, intent(in) :: ic
        real(c_double), value, intent(in) :: angle
    end function krz_traceToAngle

    ! Trace the geodesic to a particular radius.
    type(krz_TraceResult) function krz_traceToRadius(metric, ic, radius)       &
        bind(C, name="krz_traceToRadius") result (res)
        use kerrz_derived_types, only:  krz_KerrMetric, krz_InitialConditions, &
            krz_TraceResult, c_double
        type(krz_KerrMetric), value, intent(in) :: metric
        type(krz_InitialConditions), value, intent(in) :: ic
        real(c_double), value, intent(in) :: radius
    end function krz_traceToRadius

    ! Calculate a stationary orthonormal frame.
    type(krz_OrthonormalFrame) function krz_stationaryFrame(metric, x)         &
        bind(C, name="krz_stationaryFrame") result (frame)
        use kerrz_derived_types, only:  krz_KerrMetric, krz_FourVector,        &
            krz_OrthonormalFrame
        type(krz_KerrMetric), value, intent(in) :: metric
        type(krz_FourVector), value, intent(in) :: x
    end function krz_stationaryFrame

    ! Map local angles on the sky (theta = 0 points directly up, theta = pi
    ! directly down) in the orthonormal frame to initial conditions for tracing
    ! a geodesic.
    type(krz_InitialConditions) function krz_fromSkyAngles(metric, frame,      &
        theta, phi) bind(C, name="krz_fromSkyAngles") result (cond)
        use kerrz_derived_types, only:  krz_KerrMetric, krz_InitialConditions, &
            krz_OrthonormalFrame, c_double
        type(krz_KerrMetric), value, intent(in) :: metric
        type(krz_OrthonormalFrame), value, intent(in) :: frame
        real(c_double), value, intent(in) :: theta, phi
    end function krz_fromSkyAngles

    ! Trace a continuum-to-observer photon. This will use a root-solving method
    ! to find the local sky angle that maps back to the observer, and returning
    ! the angle and lensing factor in one go.
    type(krz_ContinuumLamppost) function krz_traceContinuumLamppost(metric,    &
        x_obs, height) bind(C, name="krz_traceContinuumLamppost") result (cond)
        use kerrz_derived_types, only:  krz_KerrMetric, krz_ContinuumLamppost, &
            c_double, krz_FourVector
        type(krz_KerrMetric), value, intent(in) :: metric
        type(krz_FourVector), value, intent(in) :: x_obs
        real(c_double), value, intent(in) :: height
    end function krz_traceContinuumLamppost

    ! Free resources allocated for the continuum transfer function.
    subroutine krz_ContinuumRing_deinit(cache)                                 &
        bind(C, name="krz_ContinuumRing_deinit")
        use kerrz_derived_types, only: krz_ContinuumRing
        type(krz_ContinuumRing), intent(inout) :: cache
    end subroutine krz_ContinuumRing_deinit

    ! Trace a continuum-to-observer transfer function for a ring-like corona.
    ! The resulting cache can be interpolated along the azimuthal coordinate of
    ! the ring with `krz_interpolate_continuum`. Returns a `krz_RETCODE`.
    integer(c_int) function krz_traceContinuumRing(cache, metric,              &
        x_obs, ring) bind(C, name="krz_traceContinuumRing") result (cond)
        use kerrz_derived_types, only:  krz_KerrMetric, krz_ContinuumRing,     &
            krz_FourVector, krz_RingCorona, c_int
        type(krz_ContinuumRing), intent(out) :: cache
        type(krz_KerrMetric), value, intent(in) :: metric
        type(krz_FourVector), value, intent(in) :: x_obs
        type(krz_RingCorona), value, intent(in) :: ring
    end function krz_traceContinuumRing

    end interface

contains

    double precision function mod2pi(angle) result(result_angle)
        double precision :: angle
        double precision, parameter :: pi = acos(-1.0)
        result_angle = mod(angle, 2.0 * pi)
    end function mod2pi

    double precision function deg2rad(angle) result(result_angle)
        double precision :: angle
        double precision, parameter :: pi = acos(-1.0)
        result_angle = angle * pi / 180
    end function deg2rad

end module kerrz_interface
