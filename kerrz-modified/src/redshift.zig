const std = @import("std");
const geometry = @import("geometry.zig");
const geodesic = @import("geodesic.zig");
const orbits = @import("orbits.zig");

const DualNumber = @import("zad").DualNumber;
const KerrMetric = geometry.KerrMetric;
const FourVector = geometry.FourVector;
const Keplerian = orbits.Keplerian;
const NullGeodesic = geodesic.NullGeodesic;
const TraceResult = geodesic.TraceResult;

const TEST_TOLERANCE = @import("options").test_numerical_tolerance;

/// Calculate the redshift from a Keplerian accretion disc at radius `r` and
/// poloidal angle `theta`. Here, `lambda` is the (energy normalised) angular
/// momentum of the null-geodesic.
///
/// This calculation assumes that the origin of the geodesic is stationary at
/// infinity.
pub fn keplerianRedshiftAsymptotic(
    comptime T: type,
    metric: KerrMetric(T),
    geod: NullGeodesic(T),
    result: TraceResult(T),
) T {
    const v_src = orbits.stationary(T, metric.tangentSpace(geod.x_init));
    return keplerianRedshiftResult(T, metric, geod, result, v_src);
}

/// Calculate the energyshift with respect to a Keplerian accretion disc, where
/// the observer is located at `x_obs` and has velocity `v_obs`.
///
/// When the point is below the ISCO, uses the plunging velocity components for
/// a free-fall.
pub fn keplerianRedshift(
    comptime T: type,
    metric: KerrMetric(T),
    r: T,
    theta: T,
    x_obs: FourVector(T),
    v_obs: FourVector(T),
    k_obs: FourVector(T),
    k_disc: FourVector(T),
) T {
    const ts_obs = metric.tangentSpace(x_obs);
    const ts_disc = metric.tangentSpaceAlt(r, theta);

    const v_disc = if (r.x < metric.isco.x) b: {
        // TODO: make this more explicit somehow
        // Basically because of the time-reversal nature of the ray-tracing,
        // the r-component of the four velocity is reversed.
        var vec = orbits.plungingFourVelocity(T, metric, r);
        vec.r = vec.r.neg();
        break :b vec;
    } else orbits.circularFourVelocity(T, metric, r);
    return redshift(T, ts_obs, v_obs, k_obs, ts_disc, v_disc, k_disc);
}

/// Calculate the energyshift along a geodesic, as seen by an observer at
/// `ts_obs` with velocity `v_obs`, with the photon geodesic originating at
/// `ts_em` with medium velocity `v_em`.
///
/// See also `keplerianRedshift`. Alternative function interface that uses the
/// tangent spaces and observer / disc velocities directly.
pub fn redshift(
    comptime T: type,
    ts_obs: KerrMetric(T).TangentSpace,
    v_obs: FourVector(T),
    k_obs: FourVector(T),
    ts_em: KerrMetric(T).TangentSpace,
    v_em: FourVector(T),
    k_em: FourVector(T),
) T {
    const A = T.Algebra;
    const _E_obs = k_obs.dot(ts_obs, v_obs);
    const _E_em = k_em.dot(ts_em, v_em);
    return A.abs(A.div(_E_obs, _E_em));
}

/// Alternative signature for `redshift`.
pub fn redshiftFromResult(
    comptime T: type,
    metric: KerrMetric(T),
    geod: NullGeodesic(T),
    result: TraceResult(T),
    v_medium: FourVector(T),
) T {
    const k_disc = result.velocity(metric, geod);
    const k_obs = geod.initialVelocity(metric);
    const ts_obs = metric.tangentSpace(geod.x_init);
    const ts_disc = metric.tangentSpaceAlt(result.r, result.theta);
    return redshift(
        T,
        ts_obs,
        orbits.stationary(T, ts_obs),
        k_obs,
        ts_disc,
        v_medium,
        k_disc,
    );
}

/// See `keplerianRedshift`. This is an alternative interface to the same
/// function.
pub fn keplerianRedshiftResult(
    comptime T: type,
    metric: KerrMetric(T),
    geod: NullGeodesic(T),
    result: TraceResult(T),
    v_src: FourVector(T),
) T {
    const k_disc = result.velocity(metric, geod);
    const k_obs = geod.initialVelocity(metric);
    return keplerianRedshift(
        T,
        metric,
        result.r,
        result.theta,
        geod.x_init,
        v_src,
        k_obs,
        k_disc,
    );
}
