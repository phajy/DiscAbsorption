/// References
/// - Cunningham: Cunningham et al., 1975
const std = @import("std");
const zad = @import("zad");
const geometry = @import("geometry.zig");

const TEST_TOLERANCE = @import("options").test_numerical_tolerance;

const DualNumber = zad.DualNumber;
const FourVector = geometry.FourVector;
const KerrMetric = geometry.KerrMetric;

/// Represents Keplerian circular orbits in the Kerr metric.
pub fn Keplerian(comptime T: type) type {
    const A = T.Algebra;
    return struct {
        const Self = @This();

        // TODO: there is some room for optimisation here, caching e.g. the
        // angular velocity terms instead of recalculating them a couple of
        // times.

        sigma: T,
        delta: T,
        kerr_a: T,
        r: T,
        theta: T,
        a: T,
        M: T,

        pub fn init(metric: KerrMetric(T), r: T, theta: T) Self {
            const delta = metric.delta(r);
            const sigma = metric.sigma(r, theta);
            const kerr_a = metric.kerr_a(r, theta);

            return .{
                .delta = delta,
                .sigma = sigma,
                .kerr_a = kerr_a,
                .r = r,
                .theta = theta,
                .a = metric.a,
                .M = metric.M,
            };
        }

        /// Cunningham Equation (A2c)
        fn omega(self: Self) T {
            return A.div(
                A.mult(A.mult(.promote(2.0), self.a), A.mult(self.M, self.r)),
                self.kerr_a,
            );
        }

        /// The angular velocity of a Keplerian accreting gas around a central
        /// black hole.
        /// Cunningham Equation (A7b)
        pub fn angular_velocity(self: Self) T {
            const sqrt_M = A.sqrt(self.M);
            const r_to_the_three_halfs = A.mult(self.r, A.sqrt(self.r));
            return A.div(
                sqrt_M,
                A.add(r_to_the_three_halfs, A.mult(self.a, sqrt_M)),
            );
        }

        /// The angular frequency of a Keplerian accreting gas around a central
        /// black hole.
        pub fn angular_frequency(self: Self) T {
            const r_to_the_three_halfs = A.mult(self.r, A.sqrt(self.r));
            return A.div(
                .one,
                A.add(r_to_the_three_halfs, A.mult(self.a, self.M)),
            );
        }

        /// Angular velocity of the accreting gas in a locally non-rotating
        /// reference frame. For LNRF, see Bardeen et al. 1973. For this
        /// equation, see Cunningham Equation (A7b).
        ///
        /// Note: this equation is only valid above the ISCO.
        pub fn lnrf_velocity(self: Self) T {
            const numerator = A.mult(
                self.e_phi(),
                A.sub(self.angular_velocity(), self.omega()),
            );
            return A.div(numerator, self.e_nu());
        }

        /// Cunningham Equation (A2a)
        pub fn e_nu(self: Self) T {
            return A.sqrt(A.div(A.mult(self.sigma, self.delta), self.kerr_a));
        }

        /// Cunningham Equation (A2b)
        pub fn e_phi(self: Self) T {
            return A.mult(
                A.sin(self.theta),
                A.sqrt(A.div(self.kerr_a, self.sigma)),
            );
        }
    };
}

test "keplerian" {
    const Dual = DualNumber(f64, 0);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    const k = Keplerian(Dual).init(
        metric,
        .promote(2.0),
        // TODO: should this not be 2.0?
        .promote(std.math.pi / 4.0),
    );
    try std.testing.expectApproxEqAbs(0.42793345385870535, k.e_nu().x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(1.6490708105685496, k.e_phi().x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(0.16317825469123365, k.omega().x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(0.26134040121470514, k.angular_velocity().x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(0.3782745402934963, k.lnrf_velocity().x, TEST_TOLERANCE);
}

/// Represents plunging orbits in the Kerr metric below the ISCO.
///
/// This function is only valid at `r` below the ISCO.
pub fn plungingFourVelocity(
    comptime T: type,
    metric: KerrMetric(T),
    r: T,
) FourVector(T) {
    const ts = metric.tangentSpaceAlt(r, .promote(std.math.pi / 2.0));
    return plungingFourVelocityAlt(T, metric, ts);
}

/// Alternative signature of `plungingFourVelocity`
pub fn plungingFourVelocityAlt(
    comptime T: type,
    metric: KerrMetric(T),
    ts: KerrMetric(T).TangentSpace,
) FourVector(T) {
    const A = T.Algebra;
    std.debug.assert(ts.x.r.x <= metric.isco.x);

    // Following Mummery and Balbus (2022)
    // Equation (7) and (8) for the angular momentum and energy
    const _2M_3rI = A.div(
        A.mult(.promote(2), metric.M),
        A.mult(.promote(3), metric.isco),
    );

    // TODO: These are technically constants of the spacetime, and they could
    // be calculated once and re-used.
    const gamma = A.sqrt(A.sub(.one, _2M_3rI));
    const bracket = A.sub(
        .one,
        A.div(
            A.mult(.promote(2), metric.a),
            A.mult(.promote(3), A.sqrt(A.mult(metric.M, metric.isco))),
        ),
    );
    const _J = A.mult(
        A.mult(.promote(2), A.mult(.promote(@sqrt(3.0)), metric.M)),
        bracket,
    );

    // Equation (13)
    const _Ur = A.sqrt(
        A.mult(
            _2M_3rI,
            A.powi(A.sub(A.div(metric.isco, ts.x.r), .one), 3),
        ),
    ).neg();

    // Equation (14). Nominally the full formula would be used, but a tangent
    // space is needed in order to do the velocity constraint.
    const g_inv = ts.metric_components.inverse();
    const _Uphi = A.add(
        A.mult(g_inv.tph.neg(), gamma),
        A.mult(g_inv.phph, _J),
    );

    // The time-velocity component is solved by the velcotiy normalisation
    // self.constrainTime(v, magnitude)
    return ts.constrainVector(.{
        .t = .zero,
        .r = _Ur,
        .th = .zero,
        .ph = _Uphi,
    }, 1.0);
}

test "plunging four vectors" {
    const Dual = DualNumber(f64, 0);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.5));

    {
        const v = plungingFourVelocity(Dual, metric, .promote(2.2));
        // Compared against Gradus.jl's numerical scheme
        try std.testing.expectApproxEqAbs(5.161078761977851, v.t.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.35253558698653453, v.r.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0, v.th.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.9871246622581362, v.ph.x, TEST_TOLERANCE);
    }
}

/// Calculate the four-velocity of a cirular orbit in the equatorial plane at a
/// particular radius.
///
/// Note that this function does not include sub-ISCO velocities. For that, use
/// `keplerianPlunging` or `keplerianPlungingAlt` instead.
pub fn circularFourVelocity(comptime T: type, metric: KerrMetric(T), r: T) FourVector(T) {
    const A = T.Algebra;
    const k = Keplerian(T).init(metric, r, .promote(std.math.pi / 2.0));
    // Cunningham Equation (A7a)
    const prefactor = A.mult(
        k.e_nu(),
        A.sqrt(A.sub(.one, A.powi(k.lnrf_velocity(), 2))),
    );
    const common = A.div(.one, prefactor);
    return .{
        .t = common,
        .r = .zero,
        .th = .zero,
        .ph = A.mult(common, k.angular_velocity()),
    };
}

test "circular orbit four vectors" {
    const Dual = DualNumber(f64, 0);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));

    {
        const v = circularFourVelocity(Dual, metric, .promote(5.0));
        try std.testing.expectApproxEqAbs(1.4320923259190965, v.t.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.11759339443212102, v.ph.x, TEST_TOLERANCE);
        // check normalistation condition
        const ts = metric.tangentSpaceAlt(.promote(5.0), .promote(std.math.pi / 2.0));
        const v_norm = v.dot(ts, v);
        try std.testing.expectApproxEqAbs(-1.0, v_norm.x, TEST_TOLERANCE);
    }
}

/// A prescription for a co-rotating corona velocity. This velocity vector is
/// constrained to be a time-like four-vector.
pub fn coRotating(
    comptime T: type,
    metric: KerrMetric(T),
    ts: KerrMetric(T).TangentSpace,
) FourVector(T) {
    const A = T.Algebra;
    const sin_th = A.sin(ts.x.th);
    var projected_r = A.mult(sin_th, ts.x.r);

    projected_r = A.max(metric.isco, projected_r);

    const circ_v = circularFourVelocity(
        T,
        metric,
        projected_r,
    ).scalarMult(sin_th);
    const norm = A.sqrt(A.abs(circ_v.properNorm(ts)));
    const normed_v = circ_v.scalarMult(A.div(.one, norm));
    return ts.constrainVector(normed_v, 1.0);
}

test "co-rotating velocities" {
    const Dual = DualNumber(f64, 0);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));

    {
        const ts = metric.tangentSpaceAlt(
            .promote(5.0),
            .promote(std.math.degreesToRadians(30)),
        );
        const v = coRotating(Dual, metric, ts);
        try std.testing.expectApproxEqAbs(1.6124461173765703, v.t.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.32569095608917004, v.ph.x, TEST_TOLERANCE);
    }
}

/// A prescription for a velocity profile that follows Keplerian orbits of the
/// accretion disc with the plunging region prescription. This is valid all the
/// way down to the event horizon radius.
pub fn keplerianPlunging(
    comptime T: type,
    metric: KerrMetric(T),
    ts: KerrMetric(T).TangentSpace,
) FourVector(T) {
    const A = T.Algebra;
    const sin_th = A.sin(ts.x.th);
    var projected_r = A.mult(sin_th, ts.x.r);

    projected_r = A.max(metric.horizon_radius, projected_r);

    var circ_v = if (projected_r.x < metric.isco.x)
        plungingFourVelocity(T, metric, projected_r)
    else
        circularFourVelocity(
            T,
            metric,
            projected_r,
        );

    // Ensure valid outside of the equatorial plane.
    circ_v = circ_v.scalarMult(sin_th);

    const norm = A.sqrt(A.abs(circ_v.properNorm(ts)));
    const normed_v = circ_v.scalarMult(A.div(.one, norm));
    return ts.constrainVector(normed_v, 1.0);
}

pub fn keplerianPlungingAlt(comptime T: type, metric: KerrMetric(T), r: T) FourVector(T) {
    const ts = metric.tangentSpaceAlt(r, .promote(std.math.pi / 2.0));
    return keplerianPlunging(T, metric, ts);
}

/// Return a time-like velocity vector with all spatial components set to zero.
pub fn stationary(
    comptime T: type,
    ts: KerrMetric(T).TangentSpace,
) FourVector(T) {
    const A = T.Algebra;
    const m = ts.metric_components;
    return .{
        .t = A.div(.one, A.sqrt(m.tt.neg())),
        .r = .zero,
        .ph = .zero,
        .th = .zero,
    };
}

/// The locally non-rotating frame.
pub fn lnr(
    comptime T: type,
    ts: KerrMetric(T).TangentSpace,
) FourVector(T) {
    return ts.lnrVelocity();
}

/// Return a time-like velocity vector with only a radial component. Beta is
/// the radial velocity in units of `c`. Negative values are inward directed,
/// positive values are outwards.
pub fn radialMotion(
    comptime T: type,
    ts: KerrMetric(T).TangentSpace,
    beta: T,
) FourVector(T) {
    const A = T.Algebra;
    const m = ts.metric_components;
    return .{
        // This expression comes from normalising the velocity to -1.
        .t = A.sqrt(A.div(
            A.add(.one, A.mult(A.powi(beta, 2), m.rr)).neg(),
            m.tt,
        )),
        .r = beta,
        .ph = .zero,
        .th = .zero,
    };
}

test "radial motion" {
    const Dual = DualNumber(f64, 0);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    const ts = metric.tangentSpaceAlt(
        .promote(5.0),
        .promote(std.math.degreesToRadians(0.1)),
    );
    {
        const v = radialMotion(Dual, ts, .promote(0.7));
        try std.testing.expectApproxEqAbs(-1.0, v.dot(ts, v).x, TEST_TOLERANCE);
    }
    {
        const v = radialMotion(Dual, ts, .promote(0.0));
        try std.testing.expectApproxEqAbs(-1.0, v.dot(ts, v).x, TEST_TOLERANCE);
    }
}

/// The velocity profiles combined into a single interface.
pub const VelocityProfiles = enum {
    co_rotate,
    keplerian_plunging,
    lnr,
    stationary,

    pub fn fourVector(
        self: VelocityProfiles,
        comptime T: type,
        metric: KerrMetric(T),
        ts: KerrMetric(T).TangentSpace,
    ) FourVector(T) {
        return switch (self) {
            .co_rotate => coRotating(T, metric, ts),
            .keplerian_plunging => keplerianPlunging(T, metric, ts),
            .lnr => lnr(T, ts),
            .stationary => stationary(T, ts),
        };
    }
};

/// Determine the energy and angular momentum of the circular orbit at a
/// particular radius in the spacetime.
pub fn keplerianEnergyAngularMomentum(
    comptime T: type,
    metric: KerrMetric(T),
    r: T,
) struct {
    energy: T,
    ang_mom: T,
} {
    const ts = metric.tangentSpaceAlt(r, .promote(std.math.pi / 2.0));
    const vel = keplerianPlunging(T, metric, ts);
    const mom = ts.lowerIndices(vel);
    return .{
        .energy = mom.t.neg(),
        .ang_mom = mom.ph,
    };
}

test "keplerian energy and angular momentum" {
    const Dual = zad.DualNumber(f64, 0);
    const m = KerrMetric(Dual).init(.one, .promote(0.998));
    const eam = keplerianEnergyAngularMomentum(Dual, m, m.isco);

    try std.testing.expectApproxEqAbs(0.6790058343838, eam.energy.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(1.3918131634600912, eam.ang_mom.x, TEST_TOLERANCE);
}

/// Calculate the Lorentz factor relative to a LNR observer for a Keplerian
/// circular orbit.
pub fn lorentzFactorKeplerian(
    comptime T: type,
    metric: KerrMetric(T),
    ts: KerrMetric(T).TangentSpace,
) T {
    const A = T.Algebra;

    if (ts.x.r.x > metric.isco.x) {
        @branchHint(.likely);
        const keplerian = Keplerian(T).init(metric, ts.x.r, ts.x.th);
        const v_phi_local = keplerian.lnrf_velocity();
        const factor = A.div(.one, A.sqrt(
            A.sub(.one, A.powi(v_phi_local, 2)),
        ));
        return factor;
    }

    // This includes the `v_r` component in the local Lorentz factor.
    const vec_global = plungingFourVelocity(T, metric, ts.x.r);
    const vec_local = ts.lnrBasis().apply(vec_global);
    const v_magnitude_sq = A.add(
        A.powi(A.div(vec_local.r, vec_local.t), 2),
        A.powi(A.div(vec_local.ph, vec_local.t), 2),
    );

    return A.div(
        .one,
        A.sqrt(A.sub(.one, v_magnitude_sq)),
    );
}

/// Compute the (local) Lorentz factor at a particular point in the spacetime
/// given a velocity expressed in the global Boyer--Lindquist coordinates.
///
/// This function maps the velocity vector into the tangent space and then
/// calculates the magnitude for the Lorentz factor term.
pub fn lorentzFactor(
    comptime T: type,
    ts: KerrMetric(T).TangentSpace,
    v: FourVector(T),
) T {
    const A = T.Algebra;

    const vec_local = ts.lnrBasis().apply(v);
    const v_magnitude_sq = A.add(
        A.powi(A.div(vec_local.r, vec_local.t), 2),
        A.powi(A.div(vec_local.ph, vec_local.t), 2),
    );

    return A.div(
        .one,
        A.sqrt(A.sub(.one, v_magnitude_sq)),
    );
}

test "lorentz factor" {
    const Dual = DualNumber(f64, 1);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.5));
    {
        const ts = metric.tangentSpaceAlt(.promote(2.2), .promote(std.math.pi / 2.0));
        const lf = lorentzFactorKeplerian(Dual, metric, ts);
        try std.testing.expectApproxEqAbs(
            1.859176255086245,
            lf.x,
            TEST_TOLERANCE,
        );
    }
    {
        const ts = metric.tangentSpaceAlt(.promote(6.2), .promote(std.math.pi / 2.0));
        const lf = lorentzFactorKeplerian(Dual, metric, ts);
        try std.testing.expectApproxEqAbs(
            1.1154091379015822,
            lf.x,
            TEST_TOLERANCE,
        );
    }
}
