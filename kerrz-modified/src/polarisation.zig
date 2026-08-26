const std = @import("std");
const ad = @import("zad");
const geometry = @import("geometry.zig");
const geodesic = @import("geodesic.zig");
const orbits = @import("orbits.zig");

const options = @import("options");

const TEST_TOLERANCE = options.test_numerical_tolerance;

const FourVector = geometry.FourVector;
const KerrMetric = geometry.KerrMetric;
const ComplexNumber = @import("complex.zig").ComplexNumber;

/// Computes a.x * b.y - a.y * b.x
fn commutator(comptime T: type, ax: T, by: T, ay: T, bx: T) T {
    const A = T.Algebra;
    return A.sub(A.mult(ax, by), A.mult(ay, bx));
}

test commutator {
    const Dual = ad.DualNumber(f64, 0);
    const result = commutator(
        Dual,
        .promote(3.0),
        .promote(-2.0),
        .promote(7.0),
        .promote(3.4),
    );
    try std.testing.expectApproxEqAbs(
        -29.8,
        result.x,
        TEST_TOLERANCE,
    );
}

/// The normal vector in the global coordinates of a patch of the equatorial
/// thin disc.
///
/// Calculates the time component consistently.
pub fn normalVectorEquatorial(
    comptime T: type,
    ts: KerrMetric(T).TangentSpace,
    v_medium: FourVector(T),
) FourVector(T) {
    const A = T.Algebra;
    std.debug.assert(
        std.math.approxEqAbs(T.T, ts.x.th.x, std.math.pi / 2.0, TEST_TOLERANCE),
    );

    const numerator = A.mult(ts.metric_components.thth, v_medium.th).neg();
    const denom = A.add(
        A.mult(ts.metric_components.tt, v_medium.t),
        A.mult(ts.metric_components.tph, v_medium.ph),
    );

    const t_component = A.div(numerator, denom);

    return .{
        .t = t_component,
        .r = .zero,
        // The e^mu_(theta) component of the local tetrad. Since the disc patch
        // is moving only in time and at most r and phi, this can be written
        // simply as:
        .th = .promote(1.0 / @sqrt(ts.metric_components.thth.x)),
        .ph = .zero,
    };
}

test "normalVectorEquatorial" {
    const Dual = ad.DualNumber(f64, 1);
    const m = KerrMetric(Dual).init(.one, .promote(0.998));
    const x = FourVector(Dual){
        .t = .zero,
        .ph = .zero,
        .r = .promote(10.0),
        .th = .promote(std.math.pi / 2.0),
    };
    const ts = m.tangentSpace(x);

    const v_medium = orbits.keplerianPlunging(Dual, m, ts);
    const n = normalVectorEquatorial(Dual, ts, v_medium);

    const magnitude = n.dot(ts, n);

    try std.testing.expectApproxEqAbs(
        magnitude.x,
        1,
        TEST_TOLERANCE,
    );
}

pub fn PolarisationConstant(comptime T: type) type {
    const A = T.Algebra;
    return struct {
        pub const XY = struct {
            x: T,
            y: T,

            /// Calculate the (change in the) polarisation angle.
            /// E.g. Connors et al., 1980, Equation (2).
            pub fn angle(self: XY) T {
                return A.mult(.promote(0.5), A.atan2(self.y, self.x));
            }

            /// Calculate the polarisation degree.
            /// E.g. Connors et al., 1980, Equation (2).
            pub fn degree(self: XY) T {
                return A.sqrt(A.add(A.powi(self.x, 2), A.powi(self.y, 2)));
            }
        };

        const Self = @This();
        real: T,
        imag: T,

        /// Adapt the dual number type.
        pub fn adapt(
            self: *const Self,
            comptime NewT: type,
        ) PolarisationConstant(NewT) {
            return .{
                .real = .adaptFrom(self.real),
                .imag = .adaptFrom(self.imag),
            };
        }

        /// Calculate the X and Y components of the polarisation vector from
        /// the polarisation constant for an observer at infinity.
        ///
        /// E.g. Connors et al., 1977, Equation (7), specialised for an observer
        /// at infinity.
        pub fn polarisationXYAtInfinity(
            self: Self,
            metric: KerrMetric(T),
            theta_obs: T,
            eta: T,
            lambda: T,
            sign_theta: T.T,
        ) XY {
            const sin_theta = A.sin(theta_obs);
            const cos_theta = A.cos(theta_obs);

            const alpha = A.div(lambda, sin_theta);
            const beta_squared = A.sub(
                eta,
                A.mult(
                    A.powi(cos_theta, 2),
                    A.sub(A.powi(alpha, 2), A.powi(metric.a, 2)),
                ),
            );
            const beta = A.mult(
                .promote(-sign_theta),
                A.sqrt(beta_squared),
            );

            return self.polarisationXYAtInfinityAlt(metric, theta_obs, alpha, beta);
        }

        /// Alternative signature for `polarisationXYAtInfinity`.
        pub fn polarisationXYAtInfinityAlt(
            self: Self,
            metric: KerrMetric(T),
            theta_obs: T,
            alpha: T,
            beta: T,
        ) XY {
            const sin_theta = A.sin(theta_obs);

            // Note that in Connors et al. 1977 Equation (6) defines the
            // constant as
            //
            //    k_pw = k_2 - i k_1
            //
            // and therefore everywhere in the equation below where there is
            // `imag`, it must instead be `-imag`:
            const k1 = self.imag.neg();
            const k2 = self.real;

            const _X = A.add(
                A.mult(A.sub(alpha, A.mult(metric.a, sin_theta)), k1),
                A.mult(beta, k2),
            );
            const _Y = A.sub(
                A.mult(A.sub(alpha, A.mult(metric.a, sin_theta)), k2),
                A.mult(beta, k1),
            );

            // There are many equivalent ways to calculate the normalisation.
            // In Connors, it is written in terms of `S` and `T`, which are
            // proxies for the impact parameters. The below is equivalent but
            // less fiddly:
            const norm = A.sqrt(A.add(A.powi(_X, 2), A.powi(_Y, 2)));

            const x = A.div(_X, norm);
            const y = A.div(_Y, norm);

            // Sanity checks:
            std.debug.assert(@abs(x.x) <= 1);
            std.debug.assert(@abs(y.x) <= 1);

            return .{ .x = x, .y = y };
        }
    };
}

/// Compute the polarisation invariant constant, also known as the
/// Walker-Penrose constant. Here, `k` is the (photon) wave-vector, that is,
/// the four momentum at `x`, and `f` is the polarisation vector, as defined in
/// e.g. Connors et al., 1980 Equation (6).
pub fn polarisationConstant(
    comptime T: type,
    metric: KerrMetric(T),
    r: T,
    theta: T,
    k: FourVector(T),
    f: FourVector(T),
) PolarisationConstant(T) {
    const A = T.Algebra;
    const sin_theta = A.sin(theta);
    const cos_theta = A.cos(theta);

    const term_1_1 = commutator(T, k.t, f.r, k.r, f.t);
    const term_1_2 = A.mult(
        A.mult(metric.a, A.powi(sin_theta, 2)),
        commutator(T, k.r, f.ph, k.ph, f.r),
    );

    const term_1 = A.add(term_1_1, term_1_2);

    const term_2_1 = A.mult(
        A.add(A.powi(r, 2), A.powi(metric.a, 2)),
        commutator(T, k.ph, f.th, k.th, f.ph),
    );
    const term_2_2 = A.mult(
        metric.a,
        commutator(T, k.t, f.th, k.th, f.t),
    );

    const term_2 = A.mult(sin_theta, A.sub(term_2_1, term_2_2));

    //  term_1 * r - a term_2 cos(th)
    const real_component = A.sub(
        A.mult(term_1, r),
        A.mult(metric.a, A.mult(term_2, cos_theta)),
    );
    // -term_1 a cos(th) - term_2 r
    const imag_component = A.add(
        A.mult(term_1, A.mult(metric.a, cos_theta)),
        A.mult(term_2, r),
    ).neg();

    return .{ .real = real_component, .imag = imag_component };
}

test "polarisation constant" {
    const Dual = ad.DualNumber(f64, 0);
    const FV = FourVector(Dual);
    const metric: KerrMetric(Dual) = .init(
        .promote(1.0),
        .promote(0.998),
    );

    // these components calculated to be null and orthogonal to f
    const k: FV = .{
        .t = .promote(0.3038941082242305),
        .r = .promote(0.03),
        .th = .promote(-0.052578658425038796),
        .ph = .promote(0.02),
    };

    const f: FV = .{
        .t = .promote(0.0),
        .r = .promote(0.05),
        .th = .promote(0.03),
        .ph = .promote(0.47864754555328204),
    };

    const constant = polarisationConstant(
        Dual,
        metric,
        .promote(4.0),
        .promote(std.math.degreesToRadians(30)),
        k,
        f,
    );

    try std.testing.expectApproxEqAbs(
        -0.11120654616428541,
        constant.real.x,
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        -0.87367637659737,
        constant.imag.x,
        TEST_TOLERANCE,
    );
}

/// Convert the polarisaton constant back to a polarisation vector at a
/// particular point in the spacetime.
///
/// Here, `k` is the null four-velocity of the photon.
pub fn vectorFromPolarisationConstant(
    comptime T: type,
    metric: KerrMetric(T),
    r: T,
    theta: T,
    k: FourVector(T),
    pc: PolarisationConstant(T),
) FourVector(T) {
    const A = T.Algebra;
    const ts = metric.tangentSpaceAlt(r, theta);
    const k_down = ts.lowerIndices(k);

    const sin_theta = A.sin(theta);
    const cos_theta = A.cos(theta);
    const sin_theta_squared = A.powi(sin_theta, 2);
    const a_sin_theta_squared = A.mult(metric.a, sin_theta_squared);
    const r2_a2 = A.add(A.powi(r, 2), A.powi(metric.a, 2));

    const _C1 = A.sub(k.t, A.mult(a_sin_theta_squared, k.ph));
    const _C2 = A.mult(a_sin_theta_squared, k.r);
    const _C3 = A.mult(A.sub(A.mult(r2_a2, k.ph), A.mult(metric.a, k.t)), sin_theta);
    const _C4 = A.mult(A.mult(r2_a2, k.th), sin_theta).neg();

    const a_cos_theta = A.mult(metric.a, cos_theta);

    const _R1 = A.mult(r, _C1);
    const _R2 = A.mult(a_cos_theta, _C3).neg();
    const _R3 = A.sub(A.mult(r, _C2), A.mult(a_cos_theta, _C4));

    const _J1 = A.mult(a_cos_theta, _C1).neg();
    const _J2 = A.mult(r, _C3).neg();
    const _J3 = A.add(A.mult(a_cos_theta, _C2), A.mult(r, _C4)).neg();

    const _kthph = A.div(k_down.th, k_down.ph);
    const _krph = A.div(k_down.r, k_down.ph);

    const term_R2_R3 = A.sub(_R2, A.mult(_kthph, _R3));
    const term_J2_J3 = A.sub(_J2, A.mult(_kthph, _J3));

    const ratio = A.div(term_J2_J3, term_R2_R3);

    const term_J1_J3 = A.sub(_J1, A.mult(_krph, _J3));
    const term_R1_R3 = A.sub(_R1, A.mult(_krph, _R3));

    const f_r_bracket = A.sub(term_J1_J3, A.mult(term_R1_R3, ratio));
    const f_r_constant = A.sub(pc.imag, A.mult(pc.real, ratio));

    const f_r = A.div(f_r_constant, f_r_bracket);

    const f_theta = A.div(A.sub(pc.real, A.mult(f_r, term_R1_R3)), term_R2_R3);
    const f_phi = A.div(
        A.add(A.mult(k_down.r, f_r), A.mult(k_down.th, f_theta)),
        k_down.ph,
    ).neg();

    return .{
        .t = .zero,
        .r = f_r,
        .th = f_theta,
        .ph = f_phi,
    };
}

test "polarisation constant to vector" {
    const Dual = ad.DualNumber(f64, 0);
    const FV = FourVector(Dual);
    const metric: KerrMetric(Dual) = .init(
        .promote(1.0),
        .promote(0.998),
    );

    // these components calculated to be null and orthogonal to f
    const k: FV = .{
        .t = .promote(0.3038941082242305),
        .r = .promote(0.03),
        .th = .promote(-0.052578658425038796),
        .ph = .promote(0.02),
    };
    const f: FV = .{
        .t = .promote(0.0),
        .r = .promote(0.05),
        .th = .promote(0.03),
        .ph = .promote(0.47864754555328204),
    };

    const pc: PolarisationConstant(Dual) = .{
        .real = .promote(-0.11120654616428541),
        .imag = .promote(-0.87367637659737),
    };

    const vec = vectorFromPolarisationConstant(
        Dual,
        metric,
        .promote(4.0),
        .promote(std.math.degreesToRadians(30.0)),
        k,
        pc,
    );

    try std.testing.expectApproxEqAbs(
        f.t.x,
        vec.t.x,
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        f.r.x,
        vec.r.x,
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        f.th.x,
        vec.th.x,
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        f.ph.x,
        vec.ph.x,
        TEST_TOLERANCE,
    );
}

/// Calculate an initial polarisation vector f^mu in the global coordinates.
/// The polarisation vector is calculated as the vector perpendicular to the
/// photon four-velocity k^mu and the disc normal, i.e. e_theta unit vector.
/// The polarisation vector is further normalised to have unity magnitude.
pub fn initialPolarisationVector(
    comptime T: type,
    ts: KerrMetric(T).TangentSpace,
    disc_normal: FourVector(T),
    disc_velocity: FourVector(T),
    photon_k: FourVector(T),
) FourVector(T) {
    const A = T.Algebra;
    const u_dot_k = disc_velocity.dot(ts, photon_k);
    const n_dot_k = disc_normal.dot(ts, photon_k);

    // E.g. Johannsen 2015, Equation (A2), except my convention is to not have
    // the minus sign.
    const cos_theta = A.div(n_dot_k, u_dot_k);

    const k = photon_k.toArray();
    const n = disc_normal.toArray();
    const u = disc_velocity.toArray();

    // Similar to Dovciak, Muleri, et al. 2008 Equation (14), except I use a
    // different sign convention. For my derivation, see TODO: link to blog.
    const denominator = A.sqrt(A.sub(.one, A.powi(cos_theta, 2)));
    var out: [4]T = FourVector(T).zeros.toArray();
    for (0..4) |i| {
        const term = A.mult(cos_theta, A.add(A.div(k[i], u_dot_k), u[i]));
        out[i] = A.div(A.sub(n[i], term), denominator);
    }

    const f = FourVector(T).fromArray(out);

    if (@import("builtin").mode == .Debug) {
        // Sanity check the calculation:
        const f_magnitude = f.dot(ts, f);
        std.debug.assert(
            std.math.approxEqAbs(T.T, 1.0, f_magnitude.x, TEST_TOLERANCE),
        );
        const f_perp = f.dot(ts, photon_k);
        std.debug.assert(
            std.math.approxEqAbs(T.T, 0.0, f_perp.x, TEST_TOLERANCE),
        );
    }

    return f;
}

test "initial polarisation vector" {
    const Dual = ad.DualNumber(f64, 0);
    const FV = FourVector(Dual);
    const metric: KerrMetric(Dual) = .init(.promote(1.0), .promote(0.998));
    const x: FV = .{
        .t = .zero,
        .r = .promote(3.0),
        .th = .promote(std.math.pi / 2.0),
        .ph = .zero,
    };
    const disc_u: FV = .{
        .t = .promote(1.923360437276991),
        .r = .zero,
        .th = .zero,
        .ph = .promote(0.3105122873996937),
    };
    const ts = metric.tangentSpace(x);
    const disc_normal: FV = .{
        .t = .zero,
        .r = .zero,
        .th = .promote(1.0 / @sqrt(ts.metric_components.thth.x)),
        .ph = .zero,
    };

    // This is a sanity check
    const n_dot_n = disc_normal.dot(ts, disc_normal);
    try std.testing.expectApproxEqAbs(1.0, n_dot_n.x, TEST_TOLERANCE);
    const u_dot_u = disc_u.dot(ts, disc_u);
    try std.testing.expectApproxEqAbs(-1.0, u_dot_u.x, TEST_TOLERANCE);
    const u_dot_n = disc_u.dot(ts, disc_normal);
    try std.testing.expectApproxEqAbs(0.0, u_dot_n.x, TEST_TOLERANCE);

    const k: FV = .{
        .t = .promote(1.4461370272644967),
        .r = .promote(-0.1139497538585655),
        .th = .promote(0.28867513459481287),
        .ph = .promote(0.11126501553002391),
    };

    const f = initialPolarisationVector(Dual, ts, disc_normal, disc_u, k);
    // check the properties
    const f_magnitude = f.dot(ts, f);
    try std.testing.expectApproxEqAbs(1.0, f_magnitude.x, TEST_TOLERANCE);
    const f_perp = f.dot(ts, k);
    try std.testing.expectApproxEqAbs(0.0, f_perp.x, TEST_TOLERANCE);
}

test "polarisation transport" {
    const Dual = ad.DualNumber(f64, 0);
    const FV = FourVector(Dual);
    const metric: KerrMetric(Dual) = .init(.promote(1.0), .promote(0.998));
    const x: FV = .{
        .t = .zero,
        .r = .promote(3.0),
        .th = .promote(std.math.pi / 2.0),
        .ph = .zero,
    };

    const f: FV = .{
        .t = .promote(0.8265751927029137),
        .r = .promote(0.1973667631930029),
        .th = .promote(0.16666666666666677),
        .ph = .promote(0.34510639814775684),
    };

    const ts = metric.tangentSpace(x);
    const k = geodesic.skyAnglesToVelocity(
        Dual,
        ts,
        orbits.stationary(Dual, ts),
        .promote(0.2),
        .promote(0.4),
    );

    const pc = polarisationConstant(Dual, metric, x.r, x.th, k, f);

    // now trace the geodesic some distance
    const geod = geodesic.NullGeodesic(Dual).fromNullVelocityTangentSpace(
        metric,
        ts,
        k,
    );
    const result = geod.traceBuilder(metric, .{}).atMinoTime(.promote(0.05));

    const k_new = result.velocity(metric, geod);

    // Calculate the new polarisation vector at this location
    const f_new = vectorFromPolarisationConstant(
        Dual,
        metric,
        result.r,
        result.theta,
        k_new,
        pc,
    );

    // then using this to calculate a new polarisation constant should give the
    // same constant
    const pc_new = polarisationConstant(
        Dual,
        metric,
        result.r,
        result.theta,
        k_new,
        f_new,
    );

    try std.testing.expectApproxEqAbs(pc.imag.x, pc_new.imag.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(pc.real.x, pc_new.real.x, TEST_TOLERANCE);
}

test "polarisation projections at infinity" {
    const Dual = ad.DualNumber(f64, 0);
    const FV = FourVector(Dual);
    const metric: KerrMetric(Dual) = .init(.promote(1.0), .promote(0.998));

    const x: FV = .{
        .t = .zero,
        .r = .promote(1e7),
        .th = .promote(std.math.degreesToRadians(10)),
        .ph = .zero,
    };

    const geod = geodesic.NullGeodesic(Dual).fromImpactParameters(
        metric,
        x,
        .promote(2.0),
        .promote(-3.0),
    );
    const result = geod.traceDisc(metric, .equatorial_plane, .{});

    const v_medium = orbits.keplerianPlungingAlt(Dual, metric, result.r);

    const pc = result.polarisationConstant(metric, geod, v_medium);

    const k_obs = geod.initialVelocity(metric);

    // Calculate the new polarisation vector at the origin, since we know the
    // constant at the disc.
    const f_new = vectorFromPolarisationConstant(
        Dual,
        metric,
        geod.x_init.r,
        geod.x_init.th,
        k_obs,
        pc,
    );

    const xy = pc.polarisationXYAtInfinity(
        metric,
        x.th,
        geod.eta,
        geod.lambda,
        @floatFromInt(geod.theta_sign),
    );

    // The polarisation component in the `r` direction should be zero at the
    // observer.
    try std.testing.expectApproxEqAbs(f_new.r.x, 0.0, TEST_TOLERANCE);

    // TODO: I have not checked these values, so they are currently simply
    // regression tests:
    try std.testing.expectApproxEqAbs(
        -0.16128329132425687,
        xy.x.x,
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        0.0000005683377625008126,
        f_new.ph.x,
        TEST_TOLERANCE,
    );

    const polarisation_angle = result.polarisationAngle(metric, geod, v_medium);
    try std.testing.expectApproxEqAbs(
        -1.4533308628783543,
        polarisation_angle.x,
        TEST_TOLERANCE,
    );
}
