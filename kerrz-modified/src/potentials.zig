/// References:
/// B72: Bardeen et al. 1972
const std = @import("std");
const ad = @import("zad");
const tracy = @import("tracy.zig");

const TEST_TOLERANCE = @import("options").test_numerical_tolerance;

const root = @import("root.zig");
const ComplexNumber = @import("complex.zig").ComplexNumber;
const KerrMetric = root.KerrMetric;
const FourVector = root.FourVector;

/// Compute the angular potential.
/// This is the V_theta of B72, Equation (2.10), for the case of mu = 0.
pub fn angular(comptime T: type, geometry: KerrMetric(T), eta: T, lambda: T, theta: T) T {
    const A = T.Algebra;
    const sin_theta_squared = A.powi(A.sin(theta), 2);
    const cos_theta_squared = A.sub(.one, sin_theta_squared);

    const bracket = A.sub(
        A.powi(geometry.a, 2),
        A.div(A.powi(lambda, 2), sin_theta_squared),
    );

    return A.add(eta, A.mult(cos_theta_squared, bracket));
}

test "angular potential" {
    const Dual = ad.DualNumber(f64, 0);
    const result = angular(
        Dual,
        .init(.one, .promote(0.998)),
        .promote(0.2),
        .promote(0.3),
        .promote(std.math.degreesToRadians(40.0)),
    );
    try std.testing.expectApproxEqAbs(0.6566542434829781, result.x, TEST_TOLERANCE);
}

/// Compute the radial potential.
/// This is the V_r of B72, Equation (2.10) with mu = 0 or G&L Equation (7).
pub fn radial(comptime T: type, geometry: KerrMetric(T), eta: T, lambda: T, r: T) T {
    const A = T.Algebra;

    const delta = geometry.delta(r);

    const term_1 = A.powi(A.sub(
        A.add(
            A.powi(r, 2),
            A.powi(geometry.a, 2),
        ),
        A.mult(geometry.a, lambda),
    ), 2);

    const term_2 = A.mult(
        delta,
        A.add(eta, A.powi(A.sub(lambda, geometry.a), 2)),
    );

    return A.sub(term_1, term_2);
}

/// Compute the radial potential using the roots of the potential, e.g. G&L
/// Equation (B8). The roots must all be real numbers. The roots must all be
/// real numbers (this is not checked by this function!).
///
/// Use `radialFromRootsAlt` for an alternative interface.
pub fn radialFromRoots(comptime T: type, roots: RadialRoots(T), r: T) T {
    return radialFromRootsAlt(
        T,
        roots.r1.real(),
        roots.r2.real(),
        roots.r3.real(),
        roots.r4.real(),
        r,
    );
}

/// Alternative interface for `radialFromRoots`.
pub fn radialFromRootsAlt(comptime T: type, r1: T, r2: T, r3: T, r4: T, r: T) T {
    const A = T.Algebra;
    return A.mult(
        A.mult(
            A.sub(r, r1),
            A.sub(r, r2),
        ),
        A.mult(
            A.sub(r, r3),
            A.sub(r, r4),
        ),
    );
}

test "radial potential" {
    const Dual = ad.DualNumber(f64, 1);
    const result = radial(
        Dual,
        .init(.one, .promote(0.998)),
        .promote(0.5),
        .promote(0.3),
        .derivative(10.0),
    );
    try std.testing.expectApproxEqAbs(10059.846478, result.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(4010.0944879999993, result.dx[0], TEST_TOLERANCE);
}

pub const AngularCase = enum(u8) {
    normal,
    vortical,

    /// Determine the angular case from the sign of the eta reduced constant of
    /// motion related to the Carter constant.
    pub fn fromSign(eta_sign: i32) AngularCase {
        if (eta_sign > 0) {
            return .normal;
        } else {
            return .vortical;
        }
    }
};

pub const RadialCase = enum(u8) {
    case_I,
    case_II,
    case_III,
    case_IV,
};

pub fn RadialRoots(comptime T: type) type {
    const Complex = ComplexNumber(T);
    return struct {
        const Self = @This();
        r1: Complex,
        r2: Complex,
        r3: Complex,
        r4: Complex,

        /// Adapt the dual number type.
        pub fn adapt(self: Self, comptime NewT: type) RadialRoots(NewT) {
            return .{
                .r1 = self.r1.adapt(NewT),
                .r2 = self.r2.adapt(NewT),
                .r3 = self.r3.adapt(NewT),
                .r4 = self.r4.adapt(NewT),
            };
        }

        fn construct(xi0: T, _A: T, _B: T) Self {
            const A = T.Algebra;
            const z = A.sqrt(A.mult(xi0, .promote(0.5)));

            const term_1 = A.sub(A.mult(.promote(-0.5), _A), A.powi(z, 2));
            const term_2 = A.div(_B, A.mult(.promote(4), z));

            const c1 = Complex.fromReIm(
                A.add(term_1, term_2),
                .zero,
            ).sqrt();
            const c2 = Complex.fromReIm(
                A.sub(term_1, term_2),
                .zero,
            ).sqrt();

            const r1 = c1.neg().addRe(z.neg());
            const r2 = c1.addRe(z.neg());
            const r3 = c2.neg().addRe(z);
            const r4 = c2.addRe(z);

            return .{ .r1 = r1, .r2 = r2, .r3 = r3, .r4 = r4 };
        }

        /// Determine which of the radial motion cases these roots correspond
        /// to.
        pub fn determineCase(self: *const Self, horizon: T, r: T) RadialCase {
            const error_tolerance = @sqrt(std.math.floatEps(T.T));
            _ = horizon;
            if (!self.r2.isReal()) {
                std.debug.assert(!self.r1.isReal());
                std.debug.assert(!self.r3.isReal());
                std.debug.assert(!self.r4.isReal());
                return .case_IV;
            }

            if (!self.r4.isReal()) {
                std.debug.assert(self.r1.isReal());
                std.debug.assert(self.r2.isReal());
                std.debug.assert(!self.r3.isReal());
                return .case_III;
            }

            std.debug.assert(self.r1.isReal());
            std.debug.assert(self.r2.isReal());
            std.debug.assert(self.r3.isReal());
            const r4 = self.r4.real().x;

            if (r4 <= r.x + error_tolerance) {
                return .case_II;
            }

            return .case_I;
        }
    };
}

/// Helper function that computes the cube root if the argument is real, else
/// the largest real part of the complex solutions.
fn principleCuberoot(comptime T: type, arg: T) T {
    const A = T.Algebra;

    // TODO: fix the arg == 0 case, to avoid the derivative blowing up.
    if (arg.x > 0) {
        return A.cuberoot(arg);
    }
    if (arg.x == 0) {
        return .zero;
    }

    const complex = ComplexNumber(T).fromReIm(arg, .zero);

    const r1 = complex.cuberoot();

    // TODO: there's a pretty good chance r1 is always maximal real, but it
    // would be good to check before optimizing this
    const r2 = r1.rotateArg(.promote(2 * std.math.pi / 3.0));
    const r3 = r1.rotateArg(.promote(-2 * std.math.pi / 3.0));

    const rr1 = r1.real();
    const rr2 = r2.real();
    const rr3 = r3.real();
    if (rr1.x >= rr2.x and rr1.x >= rr3.x) {
        return rr1;
    }
    if (rr2.x >= rr1.x and rr2.x >= rr3.x) {
        return rr1;
    }
    if (rr3.x >= rr1.x and rr3.x >= rr2.x) {
        return rr3;
    }
    unreachable;
}

/// Find the four roots of the radial potential.
pub fn rootsOfRadialPotential(
    comptime T: type,
    m: KerrMetric(T),
    eta: T,
    lambda: T,
) RadialRoots(T) {
    var ctx = tracy.trace(@src());
    defer ctx.end();
    const A = T.Algebra;

    // The caligraphic expressions
    // G&L: Equation (79)
    const _A = A.sub(A.sub(A.powi(m.a, 2), eta), A.powi(lambda, 2));
    // G&L: Equation (80)
    const _B = A.mult(
        A.mult(.promote(2), m.M),
        A.add(
            eta,
            A.powi(A.sub(lambda, m.a), 2),
        ),
    );
    std.debug.assert(_B.x > 0);
    // G&L: Equation (81)
    const _C = A.mult(A.powi(m.a, 2).neg(), eta);

    // G&L: Equation (85)
    const _P = A.sub(A.mult(.promote(-1.0 / 12.0), A.powi(_A, 2)), _C);
    // G&L: Equation (86)
    const _Q = A.sub(
        A.mult(
            A.mult(_A, .promote(-1.0 / 3.0)),
            A.sub(
                A.powi(A.mult(_A, .promote(1.0 / 6.0)), 2),
                _C,
            ),
        ),
        A.mult(A.powi(_B, 2), .promote(1.0 / 8.0)),
    );

    // G&L: Equation (92)
    const discriminant = A.add(
        A.mult(.promote(-4), A.powi(_P, 3)),
        A.mult(.promote(-27), A.powi(_Q, 2)),
    );

    const neg_Q_half = A.mult(_Q, .promote(-1.0 / 2.0));

    if (discriminant.x <= 0) {
        // The roots are all real, no special considerations needed
        const shared = A.sqrt(A.mult(discriminant, .promote(-1.0 / 108.0)));
        // G&L: Equation (91)
        const w_plus = principleCuberoot(T, A.add(neg_Q_half, shared));
        const w = if (w_plus.x == 0)
            principleCuberoot(T, A.sub(neg_Q_half, shared))
        else
            A.div(_P.neg(), A.mult(.promote(3), w_plus));

        const xi0 = A.add(
            A.add(w_plus, w),
            A.mult(_A, .promote(-1.0 / 3.0)),
        );

        // G&L: Required by Equation (93)
        std.debug.assert(xi0.x > 0);

        return .construct(xi0, _A, _B);
    } else {
        // Roots are complex. Calculate in in A e^(i B) format, but
        // since since the result is has manifestly no real part, B =
        // pi at this stage. All other operations are just rotating the
        // argument by thirds, so we don't need to be too clever, and
        // can exploit the fact that w+ and w- are complex conjugate,
        // such that
        //
        //    (w+) + (w-) = 2 * Re[w+] = 2 * Re[w-]
        //
        // is purely real.
        const imag = A.sqrt(A.mult(discriminant.neg(), .promote(-1.0 / 108.0)));

        // perform the cube root
        const w = ComplexNumber(T).fromReIm(
            neg_Q_half,
            imag,
        ).cuberoot().real();

        const xi0 = A.add(
            A.mult(w, .promote(2.0)),
            A.mult(_A, .promote(-1.0 / 3.0)),
        );

        // G&L: Required by Equation (93)
        std.debug.assert(xi0.x > 0);

        return .construct(xi0, _A, _B);
    }
}

test "radial roots" {
    const Dual = ad.DualNumber(f64, 1);
    const geom: KerrMetric(Dual) = .init(.one, .promote(0.998));
    const x: FourVector(Dual) = .{
        .t = .zero,
        .r = .promote(1e4),
        .th = .promote(std.math.degreesToRadians(60.0)),
        .ph = .zero,
    };

    {
        // this hits the imaginary roots branch
        const geod = root.NullGeodesic(Dual).fromConstantsOfMotion(
            x,
            .promote(-0.9998999950002677),
            .promote(-9.999995007953801),
            .promote(99.77344893700447),
        );
        const roots = rootsOfRadialPotential(Dual, geom, geod.eta, geod.lambda);
        try std.testing.expectApproxEqAbs(-14.948516116049726, roots.r1.real().x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.33730782701960926, roots.r2.real().x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(1.5039233086872947, roots.r3.real().x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(13.107284980342822, roots.r4.real().x, TEST_TOLERANCE);

        try std.testing.expectApproxEqAbs(0.0, roots.r1.imag().x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.0, roots.r2.imag().x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.0, roots.r3.imag().x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.0, roots.r4.imag().x, TEST_TOLERANCE);
    }
    {
        // now with auto-diff
        var Q = Dual.promote(99.77344893700447);
        Q.dx[0] = 1.0;
        const geod = root.NullGeodesic(Dual).fromConstantsOfMotion(
            x,
            .promote(-0.9998999950002677),
            .promote(-9.999995007953801),
            Q,
        );
        const roots = rootsOfRadialPotential(Dual, geom, geod.eta, geod.lambda);
        try std.testing.expectApproxEqAbs(-14.948516116049726, roots.r1.real().x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.33730782701960926, roots.r2.real().x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(1.5039233086872947, roots.r3.real().x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(13.107284980342822, roots.r4.real().x, TEST_TOLERANCE);

        try std.testing.expectApproxEqAbs(-0.036048891322107715, roots.r1.real().dx[0], 1e-3);
        try std.testing.expectApproxEqAbs(0.0019109405025805604, roots.r2.real().dx[0], 1e-3);
        try std.testing.expectApproxEqAbs(-0.0011222732040317672, roots.r3.real().dx[0], 1e-3);
        try std.testing.expectApproxEqAbs(0.035260224023558925, roots.r4.real().dx[0], 1e-3);
    }

    {
        const new_x: FourVector(Dual) = .{
            .t = .zero,
            .r = .promote(3),
            .th = .promote(std.math.degreesToRadians(80.0)),
            .ph = .zero,
        };
        const v: FourVector(Dual) = .{
            .t = .zero,
            .r = .promote(-0.2682964555194517),
            .th = .promote(-0.13923855392539625),
            .ph = .promote(0.23945893786208244),
        };

        // this hits the imaginary roots branch
        const geod = root.NullGeodesic(Dual).fromVelocity(
            geom,
            new_x,
            v,
        );

        try std.testing.expectApproxEqAbs(0.5792673090728395, geod.E.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(1.6554408085454932, geod.L.x, TEST_TOLERANCE);

        const roots = rootsOfRadialPotential(Dual, geom, geod.eta, geod.lambda);
        try std.testing.expectApproxEqAbs(-4.0660210976434925, roots.r1.real().x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.4155268953045834, roots.r2.real().x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(1.1754925500891664, roots.r3.real().x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(2.4750016522497424, roots.r4.real().x, TEST_TOLERANCE);

        try std.testing.expectApproxEqAbs(0.0, roots.r1.imag().x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.0, roots.r2.imag().x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.0, roots.r3.imag().x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.0, roots.r4.imag().x, TEST_TOLERANCE);
    }
}

pub fn AngularRoots(comptime T: type) type {
    return struct {
        const Self = @This();
        u_plus: T,
        u_minus: T,

        /// Adapt the dual number type.
        pub fn adapt(self: Self, comptime NewT: type) AngularRoots(NewT) {
            return .{
                .u_plus = .adaptFrom(self.u_plus),
                .u_minus = .adaptFrom(self.u_minus),
            };
        }
    };
}

/// Compute the two roots of the angular potential.
/// G&L: Equations (18) and (19).
pub fn rootsOfAngularPotential(
    comptime T: type,
    m: KerrMetric(T),
    eta: T,
    lambda: T,
) AngularRoots(T) {
    const error_tolerance = 100 * std.math.floatEps(T.T);
    var ctx = tracy.trace(@src());
    defer ctx.end();
    const A = T.Algebra;

    if (m.a.x == 0) {
        // Taking (18) and solving by hand
        const u = A.div(
            eta.neg(),
            A.sub(eta, A.powi(lambda, 2)),
        );
        return .{
            .u_plus = u,
            .u_minus = u,
        };
    }

    const discriminant = A.mult(
        .promote(1.0 / 2.0),
        A.sub(
            .one,
            A.div(
                A.add(eta, A.powi(lambda, 2)),
                A.powi(m.a, 2),
            ),
        ),
    );

    var common_squared = A.add(
        A.powi(discriminant, 2),
        A.div(eta, A.powi(m.a, 2)),
    );

    // Avoid rounding errors:
    if (common_squared.x < 0 and @abs(common_squared.x) < error_tolerance) {
        common_squared.x = 0;
    }

    const common = A.sqrt(common_squared);

    return .{
        .u_plus = A.add(discriminant, common),
        .u_minus = A.sub(discriminant, common),
    };
}

test "roots of angular potential" {
    const Dual = ad.DualNumber(f64, 0);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    const roots = rootsOfAngularPotential(Dual, metric, .promote(3.0042993693628364), .zero);
    try std.testing.expectApproxEqAbs(0.9999999999999998, roots.u_plus.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-3.016352714811222, roots.u_minus.x, TEST_TOLERANCE);
}
