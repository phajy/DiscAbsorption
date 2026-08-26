const std = @import("std");
const ad = @import("zad");
const options = @import("options");
const zfits = @import("zfits");

const SMatrix = @import("matrix.zig").SMatrix;
const orbits = @import("orbits.zig");

const TEST_TOLERANCE = options.test_numerical_tolerance;

pub fn AnglePair(comptime T: type) type {
    return struct {
        theta: T,
        phi: T,
    };
}

/// A four-vector of the Boyer-Lindquist (or other spherical) coordinates.
pub fn FourVector(comptime T: type) type {
    const _zero: T = if (@typeInfo(T) == .@"struct") T.zero else 0;
    return struct {
        const Self = @This();
        t: T,
        r: T,
        th: T,
        ph: T,

        /// Adapt the number type.
        pub fn adapt(self: Self, comptime NewT: type) FourVector(NewT) {
            return .{
                .t = .adaptFrom(self.t),
                .r = .adaptFrom(self.r),
                .th = .adaptFrom(self.th),
                .ph = .adaptFrom(self.ph),
            };
        }

        /// Convert to an array.
        pub fn toArray(self: Self) [4]T {
            return [4]T{ self.t, self.r, self.th, self.ph };
        }

        /// Get only the spatial components.
        pub fn toThreeVector(self: Self) ThreeVector(T) {
            return .{
                .v = .{ self.r, self.th, self.ph },
            };
        }

        /// Initialise from an array.
        pub fn fromArray(self: [4]T) Self {
            return .{
                .t = self[0],
                .r = self[1],
                .th = self[2],
                .ph = self[3],
            };
        }

        /// Use a numerical index to access a particular component.
        /// 0 is t, 1 is r, and so on.
        pub fn getIndex(self: Self, i: usize) T {
            const s: *const [4]T = @ptrCast(&self);
            return s.*[i];
        }

        // Used only in unit testing.
        fn testEqual(self: Self, other: Self) !void {
            try std.testing.expectApproxEqAbs(
                other.t.x,
                self.t.x,
                TEST_TOLERANCE,
            );
            try std.testing.expectApproxEqAbs(
                other.r.x,
                self.r.x,
                TEST_TOLERANCE,
            );
            try std.testing.expectApproxEqAbs(
                other.th.x,
                self.th.x,
                TEST_TOLERANCE,
            );
            try std.testing.expectApproxEqAbs(
                other.ph.x,
                self.ph.x,
                TEST_TOLERANCE,
            );
        }

        pub const zeros: Self = .{ .t = _zero, .ph = _zero, .r = _zero, .th = _zero };

        /// Sum all the components together.
        pub fn sum(self: Self) T {
            const A = T.Algebra;
            return A.add(self.t, A.add(self.r, A.add(self.th, self.ph)));
        }

        /// Elementwise compare to a scalar value.
        pub fn scalarEq(self: Self, value: T.T) [4]bool {
            return [4]bool{
                self.t.x == value,
                self.r.x == value,
                self.th.x == value,
                self.ph.x == value,
            };
        }

        /// Elementwise compare to a scalar value.
        pub fn scalarNotEq(self: Self, value: T.T) [4]bool {
            return [4]bool{
                self.t.x != value,
                self.r.x != value,
                self.th.x != value,
                self.ph.x != value,
            };
        }

        /// Add vectors together component wise.
        pub fn add(self: Self, other: Self) Self {
            const A = T.Algebra;
            return .{
                .t = A.add(self.t, other.t),
                .r = A.add(self.r, other.r),
                .th = A.add(self.th, other.th),
                .ph = A.add(self.ph, other.ph),
            };
        }

        /// Subtract vectors component wise.
        pub fn sub(self: Self, other: Self) Self {
            const A = T.Algebra;
            return .{
                .t = A.sub(self.t, other.t),
                .r = A.sub(self.r, other.r),
                .th = A.sub(self.th, other.th),
                .ph = A.sub(self.ph, other.ph),
            };
        }

        /// Multiply vectors together component wise.
        pub fn mult(self: Self, other: Self) Self {
            const A = T.Algebra;
            return .{
                .t = A.mult(self.t, other.t),
                .r = A.mult(self.r, other.r),
                .th = A.mult(self.th, other.th),
                .ph = A.mult(self.ph, other.ph),
            };
        }

        /// Multiply the vector by a scalar.
        pub fn scalarMult(self: Self, a: T) Self {
            const A = T.Algebra;
            return .{
                .t = A.mult(self.t, a),
                .r = A.mult(self.r, a),
                .th = A.mult(self.th, a),
                .ph = A.mult(self.ph, a),
            };
        }

        /// Compute the dot product of two vectors in the Kerr geometry.
        pub fn dot(self: Self, ts: KerrMetric(T).TangentSpace, other: Self) T {
            const A = T.Algebra;
            const m = ts.metric_components;

            var total: T = _zero;
            total = A.add(
                total,
                A.mult(A.add(A.mult(m.tt, self.t), A.mult(m.tph, self.ph)), other.t),
            );
            total = A.add(
                total,
                A.mult(A.mult(m.rr, self.r), other.r),
            );
            total = A.add(
                total,
                A.mult(A.mult(m.thth, self.th), other.th),
            );
            total = A.add(
                total,
                A.mult(A.add(A.mult(m.phph, self.ph), A.mult(m.tph, self.t)), other.ph),
            );
            return total;
        }

        /// Compute the norm in the Kerr geometry.
        pub fn properNorm(self: Self, m: KerrMetric(T).TangentSpace) T {
            return self.dot(m, self);
        }
    };
}

/// A conventional three-vector implementation.
pub fn ThreeVector(comptime T: type) type {
    return struct {
        const Self = @This();
        /// The convetion is either (x, y, z), (r, theta, phi), or (z, rho,
        /// theta).
        v: [3]T,

        /// Initialise the three-vector from components.
        pub fn init(v1: T, v2: T, v3: T) Self {
            return .{ .v = [3]T{ v1, v2, v3 } };
        }

        /// Turn into a FourVector with the time component given by the
        /// argument.
        pub fn toFourVector(self: Self, time: T) FourVector(T) {
            return .{
                .t = time,
                .r = self.v[0],
                .th = self.v[1],
                .ph = self.v[2],
            };
        }

        /// Multiply vector by a scalar.
        pub fn scalarMult(self: Self, scalar: T) Self {
            const A = T.Algebra;
            return .init(
                A.mult(scalar, self.v[0]),
                A.mult(scalar, self.v[1]),
                A.mult(scalar, self.v[2]),
            );
        }

        /// Add two vectors together.
        pub fn add(self: Self, other: Self) Self {
            const A = T.Algebra;
            return .init(
                A.add(other.v[0], self.v[0]),
                A.add(other.v[1], self.v[1]),
                A.add(other.v[2], self.v[2]),
            );
        }

        /// Return a unit-normalised version of this vector.
        pub fn toUnitVector(self: Self) Self {
            const A = T.Algebra;
            const inv_norm = A.div(
                .one,
                A.sqrt(self.dot(self)),
            );
            return scalarMult(self, inv_norm);
        }

        /// Cross product.
        pub fn cross(self: Self, other: Self) Self {
            const A = T.Algebra;
            const r0 = A.sub(A.mult(self.v[1], other.v[2]), A.mult(self.v[2], other.v[1]));
            const r1 = A.sub(A.mult(self.v[2], other.v[0]), A.mult(self.v[0], other.v[2]));
            const r2 = A.sub(A.mult(self.v[0], other.v[1]), A.mult(self.v[1], other.v[0]));
            return .init(r0, r1, r2);
        }

        /// Dot product.
        pub fn dot(self: Self, other: Self) T {
            const A = T.Algebra;
            var out: T = .zero;
            inline for (0..3) |i| {
                out = A.add(out, A.mult(self.v[i], other.v[i]));
            }
            return out;
        }
    };
}

test "three vector" {
    const Dual = ad.DualNumber(f64, 0);

    const v1: ThreeVector(Dual) = .init(.promote(2.0), .promote(-3.0), .promote(0.5));
    const v2: ThreeVector(Dual) = .init(.promote(-9.3), .promote(1.0), .promote(2.5));

    try std.testing.expectApproxEqAbs(-20.35, v1.dot(v2).x, TEST_TOLERANCE);

    const cross = v1.cross(v2);

    try std.testing.expectApproxEqAbs(-8.0, cross.v[0].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-9.65, cross.v[1].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-25.9, cross.v[2].x, TEST_TOLERANCE);

    const add = v1.add(v2);

    try std.testing.expectApproxEqAbs(-7.3, add.v[0].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-2.0, add.v[1].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(3.0, add.v[2].x, TEST_TOLERANCE);

    const mult = v1.scalarMult(.promote(-2));

    try std.testing.expectApproxEqAbs(-4.0, mult.v[0].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(6.0, mult.v[1].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-1.0, mult.v[2].x, TEST_TOLERANCE);
}

/// For transforming a spherical (tangent) vector at a given point to a
/// cartesian (tangent) vector. See also `SphericalFromCartesian`.
pub fn CartesianFromSpherical(comptime T: type) type {
    return struct {
        const Self = @This();
        jacobian: SMatrix(T, 3, 3),

        fn cartesianFromSphericalJacobian(r: T, th: T, ph: T) SMatrix(T, 3, 3) {
            const A = T.Algebra;
            const sin_th = A.sin(th);
            const cos_th = A.cos(th);
            const sin_ph = A.sin(ph);
            const cos_ph = A.cos(ph);

            return .{
                .v = [9]T{
                    // Row 1
                    A.mult(r, A.mult(sin_th, cos_ph)),
                    A.mult(r, A.mult(sin_th, sin_ph)),
                    A.mult(r, cos_th),
                    // Row 2
                    A.mult(r, A.mult(cos_th, cos_ph)),
                    A.mult(r, A.mult(cos_th, sin_ph)),
                    A.mult(r, sin_th).neg(),
                    // Row 3
                    A.mult(r, sin_ph).neg(),
                    A.mult(r, cos_ph),
                    .zero,
                },
            };
        }

        fn init(r: T, th: T, ph: T) Self {
            return .{
                .jacobian = cartesianFromSphericalJacobian(r, th, ph),
            };
        }

        /// Transform a spherical vector to a cartesian vector.
        pub fn transform(_: Self, vec: ThreeVector(T)) ThreeVector(T) {
            return transformAlt(vec);
        }

        /// Transform a spherical vector to a cartesian vector.
        pub fn transformAlt(vec: ThreeVector(T)) ThreeVector(T) {
            const A = T.Algebra;
            const sin_th = A.sin(vec.v[1]);
            const cos_th = A.cos(vec.v[1]);
            const sin_ph = A.sin(vec.v[2]);
            const cos_ph = A.cos(vec.v[2]);
            return .{ .v = [3]T{
                A.mult(vec.v[0], A.mult(sin_th, cos_ph)),
                A.mult(vec.v[0], A.mult(sin_th, sin_ph)),
                A.mult(vec.v[0], cos_th),
            } };
        }

        /// Transform a spherical vector to a cartesian vector at the location
        /// specified. Mathematically, this is achieved by also multiplying by
        /// the Jacobian of the transformation, which maps how the unit vector
        /// change in the coordinate transformation.
        pub fn apply(self: Self, vec: ThreeVector(T)) ThreeVector(T) {
            const v_prime = self.transform(vec);
            return .{
                .v = self.jacobian.multR(v_prime.v),
            };
        }
    };
}

/// Get a transformation which maps spherical vectors at some point `(r, theta,
/// phi)` in cartesian polar coordinates to a cartesian vector.
pub fn cartesianFromSpherical(comptime T: type, th: T, ph: T) CartesianFromSpherical(T) {
    return .init(.one, th, ph);
}

test "cartesian from spherical" {
    const Dual = ad.DualNumber(f64, 0);

    const xfm = cartesianFromSpherical(
        Dual,
        .promote(std.math.degreesToRadians(30)),
        .promote(std.math.degreesToRadians(60)),
    );

    const p: ThreeVector(Dual) = .{
        .v = [3]Dual{ .one, .promote(0.2), .promote(-0.8) },
    };

    {
        const result = xfm.transform(p);
        try std.testing.expectApproxEqAbs(0.13841425570643057, result.v[0].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.1425166545207693, result.v[1].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.9800665778412416, result.v[2].x, TEST_TOLERANCE);
    }

    {
        const result = xfm.apply(p);
        try std.testing.expectApproxEqAbs(0.8216545960985236, result.v[0].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.5369856489673557, result.v[1].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.19112858894806875, result.v[2].x, TEST_TOLERANCE);
    }
}

/// For transforming a cartesian (tangent) vector at a given point to a
/// spherical (tangent) vector. See also `CartesianFromSpherical`.
pub fn SphericalFromCartesian(comptime T: type) type {
    return struct {
        const Self = @This();
        jacobian: SMatrix(T, 3, 3),

        fn sphericalFromCartesianJacobian(r: T, th: T, ph: T) SMatrix(T, 3, 3) {
            std.debug.assert(r.x == 1.0);

            const A = T.Algebra;
            const sin_th = A.sin(th);
            const cos_th = A.cos(th);
            const sin_ph = A.sin(ph);
            const cos_ph = A.cos(ph);

            return .{
                .v = [9]T{
                    // Row 1
                    A.mult(r, A.mult(sin_th, cos_ph)),
                    A.mult(r, A.mult(cos_th, cos_ph)),
                    A.mult(r, sin_ph).neg(),
                    // Row 2
                    A.mult(r, A.mult(sin_th, sin_ph)),
                    A.mult(r, A.mult(cos_th, sin_ph)),
                    A.mult(r, cos_ph),
                    // Row 3
                    A.mult(r, cos_th),
                    A.mult(r, sin_th).neg(),
                    .zero,
                },
            };
        }

        fn init(r: T, th: T, ph: T) Self {
            return .{
                .jacobian = sphericalFromCartesianJacobian(r, th, ph),
            };
        }

        /// Transform a cartesian vector to a spherical vector. The input
        /// vector should be `x, y, z` and the output vector is `r, theta,
        /// phi`.
        pub fn transform(_: Self, vec: ThreeVector(T)) ThreeVector(T) {
            return transformAlt(vec);
        }

        pub fn transformAlt(vec: ThreeVector(T)) ThreeVector(T) {
            const A = T.Algebra;

            const r = A.sqrt(A.add(
                A.add(
                    A.powi(vec.v[0], 2),
                    A.powi(vec.v[1], 2),
                ),
                A.powi(vec.v[2], 2),
            ));
            const theta = A.acos(A.div(vec.v[2], r));
            const phi = A.atan2(vec.v[1], vec.v[0]);

            return .{ .v = [3]T{ r, theta, phi } };
        }

        /// Transform a cartesian vector to a spherical vector at the location
        /// specified. Mathematically, this is achieved by also multiplying by
        /// the Jacobian of the transformation, which maps how the unit vector
        /// change in the coordinate transformation.
        pub fn apply(self: Self, vec: ThreeVector(T)) ThreeVector(T) {
            const v_prime = self.jacobian.multR(vec.v);
            return self.transform(.{ .v = v_prime });
        }
    };
}

test "spherical from cartesian" {
    const Dual = ad.DualNumber(f64, 0);

    const xfm = sphericalFromCartesian(
        Dual,
        .promote(std.math.degreesToRadians(30)),
        .promote(std.math.degreesToRadians(60)),
    );

    {
        const p: ThreeVector(Dual) = .{ .v = [3]Dual{
            .promote(3.0),
            .promote(-0.4),
            .promote(0.9),
        } };
        const result = xfm.transform(p);
        try std.testing.expectApproxEqAbs(3.157530680769389, result.v[0].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(1.2817556168551425, result.v[1].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.13255153229667402, result.v[2].x, TEST_TOLERANCE);
    }

    {
        const p: ThreeVector(Dual) = .{ .v = [3]Dual{
            .promote(0.8216545960985236),
            .promote(-0.5369856489673557),
            .promote(-0.19112858894806875),
        } };

        const cart = xfm.transform(p);
        try std.testing.expectApproxEqAbs(1.0, cart.v[0].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(1.7631081301296758, cart.v[1].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.5788610557941714, cart.v[2].x, TEST_TOLERANCE);

        const result = xfm.apply(p);
        try std.testing.expectApproxEqAbs(1.0, result.v[0].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.2, result.v[1].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.8, result.v[2].x, TEST_TOLERANCE);
    }
}

test "cartesian to spherical to cartesian" {
    const Dual = ad.DualNumber(f64, 0);

    const x: FourVector(Dual) = .{
        .t = .zero,
        .r = .promote(13.0),
        .th = .promote(std.math.degreesToRadians(20)),
        .ph = .zero,
    };

    const cartFromSpher = cartesianFromSpherical(Dual, x.th, x.ph);
    const spherFromCart = sphericalFromCartesian(Dual, x.th, x.ph);

    var v: ThreeVector(Dual) = .init(
        .promote(1.0),
        .promote(0.5),
        .promote(0.1),
    );
    v = v.toUnitVector();

    // Sanity check:
    const v_len = @sqrt(v.dot(v).x);
    try std.testing.expectApproxEqAbs(1, v_len, TEST_TOLERANCE);

    // First test the regular transformation formulae:
    const t_spher = spherFromCart.transform(v);
    const t_cart = cartFromSpher.transform(t_spher);

    // Now should have the same vector back.
    try std.testing.expectApproxEqAbs(v.v[0].x, t_cart.v[0].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(v.v[1].x, t_cart.v[1].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(v.v[2].x, t_cart.v[2].x, TEST_TOLERANCE);

    // Now with the Jacobian matrix that accounts for changes in the unit
    // vectors:
    const v_spher = spherFromCart.apply(v);

    const v_cart = cartFromSpher.apply(v_spher);

    const len = @sqrt(v_cart.dot(v_cart).x);
    try std.testing.expectApproxEqAbs(1, len, TEST_TOLERANCE);

    // Now should have the same vector back:
    try std.testing.expectApproxEqAbs(v.v[0].x, v_cart.v[0].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(v.v[1].x, v_cart.v[1].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(v.v[2].x, v_cart.v[2].x, TEST_TOLERANCE);
}

/// Get a transformation which maps cartesian (tangent) vectors at some point
/// `(r, theta, phi)` in spherical polar coordinates to a (tangent) spherical
/// vector.
pub fn sphericalFromCartesian(comptime T: type, th: T, ph: T) SphericalFromCartesian(T) {
    return .init(.one, th, ph);
}

pub fn KerrMetric(comptime T: type) type {
    return struct {
        const Self = @This();

        pub const MetricComponents = GeneralMetricComponents(T, .regular);
        pub const InverseMetricComponents = GeneralMetricComponents(T, .inverse);

        M: T,
        a: T,
        /// The positive horizon root
        horizon_radius: T,
        /// The negative horizon root.
        horizon_radius_negative: T,
        /// The innermost stable circular orbit radius.
        isco: T,

        pub fn init(M: T, a: T) Self {
            return .{
                .M = M,
                .a = a,
                .horizon_radius = horizonRadius(M, a, .positive),
                .horizon_radius_negative = horizonRadius(M, a, .negative),
                .isco = innermostStableCircularOrbit(M, a),
            };
        }

        fn iscoZ1(M: T, a: T) T {
            const A = T.Algebra;
            const a_over_M = A.div(a, M);

            const bracket = A.add(
                A.cuberoot(A.add(.one, a_over_M)),
                A.cuberoot(A.sub(.one, a_over_M)),
            );

            const prefactor = A.cuberoot(A.sub(.one, A.powi(a_over_M, 2)));

            return A.add(.one, A.mult(prefactor, bracket));
        }

        fn iscoZ2(M: T, a: T, z1: T) T {
            const A = T.Algebra;
            const term_1 = A.mult(.promote(3), A.powi(A.div(a, M), 2));
            const term_2 = A.powi(z1, 2);
            return A.sqrt(A.add(term_1, term_2));
        }

        fn innermostStableCircularOrbit(M: T, a: T) T {
            const A = T.Algebra;
            const z1 = iscoZ1(M, a);
            const z2 = iscoZ2(M, a, z1);

            const discr = A.mult(
                A.sub(.promote(3), z1),
                A.add(A.add(.promote(3), z1), A.mult(.promote(2), z2)),
            );

            const prefactor = A.mult(M, A.add(.promote(3), z2));

            if (a.x > 0) {
                return A.sub(prefactor, A.sqrt(discr));
            } else {
                return A.add(prefactor, A.sqrt(discr));
            }
        }

        /// Compute Carter's constant at a particular point in the spacetime for a
        /// given four-momentum vector.
        pub fn carterConstant(m: Self, x: FourVector(T), p: FourVector(T)) T {
            const A = T.Algebra;
            return A.sub(
                A.powi(p.th, 2),
                A.mult(
                    A.powi(A.cos(x.th), 2),
                    A.sub(
                        A.powi(A.mult(m.a, p.t), 2),
                        A.powi(A.mult(A.csc(x.th), p.ph), 2),
                    ),
                ),
            );
        }

        /// Calculate the event horizon radius for the spacetime. The `root`
        /// parameter selects either the positive or the negative root.
        fn horizonRadius(
            M: T,
            a: T,
            comptime root: enum { positive, negative },
        ) T {
            const A = T.Algebra;
            const det = A.sqrt(A.sub(A.powi(M, 2), A.powi(a, 2)));
            if (root == .negative) {
                return A.sub(M, det);
            } else {
                return A.add(M, det);
            }
        }

        /// Compute the Delta radial term.
        /// Bardeen et al. 1972, Equation (2.3)
        pub fn delta(self: Self, r: T) T {
            const A = T.Algebra;
            return A.add(
                A.add(
                    A.powi(r, 2),
                    A.mult(A.mult(.promote(-2), self.M), r),
                ),
                A.powi(self.a, 2),
            );
        }

        /// Compute the Sigma term.
        /// Bardeen et al. 1972, Equation (2.3)
        pub fn sigma(self: Self, r: T, theta: T) T {
            const A = T.Algebra;
            return A.add(
                A.powi(r, 2),
                A.powi(A.mult(self.a, A.cos(theta)), 2),
            );
        }

        /// Compute the A term.
        /// Bardeen et al. 1972, Equation (2.3)
        pub fn kerr_a(self: Self, r: T, theta: T) T {
            const A = T.Algebra;

            const a_squared = A.powi(self.a, 2);
            const term_1 = A.powi(A.add(A.powi(r, 2), a_squared), 2);
            const delta_term = self.delta(r);

            return A.sub(
                term_1,
                A.mult(
                    A.mult(a_squared, delta_term),
                    A.powi(A.sin(theta), 2),
                ),
            );
        }

        /// Represent a tangent space at a particular point `r` and `theta` in
        /// the spacetime.
        pub const TangentSpace = struct {
            x: FourVector(T),
            metric_components: MetricComponents,

            /// Adapt the number type.
            pub fn adapt(
                self: TangentSpace,
                comptime NewT: type,
            ) KerrMetric(NewT).TangentSpace {
                return .{
                    .x = self.x.adapt(NewT),
                    .metric_components = self.metric_components.adapt(NewT),
                };
            }

            /// Calculate a local (tetrad) frame local to the point `x_src`
            /// with velocity `v_src`, conventionally represented with
            /// e^mu_(nu). This can be used to map local vectors to the global
            /// coordinates. It is the inverse of `localBasis`.
            pub fn localFrame(
                self: TangentSpace,
                v_src: FourVector(T),
            ) TetradFrame(T) {
                // TODO: this should also be tested for when v_src has all
                // components non-zero, as I worry it could not be solving
                // correctly.
                const basis = orthonormalBasis(T, self, v_src);
                return .{ .m = .{ .v = @bitCast(basis) } };
            }

            /// Calculate the local basis at the point `x_src` moving with
            /// velocity `v_src`, conventionally represented with e^(mu)_nu.
            /// This can be used to project global vectors to the local
            /// coordinates. It is the inverse of `localFrame`.
            pub fn localBasis(
                self: TangentSpace,
                v_src: FourVector(T),
            ) TetradFrame(T) {
                const frame = self.localFrame(v_src);
                // multiply by g_(mu, nu) and transpose all in one go
                var basis = multLeftWithMatrixAndTranspose(T, self, frame.m);
                // the last column is negated here but match what is needed for the
                // inverse tetrad frame
                for (0..4) |row| {
                    const index = 3 + 4 * row;
                    basis.v[index] = basis.v[index].neg();
                }
                return .{ .m = basis };
            }

            /// Constrain a vector to have a given magnitude at a point `x` in the
            /// spacetime.
            /// This currently only gives the positive root, but the negative
            /// (backwards in time) root is also physical.
            pub fn constrainVector(
                self: TangentSpace,
                v: FourVector(T),
                magnitude: T.T,
            ) FourVector(T) {
                return .{
                    .t = self.constrainTime(v, magnitude),
                    .r = v.r,
                    .th = v.th,
                    .ph = v.ph,
                };
            }

            fn constrainTime(
                self: TangentSpace,
                v: FourVector(T),
                magnitude: T.T,
            ) T {
                const A = T.Algebra;
                const m = self.metric_components;

                const term_1 = A.mult(A.mult(m.tt, m.rr), A.powi(v.r, 2));
                const term_2 = A.mult(A.mult(m.tt, m.thth), A.powi(v.th, 2));
                const term_3 = A.mult(m.tt, .promote(magnitude * magnitude));
                const term_4 = A.mult(
                    A.sub(A.mult(m.tt, m.phph), A.powi(m.tph, 2)),
                    A.powi(v.ph, 2),
                );

                const discr = A.add(A.add(term_1, term_2), A.add(term_3, term_4)).neg();

                return A.div(A.add(A.mult(m.tph, v.ph), A.sqrt(discr)), m.tt).neg();
            }

            /// Lower the indices of a four-vector. This is only valid in the tangent
            /// space of the point at `x`.
            pub fn lowerIndices(
                self: TangentSpace,
                v: FourVector(T),
            ) FourVector(T) {
                const m = self.metric_components;
                const A = T.Algebra;
                return .{
                    .t = A.add(A.mult(m.tt, v.t), A.mult(m.tph, v.ph)),
                    .r = A.mult(m.rr, v.r),
                    .th = A.mult(m.thth, v.th),
                    .ph = A.add(A.mult(m.phph, v.ph), A.mult(m.tph, v.t)),
                };
            }

            /// Calculate the locally non-rotating frame for this
            /// tangent space.
            ///
            /// Frame has the latin indices down, i.e. e^mu_(i).
            pub fn lnrFrame(self: TangentSpace) TetradFrame(T) {
                const A = T.Algebra;
                var frame: [4]FourVector(T) = undefined;

                const g_inv = self.metric_components.inverse();
                const omega = A.div(
                    self.metric_components.tph.neg(),
                    self.metric_components.phph,
                );

                // Baker PhD thesis Equation (55)
                frame[0] = .{
                    .t = A.sqrt(g_inv.tt.neg()),
                    .r = .zero,
                    .th = .zero,
                    .ph = A.mult(omega, A.sqrt(g_inv.tt.neg())),
                };
                frame[1] = .{
                    .t = .zero,
                    .r = A.sqrt(g_inv.rr),
                    .th = .zero,
                    .ph = .zero,
                };
                frame[2] = .{
                    .t = .zero,
                    .r = .zero,
                    .th = A.sqrt(g_inv.thth),
                    .ph = .zero,
                };
                frame[3] = .{
                    .t = .zero,
                    .r = .zero,
                    .th = .zero,
                    .ph = A.sqrt(
                        A.add(g_inv.phph, A.mult(A.powi(omega, 2), g_inv.tt.neg())),
                    ),
                };

                return .{ .m = .{ .v = @bitCast(frame) } };
            }

            /// Calculate the locally non-rotating frame for this
            /// tangent space.
            ///
            /// Basis has the latin indices up, i.e. e^(i)_mu.
            pub fn lnrBasis(self: TangentSpace) TetradFrame(T) {
                const A = T.Algebra;
                var frame: [4]FourVector(T) = undefined;

                const g = self.metric_components;
                const omega = A.div(
                    g.tph.neg(),
                    g.phph,
                );

                // Baker PhD thesis Equation (55)
                // Note: this is implicitly transposed as well, to match the
                // rest of the implementation for the bases.
                frame[0] = .{
                    .t = A.sqrt(A.add(g.tt, A.mult(omega, g.tph)).neg()),
                    .r = .zero,
                    .th = .zero,
                    .ph = A.mult(omega.neg(), A.sqrt(g.phph)),
                };
                frame[1] = .{
                    .t = .zero,
                    .r = A.sqrt(g.rr),
                    .th = .zero,
                    .ph = .zero,
                };
                frame[2] = .{
                    .t = .zero,
                    .r = .zero,
                    .th = A.sqrt(g.thth),
                    .ph = .zero,
                };
                frame[3] = .{
                    .t = .zero,
                    .r = .zero,
                    .th = .zero,
                    .ph = A.sqrt(g.phph),
                };

                return .{ .m = .{ .v = @bitCast(frame) } };
            }

            /// Compute the co-rotation velocity vector of the non-rotating
            /// frame.
            pub fn lnrVelocity(self: TangentSpace) FourVector(T) {
                const A = T.Algebra;
                const v = A.div(
                    self.metric_components.tph.neg(),
                    self.metric_components.phph,
                );
                return .{ .t = .one, .r = .zero, .th = .zero, .ph = v };
            }
        };

        /// Obtain a tangent space at a particular point in the spacetime
        /// specified by `r` and `theta`. Sets the `t` and `phi` components to
        /// zero.
        pub fn tangentSpaceAlt(self: Self, r: T, theta: T) TangentSpace {
            return self.tangentSpace(.{ .t = .zero, .r = r, .th = theta, .ph = .zero });
        }

        /// Obtain a tangent space at a particular point in the spacetime
        /// specified by a four vector.
        pub fn tangentSpace(self: Self, x: FourVector(T)) TangentSpace {
            return .{
                .x = x,
                .metric_components = MetricComponents.initPoint(self, x.r, x.th),
            };
        }

        /// Adapt the number type.
        pub fn adapt(self: Self, comptime NewT: type) KerrMetric(NewT) {
            return .{
                .M = .adaptFrom(self.M),
                .a = .adaptFrom(self.a),
                .isco = .adaptFrom(self.isco),
                .horizon_radius = .adaptFrom(self.horizon_radius),
                .horizon_radius_negative = .adaptFrom(self.horizon_radius_negative),
            };
        }

        /// Calculate the radius of the photon orbit.
        pub fn photonOrbit(self: *const Self) T {
            const A = T.Algebra;
            const _angle = A.mult(
                .promote(2.0 / 3.0),
                A.acos(A.div(self.a, self.M)),
            );
            return A.mult(
                A.mult(.promote(2), self.M),
                A.add(.one, A.cos(_angle)),
            );
        }

        /// Write the spacetime information to a HDU.
        pub fn addToHdu(self: *const Self, hdu: *zfits.Hdu) !void {
            try hdu.appendHeaderRecord(
                "MASS",
                .{
                    .value = .{ .float = @floatCast(self.M.x) },
                    .comment = "The mass scale of the black hole",
                },
            );
            try hdu.appendHeaderRecord(
                "SPIN",
                .{
                    .value = .{ .float = @floatCast(self.a.x) },
                    .comment = "The dimensionless spin of black hole",
                },
            );
        }
    };
}

test "spacetime constraints" {
    const Dual = ad.DualNumber(f64, 0);
    const FV = FourVector(Dual);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));

    const x: FV = .{
        .t = .zero,
        .r = .promote(10.0),
        .th = .promote(std.math.degreesToRadians(20)),
        .ph = .zero,
    };

    const ts = metric.tangentSpace(x);

    {
        const v: FV = .{
            .t = .promote(0.3),
            .r = .promote(0.1),
            .th = .zero,
            .ph = .promote(-0.3),
        };
        const null_t = ts.constrainTime(v, 0);
        try std.testing.expectApproxEqAbs(1.167163878033404, null_t.x, TEST_TOLERANCE);

        const null_v = ts.constrainVector(v, 0);
        try std.testing.expectApproxEqAbs(1.167163878033404, null_v.t.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(v.r.x, null_v.r.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(v.th.x, null_v.th.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(v.ph.x, null_v.ph.x, TEST_TOLERANCE);
    }

    {
        const v: FourVector(Dual) = .{
            .t = .one,
            .r = .promote(0.743),
            .th = .promote(0.056),
            .ph = .promote(-0.010),
        };
        const null_t = ts.constrainTime(v, 0);
        try std.testing.expectApproxEqAbs(1.1199560102989092, null_t.x, TEST_TOLERANCE);
    }
}

test "indices" {
    const Dual = ad.DualNumber(f64, 0);
    const FV = FourVector(Dual);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));

    const x: FV = .{
        .t = .zero,
        .r = .promote(10.0),
        .th = .promote(std.math.degreesToRadians(20)),
        .ph = .zero,
    };

    const v: FV = .{
        .t = .promote(0.3),
        .r = .promote(0.1),
        .th = .zero,
        .ph = .promote(-0.3),
    };

    const ts = metric.tangentSpace(x);
    const lowered = ts.lowerIndices(v);
    try std.testing.expectApproxEqAbs(-0.2335795344139173, lowered.t.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(0.12454872917528392, lowered.r.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(0, lowered.th.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-3.5520406328569174, lowered.ph.x, TEST_TOLERANCE);
}

test "lnrf" {
    const Dual = ad.DualNumber(f64, 0);
    const FV = FourVector(Dual);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    const x: FV = .{
        .t = .zero,
        .r = .promote(5.0),
        .th = .promote(std.math.degreesToRadians(70)),
        .ph = .zero,
    };
    const ts = metric.tangentSpace(x);
    const f1 = ts.lnrFrame();

    try std.testing.expectApproxEqAbs(1.2833732687063686, f1.m.v[0].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(0.0, f1.m.v[1].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(0.0, f1.m.v[2].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(0.0193556, f1.m.v[3].x, TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(0.207327, f1.m.v[15].x, TEST_TOLERANCE);

    // Check that the numerical way of calculating this also works
    const v_lnr = ts.lnrVelocity();
    const f2 = ts.localFrame(v_lnr);

    for (f1.m.v, f2.m.v) |v1, v2| {
        try std.testing.expectApproxEqAbs(v1.x, v2.x, TEST_TOLERANCE);
    }

    // And the same for the basis
    const b1 = ts.lnrBasis();
    const b2 = ts.localBasis(v_lnr);

    for (b1.m.v, b2.m.v) |w1, w2| {
        // TODO: remedy sign problems
        try std.testing.expectApproxEqAbs(@abs(w1.x), @abs(w2.x), TEST_TOLERANCE);
    }

    // Check that very close to the black hole doesn't panic:
    const x2: FV = .{
        .t = .zero,
        .r = .promote(1.07),
        .th = .promote(std.math.degreesToRadians(90)),
        .ph = .zero,
    };
    const ts2 = metric.tangentSpace(x2);
    const frame = ts2.localBasis(orbits.keplerianPlunging(Dual, metric, ts2));
    try std.testing.expectApproxEqAbs(-0.6790058343832897, frame.m.v[0].x, TEST_TOLERANCE);
}

/// Container for the static, axis-symmetric metric components.
///
/// The `form` comptime parameter is used to denote whether these are the
/// inverse or regular metric components.
pub fn GeneralMetricComponents(
    comptime T: type,
    comptime form: enum { regular, inverse },
) type {
    return struct {
        pub const component_form = form;
        pub const InverseMetricComponents = GeneralMetricComponents(
            T,
            if (form == .regular) .inverse else .regular,
        );

        const Self = @This();

        pub fn toMatrix(self: *const Self) SMatrix(T, 4, 4) {
            var values: [16]T = undefined;
            for (&values) |*v| v.* = .zero;

            // Assign the non-zero values:
            values[0] = self.tt;
            values[5] = self.rr;
            values[10] = self.thth;
            values[15] = self.phph;

            // And the off-diagonal
            values[3] = self.tph;
            values[12] = self.tph;

            return .{ .v = values };
        }

        /// adapt the number type.
        pub fn adapt(
            self: Self,
            comptime NewT: type,
        ) GeneralMetricComponents(NewT, form) {
            return .{
                .tt = .adaptFrom(self.tt),
                .rr = .adaptFrom(self.rr),
                .thth = .adaptFrom(self.thth),
                .phph = .adaptFrom(self.phph),
                .tph = .adaptFrom(self.tph),
            };
        }

        tt: T,
        rr: T,
        thth: T,
        phph: T,
        tph: T,

        /// Calculate the inverse metric components.
        pub fn inverse(self: Self) InverseMetricComponents {
            const A = T.Algebra;
            // This formula is easy to calculate: since the only non-zero
            // off-diagonal is g_(t ph), the rr and thth case are simple
            // reciprocals. For the g_(t ph) case, simply invert the
            // residual 2x2 matrix.
            const determinant = A.sub(
                A.mult(self.tt, self.phph),
                A.powi(self.tph, 2),
            );

            const tt_inv = A.div(self.phph, determinant);
            const tph_inv = A.div(self.tph.neg(), determinant);
            const phph_inv = A.div(self.tt, determinant);
            const rr_inv = A.div(.one, self.rr);
            const thth_inv = A.div(.one, self.thth);

            return .{
                .tt = tt_inv,
                .rr = rr_inv,
                .thth = thth_inv,
                .phph = phph_inv,
                .tph = tph_inv,
            };
        }

        /// Compute the metric components as (t^2, r^2, theta^2, phi^2,
        /// tphi). This is not a public method and the metric components
        /// should be constructed via `tangentSpace` instead.
        fn initPoint(metric: KerrMetric(T), r: T, theta: T) Self {
            const A = T.Algebra;
            const sin_theta_squared = A.powi(A.sin(theta), 2);

            const R: T = A.mult(.promote(2), metric.M);

            // TODO: instead of using `sigma` here, it's faster to compute it all
            // in one go reusing the trig results.
            const _Sigma = metric.sigma(r, theta);
            const _Sigma_inv = A.div(.one, _Sigma);
            const gamma = A.mult(A.mult(A.mult(sin_theta_squared, R), r), metric.a);

            const tt = A.sub(.one, A.mult(R, A.mult(r, _Sigma_inv))).neg();
            const rr = A.div(_Sigma, metric.delta(r));
            const thth = _Sigma;
            const phph = A.mult(
                sin_theta_squared,
                A.add(
                    A.powi(r, 2),
                    A.add(
                        A.powi(metric.a, 2),
                        A.mult(gamma, A.mult(metric.a, _Sigma_inv)),
                    ),
                ),
            );
            const tph = A.mult(gamma, _Sigma_inv).neg();

            return .{
                .tt = tt,
                .rr = rr,
                .thth = thth,
                .phph = phph,
                .tph = tph,
            };
        }
    };
}

test "inverse metric components" {
    const Dual = ad.DualNumber(f64, 0);
    const MC = GeneralMetricComponents(Dual, .regular);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    const r = Dual.promote(5.0);
    const th = Dual.promote(std.math.degreesToRadians(70));

    const metric_components = MC.initPoint(metric, r, th);

    try std.testing.expectApproxEqAbs(
        -0.6018555178833805,
        metric_components.tt.x,
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        1.570174046920585,
        metric_components.rr.x,
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        25.116510335237862,
        metric_components.thth.x,
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        23.264253165342776,
        metric_components.phph.x,
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        -0.35086728425006897,
        metric_components.tph.x,
        TEST_TOLERANCE,
    );

    const inv_metric_components = metric_components.inverse();

    try std.testing.expectApproxEqAbs(
        -1.647046946830069,
        inv_metric_components.tt.x,
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        0.6368720728515375,
        inv_metric_components.rr.x,
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        0.039814448211661954,
        inv_metric_components.thth.x,
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        0.042609762115183136,
        inv_metric_components.phph.x,
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        -0.02484046597840225,
        inv_metric_components.tph.x,
        TEST_TOLERANCE,
    );
}

/// A tetrad transformation that takes local vectors and maps them to global
/// vectors.
pub fn TetradFrame(comptime T: type) type {
    return struct {
        const Self = @This();
        m: SMatrix(T, 4, 4),

        /// Map a vector to the global coordinates.
        pub fn apply(self: Self, vec: FourVector(T)) FourVector(T) {
            return .fromArray(self.m.multR_T(vec.toArray()));
        }
    };
}

test "kerr metric" {
    const Dual = ad.DualNumber(f64, 0);
    const geom: KerrMetric(Dual) = .init(.one, .promote(0.998));
    try std.testing.expectApproxEqAbs(
        10119.520480000001,
        geom.kerr_a(.promote(10.0), .promote(std.math.pi / 2.0)).x,
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        4.498002,
        geom.sigma(.promote(2.0), .promote(std.math.pi / 4.0)).x,
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        0.996004,
        geom.delta(.promote(2.0)).x,
        TEST_TOLERANCE,
    );

    // metric components
    const comp = KerrMetric(Dual).MetricComponents.initPoint(
        geom,
        .promote(2.0),
        .promote(std.math.degreesToRadians(30)),
    );

    try std.testing.expectApproxEqAbs(-0.157363077293189, comp.tt.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(4.76604812832077, comp.rr.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(4.747003, comp.thth.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(1.3014553590977294, comp.phph.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-0.2102379122153493, comp.tph.x, TEST_TOLERANCE);
}

test "isco" {
    const Dual = ad.DualNumber(f64, 0);
    const M = KerrMetric(Dual);
    try std.testing.expectApproxEqAbs(6.0, M.init(.one, .zero).isco.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(
        1.2369706551751847,
        M.init(.one, .promote(0.998)).isco.x,
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        8.99437445480357,
        M.init(.one, .promote(-0.998)).isco.x,
        TEST_TOLERANCE,
    );
}

test "tetrad frame" {
    const Dual = ad.DualNumber(f64, 0);
    const FV = FourVector(Dual);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));

    const x: FV = .{
        .t = .zero,
        .r = .promote(10.0),
        .th = .promote(std.math.degreesToRadians(20)),
        .ph = .zero,
    };

    const ts = metric.tangentSpace(x);
    {
        const v: FV = .{
            .t = .promote(0.3),
            .r = .promote(0.1),
            .th = .zero,
            .ph = .promote(-0.3),
        };
        const local_frame = ts.localFrame(v);
        const v_test: FV = .{
            .t = .promote(0.1),
            .r = .promote(-0.2),
            .th = .promote(0.2),
            .ph = .zero,
        };
        const result = local_frame.apply(v_test);
        try std.testing.expectApproxEqAbs(-0.298533827094303, result.t.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.30685743929617215, result.r.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.01991262654596481, result.th.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.019393408899399094, result.ph.x, TEST_TOLERANCE);
    }

    // Using a stationary case calculated with Gradus.jl
    // This checks that the ordering of the matrix is correct form the permutations.
    {
        const v: FV = .{
            .t = .promote(1.1168175629588084),
            .r = .zero,
            .th = .zero,
            .ph = .zero,
        };
        const local_frame = ts.localFrame(v);

        try std.testing.expectApproxEqAbs(1.1168175629588084, local_frame.m.v[0].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.0, local_frame.m.v[1].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.0, local_frame.m.v[2].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.0, local_frame.m.v[3].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.0, local_frame.m.v[4].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.89604609128034, local_frame.m.v[5].x, TEST_TOLERANCE);
    }
}

test "tetrad basis" {
    const Dual = ad.DualNumber(f64, 0);
    const FV = FourVector(Dual);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));

    const x: FV = .{
        .t = .zero,
        .r = .promote(10.0),
        .th = .promote(std.math.degreesToRadians(20)),
        .ph = .zero,
    };

    const ts = metric.tangentSpace(x);
    {
        const v: FV = .{
            .t = .promote(0.3),
            .r = .promote(0.1),
            .th = .zero,
            .ph = .promote(-0.3),
        };

        const basis = ts.localBasis(v);

        // row 1, i.e. t components
        try std.testing.expectApproxEqAbs(-0.2326515695280742, basis.m.v[0].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-1.315308110718852, basis.m.v[1].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-7.649336575637465e-17, basis.m.v[2].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-1.6080751944318838, basis.m.v[3].x, TEST_TOLERANCE);

        // row 2, i.e. r components
        try std.testing.expectApproxEqAbs(0.12405392192454977, basis.m.v[0 + 4].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(1.9729621660782786, basis.m.v[1 + 4].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.0, basis.m.v[2 + 4].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(1.6317113079734158, basis.m.v[3 + 4].x, TEST_TOLERANCE);

        // row 3, i.e. theta components
        try std.testing.expectApproxEqAbs(0.0, basis.m.v[0 + 2 * 4].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.0, basis.m.v[1 + 2 * 4].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(10.043878417462158, basis.m.v[2 + 2 * 4].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.0, basis.m.v[3 + 2 * 4].x, TEST_TOLERANCE);

        // row 4, i.e. phi components
        try std.testing.expectApproxEqAbs(-3.537929084134778, basis.m.v[0 + 3 * 4].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.6576540553594256, basis.m.v[1 + 3 * 4].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-3.824668287818732e-17, basis.m.v[2 + 3 * 4].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-1.064171425107412, basis.m.v[3 + 3 * 4].x, TEST_TOLERANCE);
    }
}

/// Project a vector onto a basis in the KerrMetric.
pub fn projectOntoBasis(
    comptime T: type,
    ts: KerrMetric(T).TangentSpace,
    basis: []const FourVector(T),
    v: FourVector(T),
) FourVector(T) {
    const A = T.Algebra;

    var out: FourVector(T) = .zeros;
    for (basis) |basis_vector| {
        // Calculate the projection fraction of the vector onto the basis vector.
        const proj = A.div(
            v.dot(ts, basis_vector),
            basis_vector.properNorm(ts),
        );
        out = out.add(basis_vector.scalarMult(proj));
    }

    return out;
}

/// Perform an iteration of the Gram-Schmidt orthonormalisation routine of the
/// vector `v` onto the basis `basis`.
fn orthonormaliseVector(
    comptime T: type,
    ts: KerrMetric(T).TangentSpace,
    basis: []const FourVector(T),
    v: FourVector(T),
) FourVector(T) {
    const error_tolerance = @sqrt(std.math.floatEps(T.T));
    const A = T.Algebra;

    var new = v;
    var p = projectOntoBasis(T, ts, basis, new);

    // TODO: this could be optimised by calculating all the basis norms once
    // outside of the loop.
    var itt_count: usize = 0;
    while (@abs(p.sum().x) > error_tolerance) {
        new = new.sub(p);
        p = projectOntoBasis(T, ts, basis, new);
        itt_count += 1;
        if (itt_count > 100) {
            @panic("Orthonormalisation failed. Please report this issue.");
        }
    }

    new = new.sub(p);
    const norm = A.sqrt(A.abs(new.properNorm(ts)));
    return new.scalarMult(A.div(.one, norm));
}

test "orthonormalising" {
    const Dual = ad.DualNumber(f64, 0);
    const FV = FourVector(Dual);

    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    const ts = metric.tangentSpaceAlt(
        .promote(2.0),
        .promote(std.math.degreesToRadians(30.0)),
    );

    const v: FV = .{
        .t = .promote(0.3),
        .r = .promote(0.1),
        .th = .promote(0.6),
        .ph = .promote(-0.3),
    };

    // Check the vector norm is working.
    const proj = v.properNorm(ts);
    try std.testing.expectApproxEqAbs(
        1.897392690844379,
        proj.x,
        TEST_TOLERANCE,
    );

    const b1: FV = .{
        .t = .promote(0.07886598517963955),
        .r = .promote(0.3154639407185582),
        .th = .promote(0.07886598517963955),
        .ph = .promote(0.6309278814371164),
    };

    // Check the dot product is working correctly.
    const dotprod = v.dot(ts, b1);
    try std.testing.expectApproxEqAbs(
        0.09009805123737324,
        dotprod.x,
        TEST_TOLERANCE,
    );

    const b2: FV = .{
        .t = .promote(0.4104217983512219),
        .r = .promote(0.10238781810106898),
        .th = .promote(0.4104217983512219),
        .ph = .promote(-0.30832415556580167),
    };
    const b3: FV = .{
        .t = .promote(0.2117555218137058),
        .r = .promote(-0.31593776792548733),
        .th = .promote(0.2117555218137058),
        .ph = .promote(0.5300612184821896),
    };

    // Check projection is working.
    const projected_v = projectOntoBasis(Dual, ts, &.{ b1, b2, b3 }, v);
    try std.testing.expectApproxEqAbs(0.61259911539648, projected_v.t.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(0.09999999999999995, projected_v.r.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(0.61259911539648, projected_v.th.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-0.24950254350108225, projected_v.ph.x, TEST_TOLERANCE);
    {
        const out = orthonormaliseVector(Dual, ts, &.{ b1, b2, b3 }, v);
        try std.testing.expectApproxEqAbs(-2.3337124193431107, out.t.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.0, out.r.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.09405884605978454, out.th.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.3769893629650683, out.ph.x, TEST_TOLERANCE);
    }
    {
        const out = orthonormaliseVector(Dual, ts, &.{ b1, b2 }, v);
        try std.testing.expectApproxEqAbs(-1.581204756998643, out.t.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.4032936062148778, out.r.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.19568130853202714, out.th.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.37752745459428977, out.ph.x, TEST_TOLERANCE);
    }
}

/// Permute the spatial vector inplace.
/// That is, maps (_, a, b, c) -> (_, c, a, b)
fn tetradPermute(comptime T: type, out: *[4]T) void {
    const o1 = out[1];
    out[1] = out[3];
    out[3] = out[2];
    out[2] = o1;
}

fn stateVec(comptime T: type, state: [4]bool) FourVector(T) {
    return .{
        .t = .promote(@floatFromInt(@intFromBool(state[0]))),
        .r = .promote(@floatFromInt(@intFromBool(state[1]))),
        .th = .promote(@floatFromInt(@intFromBool(state[2]))),
        .ph = .promote(@floatFromInt(@intFromBool(state[3]))),
    };
}

fn updateState(state: *[4]bool) void {
    const previous = state.*;
    tetradPermute(bool, state);
    for (state, previous) |*s, p| {
        s.* |= p;
    }
}

fn searchSortedFirstBool(v: []const bool, x: bool) usize {
    var len = v.len;
    var lo: usize = 0;
    var hi: usize = 0;
    while (len != 0) {
        const half = @divTrunc(len, 2);
        const m = lo + half;
        if (v[m] != x) {
            lo = m + 1;
            len -= half + 1;
        } else {
            hi = m;
            len = half;
        }
    }
    return lo;
}

/// Create an orthonormal basis around a particular vector using the
/// Gram-Schmidt procedure. Returns a set of four FourVector(T)s that form the
/// tetrad basis. See also Secion 1.4 of Baker 2025 and the implementation in
/// Gradus.jl.
pub fn orthonormalBasis(
    comptime T: type,
    m: KerrMetric(T).TangentSpace,
    v: FourVector(T),
) [4]FourVector(T) {
    const A = T.Algebra;

    // The t-like directed basis vector is the initial vector.
    const b1 = v.scalarMult(
        A.div(.one, A.sqrt(A.abs(v.properNorm(m)))),
    );

    var state = b1.scalarNotEq(0);

    // Make sure there is an initial direction, else force it to be in the t
    // and phi directions.
    var count: usize = 0;
    for (state) |s| {
        count += @intFromBool(s);
    }
    if (count == 1) {
        state = [4]bool{ true, false, false, true };
    }

    // Calculate the other vectors
    const b2 = orthonormaliseVector(T, m, &.{b1}, stateVec(T, state));

    updateState(&state);
    const b3 = orthonormaliseVector(T, m, &.{ b1, b2 }, stateVec(T, state));

    updateState(&state);
    const b4 = orthonormaliseVector(T, m, &.{ b1, b2, b3 }, stateVec(T, state));

    var out = [4]FourVector(T){ b1, b2, b3, b4 };

    // Sort the vectors according to whichever has the greatest (r, theta, phi)
    // components, in that order:
    const i_r = indexGreatestComponent(T, out[1..], 1, &.{});
    const i_th = indexGreatestComponent(T, out[1..], 2, &.{i_r});
    const i_ph = indexGreatestComponent(T, out[1..], 3, &.{ i_r, i_th });

    return .{
        out[0],
        ensurePositive(T, out[i_r + 1], 1),
        ensurePositive(T, out[i_th + 1], 2),
        ensurePositive(T, out[i_ph + 1], 3),
    };
}

/// Ensure the vector is positive in `component`. Returns the negative version
/// of the vector if that component is negative.
fn ensurePositive(
    comptime T: type,
    vec: FourVector(T),
    component: usize,
) FourVector(T) {
    if (vec.getIndex(component).x < 0) return vec.scalarMult(.promote(-1));
    return vec;
}

/// Return the vector from `vecs` that has the greatest (magnitude) value of
/// the `component`, where `component` is the index according to
/// `FourVector(T).getIndex`.
///
/// Will skip the vectors given in `skip`.
fn indexGreatestComponent(
    comptime T: type,
    vecs: []const FourVector(T),
    component: usize,
    skip: []const usize,
) usize {
    var index: usize = 0;
    var max_value: T.T = 0;
    outer: for (vecs, 0..) |v, i| {
        for (skip) |s| if (i == s) continue :outer;
        const value = @abs(v.getIndex(component).x);
        if (value > max_value) {
            index = i;
            max_value = value;
        }
    }
    return index;
}

test "orthonormalBasis" {
    const Dual = ad.DualNumber(f64, 0);
    const FV = FourVector(Dual);

    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    const ts = metric.tangentSpaceAlt(
        .promote(3.0),
        .promote(std.math.degreesToRadians(30.0)),
    );

    const v: FV = .{
        .t = .promote(0.1),
        .r = .promote(0.4),
        .th = .promote(0.0),
        .ph = .promote(0.5),
    };

    const basis = orthonormalBasis(Dual, ts, v);

    const b1: FV = .{
        .t = .promote(0.09973122116861993),
        .r = .promote(0.39892488467447973),
        .th = .promote(0.0),
        .ph = .promote(0.4986561058430996),
    };
    const b2: FV = .{
        .t = .promote(1.1716729164092647),
        .r = .promote(0.6402421240876959),
        .th = .promote(0),
        .ph = .promote(-0.3906880184594512),
    };
    const b3: FV = .{
        .t = .promote(0),
        .r = .promote(0),
        .th = .promote(0.32030553989136906),
        .ph = .promote(0),
    };
    const b4: FV = .{
        .t = .promote(-1.9805609515046438),
        .r = .promote(-0.39884640866059123),
        .th = .promote(0.0),
        .ph = .promote(0.12839177228742615),
    };

    try b1.testEqual(basis[0]);
    try b2.testEqual(basis[1]);
    try b3.testEqual(basis[2]);
    try b4.testEqual(basis[3]);
}

test "orthonormalBasis problematic" {
    const Dual = ad.DualNumber(f64, 0);

    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    const ts = metric.tangentSpaceAlt(
        .promote(47.2),
        .promote(1.541125729510993),
    );

    const v = @import("orbits.zig").coRotating(Dual, metric, ts);

    const basis = orthonormalBasis(Dual, ts, v);
    try basis[0].testEqual(v);
    try basis[1].testEqual(.{
        .t = .zero,
        .r = .promote(0.9788124505143562),
        .th = .zero,
        .ph = .zero,
    });
    try basis[2].testEqual(.{
        .t = .zero,
        .r = .zero,
        .th = .promote(0.021186436509943223),
        .ph = .zero,
    });
    try basis[3].testEqual(.{
        .t = .promote(0.15232778422909846),
        .r = .zero,
        .th = .zero,
        .ph = .promote(0.021427865565125003),
    });
}

/// Multiply a matrix `mat` by the metric matrix from the left and transpose
/// T[(g * mat)]. This exploits the sparsity pattern present in the matrix to
/// be slightly more efficient than a raw 4x4 matrix multiplication.
/// A tangent space should be supplied so that the metric components are
/// already calculated and can be reused outside of this calculation.
fn multLeftWithMatrixAndTranspose(
    comptime T: type,
    ts: KerrMetric(T).TangentSpace,
    mat: SMatrix(T, 4, 4),
) SMatrix(T, 4, 4) {
    const A = T.Algebra;
    const N = 4;
    const data = mat.v;
    var output: [16]T = .{T.zero} ** 16;

    const m = ts.metric_components;

    // the loop is over the columns of the input matrix
    for (0..4) |col| {
        // these are the input indices, which refer to the ith column of the
        // jth row.
        const _j1 = N * col;
        const _j2 = _j1 + 1;
        const _j3 = _j1 + 2;
        const _j4 = _j1 + 3;

        // these are the output indices, which are transposed
        const _1i = col;
        const _2i = _1i + 1 * N;
        const _3i = _1i + 2 * N;
        const _4i = _1i + 3 * N;

        output[_1i] = A.add(A.mult(data[_j1], m.tt), A.mult(data[_j4], m.tph));
        output[_2i] = A.mult(data[_j2], m.rr);
        output[_3i] = A.mult(data[_j3], m.thth);
        output[_4i] = A.add(
            A.mult(data[_j1], m.tph),
            A.mult(data[_j4], m.phph),
        );
    }

    return .{ .v = output };
}

test multLeftWithMatrixAndTranspose {
    const Dual = ad.DualNumber(f64, 0);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    const ts = metric.tangentSpaceAlt(
        .promote(10.0),
        .promote(std.math.degreesToRadians(20)),
    );

    const sample_matrix: SMatrix(Dual, 4, 4) = .{
        .v = .{
            // row 1 / vec 1
            .promote(0.1335885770022331),
            .promote(0.04452952566741103),
            .promote(0.08905905133482206),
            .promote(-0.1335885770022331),

            // row 2 / vec 2
            .promote(0.0022901668741456805),
            .promote(0.0875464361591424),
            .promote(0.04491830151664404),
            .promote(0.25805897472913586),

            // row 3
            .promote(-0.16086881948281193),
            .promote(0.8989452450240591),
            .promote(-0.010417824015137904),
            .promote(-0.017074920248290843),

            // row 4
            .promote(-1.1361957919350887),
            .promote(0.12186555293816342),
            .promote(-0.012036697804311287),
            .promote(0.01061899109600168),
        },
    };

    const out = multLeftWithMatrixAndTranspose(Dual, ts, sample_matrix);
    // row 1
    try std.testing.expectApproxEqAbs(-0.10401185873066449, out.v[0].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-0.007808954360861072, out.v[1].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(0.12937075734036266, out.v[2].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(0.9106919849214108, out.v[3].x, TEST_TOLERANCE);

    // row 2
    try std.testing.expectApproxEqAbs(0.055460958326542305, out.v[0 + 4].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(0.1090379736744631, out.v[1 + 4].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(0.15178199748691543, out.v[3 + 4].x, TEST_TOLERANCE);

    // row 3
    try std.testing.expectApproxEqAbs(8.984232004920909, out.v[0 + 2 * 4].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(4.531335513280168, out.v[1 + 2 * 4].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-1.2142559798946768, out.v[3 + 2 * 4].x, TEST_TOLERANCE);

    // row 4
    try std.testing.expectApproxEqAbs(-1.5817068453248901, out.v[0 + 3 * 4].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(3.0494273790507394, out.v[1 + 3 * 4].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(0.15178199748705642, out.v[3 + 3 * 4].x, TEST_TOLERANCE);
}

/// Project a ray starting at `x` with tangent vector `v` (both in spherical
/// coordinates) to the celestial sphere at infinity.
///
/// Assumes flat space.
pub fn projectRayToInfinitySpherical(
    comptime T: type,
    x: ThreeVector(T),
    v: ThreeVector(T),
) AnglePair(T) {
    const A = T.Algebra;
    const xfm = cartesianFromSpherical(T, x.v[1], x.v[2]);
    const d = xfm.apply(v).toUnitVector();

    const theta_inf = A.atan2(d.v[1], d.v[0]).neg();
    const phi_inf = A.acos(d.v[0]);

    return .{
        .phi = phi_inf,
        .theta = theta_inf,
    };
}
