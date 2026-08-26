const std = @import("std");
const root = @import("root.zig");
const ad = @import("zad");
const elliptic_integrals = @import("elliptic-integrals.zig");
const potentials = @import("potentials.zig");
const geodesic = @import("geodesic.zig");
const tracy = @import("tracy.zig");

const KerrMetric = root.KerrMetric;
const RadialRoots = potentials.RadialRoots;
const AngularRoots = potentials.AngularRoots;
const RadialCase = potentials.RadialCase;
const AngularCase = potentials.AngularCase;
const ComplexNumber = @import("complex.zig").ComplexNumber;

const TEST_TOLERANCE = @import("options").test_numerical_tolerance;

/// Contains both the half libration and full value of a particular
/// antiderivative.
///
/// TODO: Is this still needed or can functions just return G_theta_final?
pub fn PartialAngularPotentials(comptime T: type) type {
    return struct {
        G_theta_final: T,
    };
}

/// Cast a dual number from one numerical precision to another.
/// TODO: this should be moved to zad
fn castDual(comptime To: type, x: anytype) ad.DualNumber(To, @TypeOf(x).N) {
    var out: ad.DualNumber(To, @TypeOf(x).N) = undefined;
    // Cast down to f64
    inline for (&out.dx, &x.dx) |*t1, t2| {
        t1.* = @floatCast(t2);
    }
    out.x = @floatCast(x.x);
    return out;
}

fn angularAmplitude_normal(comptime T: type, theta: T, sqrt_u_plus: T) T {
    const A = T.Algebra;
    // G&L Equation (28), rearranged
    const amp = A.asin(A.div(
        A.cos(theta),
        sqrt_u_plus,
    ));
    // Handle 'exactly' face-on case, hopefully without losing derivative
    // information.
    if (std.math.isNan(amp.x)) {
        return .{ .x = std.math.asin(@as(T.T, 1)), .dx = amp.dx };
    } else {
        return amp;
    }
}

fn angularAmplitude_vortical(comptime T: type, theta: T, u_plus: T, u_minus: T) T {
    const A = T.Algebra;
    // G&L Equation (59)
    const amp = A.asin(A.sqrt(
        A.div(
            A.sub(A.powi(A.cos(theta), 2), u_minus),
            A.sub(u_plus, u_minus),
        ),
    ));
    if (std.math.isNan(amp.x)) {
        return .{ .x = std.math.asin(@as(T.T, 1)), .dx = amp.dx };
    } else {
        return amp;
    }
}

/// A cache for all values related to calculating a particular angular case.
pub fn AngularCaseCache(comptime T: type) type {
    return struct {
        const Self = @This();

        /// The G_theta angular potential at the start of the geodesic. This is
        /// calculated once when the initial conditions of the geodesic are
        /// known.
        G_theta_init: T,

        /// The half libration value of G_theta, for advancing between turning
        /// points.
        G_theta_half_libration: T,

        /// Case specific cache to save some more expensive computations.
        case_specific: union(enum) {
            normal: struct {
                // 1 / sqrt(-u_minus a^2)
                prefactor: T,
                sqrt_u_plus: T,
                u_ratio: T,
            },
            vortical: struct {
                // 1 / sqrt(u_minus a^2)
                prefactor: T,
                u_plus: T,
                u_minus: T,
                h: T,
            },
        },

        /// Adapt the dual number type.
        pub fn adapt(self: Self, comptime NewT: type) AngularCaseCache(NewT) {
            return .{
                .G_theta_init = .adaptFrom(self.G_theta_init),
                .G_theta_half_libration = .adaptFrom(self.G_theta_half_libration),
                .case_specific = switch (self.case_specific) {
                    .normal => |cs| .{ .normal = .{
                        .prefactor = .adaptFrom(cs.prefactor),
                        .sqrt_u_plus = .adaptFrom(cs.sqrt_u_plus),
                        .u_ratio = .adaptFrom(cs.u_ratio),
                    } },
                    .vortical => |cs| .{ .vortical = .{
                        .prefactor = .adaptFrom(cs.prefactor),
                        .u_plus = .adaptFrom(cs.u_plus),
                        .u_minus = .adaptFrom(cs.u_minus),
                        .h = .adaptFrom(cs.h),
                    } },
                },
            };
        }
    };
}

/// Compute the angular integral in theta, denoted G_theta, along with its
/// half-libration value, for the initial angle.
pub fn angularCache(
    comptime T: type,
    m: KerrMetric(T),
    roots: AngularRoots(T),
    case: potentials.AngularCase,
    theta_init: T,
) AngularCaseCache(T) {
    var ctx = tracy.trace(@src());
    defer ctx.end();
    return switch (case) {
        .normal => angularCache_normal(T, m, roots, theta_init),
        .vortical => angularCache_vortical(T, m, roots, theta_init),
    };
}

/// Compute the angular integral in theta, denoted G_theta, along with its
/// half-libration value, for the initial angle for the vortical motion case.
///
/// G&L: Equation (56) and (60).
fn angularCache_vortical(
    comptime T: type,
    m: KerrMetric(T),
    roots: AngularRoots(T),
    theta_init: T,
) AngularCaseCache(T) {
    const A = T.Algebra;

    const h = T.promote(std.math.sign(@cos(theta_init.x)));
    const u_plus = roots.u_plus;
    const u_minus = roots.u_minus;

    const u_ratio = A.div(u_plus, u_minus);
    const one_sub_u_ratio = A.sub(.one, u_ratio);

    // Equation (59)
    const amplitude_init = angularAmplitude_vortical(T, theta_init, u_plus, u_minus);
    std.debug.assert(!std.math.isNan(amplitude_init.x));

    // This is a common prefactor used in cal_G_theta at both the start and
    // end.
    const prefactor = A.div(
        .one,
        A.sqrt(A.mult(A.powi(m.a, 2), u_minus)),
    );

    // This is K in the text
    const half = elliptic_integrals.Complete.firstKind(T, one_sub_u_ratio);

    // This is F in the text
    const _G_theta = elliptic_integrals.Incomplete.firstKind(
        T,
        amplitude_init,
        one_sub_u_ratio,
    );

    // Equation (56)
    const G_theta_init = A.mult(A.mult(h, prefactor.neg()), _G_theta);

    // Equation (60)
    const half_libration = A.mult(prefactor, half);

    return .{
        .G_theta_half_libration = half_libration,
        .G_theta_init = G_theta_init,
        .case_specific = .{
            .vortical = .{
                .prefactor = prefactor,
                .u_plus = u_plus,
                .u_minus = u_minus,
                .h = h,
            },
        },
    };
}

/// Compute the angular integral in theta, denoted G_theta, along with its
/// half-libration value, for the initial angle for the normal motion case.
///
/// G&L: Equation (29) and (33).
fn angularCache_normal(
    comptime T: type,
    m: KerrMetric(T),
    roots: AngularRoots(T),
    theta_init: T,
) AngularCaseCache(T) {
    const A = T.Algebra;

    const u_plus = roots.u_plus;
    const u_minus = roots.u_minus;

    const u_ratio = A.div(u_plus, u_minus);

    const sqrt_u_plus = A.sqrt(u_plus);

    // Equation (45)
    const amplitude_init = angularAmplitude_normal(T, theta_init, sqrt_u_plus);
    std.debug.assert(!std.math.isNan(amplitude_init.x));

    // This is a common prefactor used in cal_G_theta at both the start and
    // end.
    const prefactor = A.div(
        .one,
        A.sqrt(A.mult(A.powi(m.a, 2).neg(), u_minus)),
    );

    // This is K in the text
    const half = elliptic_integrals.Complete.firstKind(T, u_ratio);

    // This is F in the text
    const _G_theta = elliptic_integrals.Incomplete.firstKind(
        T,
        amplitude_init,
        u_ratio,
    );

    const G_theta_init = A.mult(prefactor.neg(), _G_theta);

    const half_libration = A.mult(
        A.mult(.promote(2.0), prefactor),
        half,
    );

    return .{
        .G_theta_half_libration = half_libration,
        .G_theta_init = G_theta_init,
        .case_specific = .{
            .normal = .{
                .prefactor = prefactor,
                .sqrt_u_plus = sqrt_u_plus,
                .u_ratio = u_ratio,
            },
        },
    };
}

/// Compute the G_theta antiderivative at `theta_final` given the initial
/// values calculated in `AngularCase`.
pub fn angularFromCache(
    comptime T: type,
    case: potentials.AngularCase,
    cache: AngularCaseCache(T),
    theta_final: T,
) PartialAngularPotentials(T) {
    const A = T.Algebra;
    var parameter: T = .zero;
    var amplitude: T = .zero;
    var prefactor: T = .zero;

    switch (case) {
        .normal => {
            const sp = cache.case_specific.normal;
            amplitude = angularAmplitude_normal(T, theta_final, sp.sqrt_u_plus);
            parameter = sp.u_ratio;
            prefactor = sp.prefactor;
        },
        .vortical => {
            const sp = cache.case_specific.vortical;
            amplitude = angularAmplitude_vortical(T, theta_final, sp.u_plus, sp.u_minus);
            parameter = A.sub(.one, A.div(sp.u_plus, sp.u_minus));
            // TODO: check this is correct sign?
            prefactor = A.mult(sp.h, sp.prefactor);
        },
    }
    std.debug.assert(!std.math.isNan(amplitude.x));

    // G&L Equation (29) for normal
    // G&L Equation (56) for vortical
    const cal_G_theta_final = elliptic_integrals.Incomplete.firstKind(
        T,
        amplitude,
        parameter,
    );

    return .{
        .G_theta_final = A.mult(prefactor.neg(), cal_G_theta_final),
    };
}

pub fn PrincipalAngularValues(comptime T: type) type {
    return struct {
        theta: T,
    };
}

/// Compute the angular values for a particular case of angular roots.  The
/// `theta` angle is the initial inclination of the null-geodesic at the origin
/// of the trajectory. `theta` sign should be the sign of the poloidal
/// component of the four-momentum of the geodesic.
pub fn angularValues(
    comptime T: type,
    case: potentials.AngularCase,
    cache: AngularCaseCache(T),
    mino_time: T,
    sign_theta: T.T,
) PrincipalAngularValues(T) {
    return switch (case) {
        .normal => angularValues_normal(T, cache, mino_time, sign_theta),
        .vortical => angularValues_vortical(T, cache, mino_time, sign_theta),
    };
}

fn angularValues_vortical(
    comptime T: type,
    cache: AngularCaseCache(T),
    mino_time: T,
    sign_theta: T.T,
) PrincipalAngularValues(T) {
    const A = T.Algebra;

    const sign: T = .promote(sign_theta);

    const sp = cache.case_specific.vortical;

    // Equation (66)
    const argument = A.div(
        A.add(mino_time, A.mult(sign, cache.G_theta_init)),
        sp.prefactor,
    );
    const _Upsilon = elliptic_integrals.Jacobi.am(
        T,
        argument,
        A.sub(.one, A.div(sp.u_plus, sp.u_minus)),
    );

    // Equation (69)
    const _sin_Upsilon_squared = A.powi(A.sin(_Upsilon), 2);
    const cos_theta = A.mult(
        sp.h,
        A.sqrt(A.add(
            sp.u_minus,
            A.mult(
                A.sub(sp.u_plus, sp.u_minus),
                _sin_Upsilon_squared,
            ),
        )),
    );

    return .{ .theta = A.acos(cos_theta) };
}

fn angularValues_normal(
    comptime T: type,
    cache: AngularCaseCache(T),
    mino_time: T,
    sign_theta: T.T,
) PrincipalAngularValues(T) {
    const A = T.Algebra;

    const sign: T = .promote(sign_theta);

    const sp = cache.case_specific.normal;

    // Equation (38)
    const term_1 = A.mult(
        A.div(.one, sp.prefactor),
        A.add(mino_time, A.mult(sign, cache.G_theta_init)),
    );

    const cos_theta_sqrt_u_plus = A.mult(
        sign.neg(),
        elliptic_integrals.Jacobi.sn(T, term_1, sp.u_ratio),
    );

    const theta = A.acos(A.mult(cos_theta_sqrt_u_plus, sp.sqrt_u_plus));

    return .{ .theta = theta };
}

/// A cache for all values related to calculating a particular radial case.
pub fn RadialCaseCache(comptime T: type) type {
    return struct {
        const Self = @This();
        /// The radial coordinate at which this cache has been evaluated. This
        /// is usually the starting point of the geodesic.
        r: T,

        /// The values of the `k` parameter used in the Jacobi elliptic
        /// integrals (X|k). It depends only on the value of the radial roots.
        ///
        /// E.g. G&L Equation (B16)
        k: T,

        /// The evaluated incomplete integral of the first kind with `x` as the
        /// amplitude and `k` as the parameter.
        ///
        /// E.g. G&L Equation (B16)
        I_0: T,

        /// The small `x` used in integrating the radial coordinates. This is
        /// calculated from the initial positoin and the radial roots.
        ///
        /// E.g. G&L Equation (B15)
        x: T,

        /// Case specific cache to save some more expensive computations.
        case_specific: union(enum) {
            case_I: struct {
                coeff: T,
                r1: T,
                r2: T,
                r31: T,
                r32: T,
            },
            case_II: struct {
                coeff: T,
                r3: T,
                r4: T,
                r31: T,
                r41: T,
            },
            case_III: struct {
                coeff: T,
                t1: T,
                t2: T,
                d1: T,
                d2: T,
            },
            case_IV: struct {
                /// This is 2 / (C + D)
                coeff: T,
                g0: T,
                a2: T,
                b1: T,
                /// This is arctan{x_4(r_s)} + arctan(g0)
                phi: T,
            },
        },

        /// Adapt the dual number type.
        pub fn adapt(self: Self, comptime NewT: type) RadialCaseCache(NewT) {
            return .{
                .r = .adaptFrom(self.r),
                .k = .adaptFrom(self.k),
                .I_0 = .adaptFrom(self.I_0),
                .x = .adaptFrom(self.x),
                .case_specific = switch (self.case_specific) {
                    .case_I => |cs| .{
                        .case_I = .{
                            .coeff = .adaptFrom(cs.coeff),
                            .r1 = .adaptFrom(cs.r1),
                            .r2 = .adaptFrom(cs.r2),
                            .r31 = .adaptFrom(cs.r31),
                            .r32 = .adaptFrom(cs.r32),
                        },
                    },
                    .case_II => |cs| .{
                        .case_II = .{
                            .coeff = .adaptFrom(cs.coeff),
                            .r3 = .adaptFrom(cs.r3),
                            .r4 = .adaptFrom(cs.r4),
                            .r31 = .adaptFrom(cs.r31),
                            .r41 = .adaptFrom(cs.r41),
                        },
                    },
                    .case_III => |cs| .{
                        .case_III = .{
                            .coeff = .adaptFrom(cs.coeff),
                            .t1 = .adaptFrom(cs.t1),
                            .t2 = .adaptFrom(cs.t2),
                            .d1 = .adaptFrom(cs.d1),
                            .d2 = .adaptFrom(cs.d2),
                        },
                    },
                    .case_IV => |cs| .{
                        .case_IV = .{
                            .coeff = .adaptFrom(cs.coeff),
                            .g0 = .adaptFrom(cs.g0),
                            .a2 = .adaptFrom(cs.a2),
                            .b1 = .adaptFrom(cs.b1),
                            .phi = .adaptFrom(cs.phi),
                        },
                    },
                },
            };
        }
    };
}

pub fn PrincipalRadialValues(comptime T: type) type {
    return struct {
        const Self = @This();

        /// The radial coordinate.
        r: T,

        /// The X potential, e.g. G&L Equation (B26).
        X: T,

        /// The Jacobi integrals (X|k)
        integrals: elliptic_integrals.Jacobi.Integrals(T),

        /// Adapt the dual number type.
        pub fn adapt(self: Self, comptime NewT: type) PrincipalRadialValues(NewT) {
            return .{
                .r = .adaptFrom(self.r),
                .X = .adaptFrom(self.X),
                .integrals = .{
                    .cn = .adaptFrom(self.integrals.cn),
                    .dn = .adaptFrom(self.integrals.dn),
                    .sc = .adaptFrom(self.integrals.sc),
                    .sn = .adaptFrom(self.integrals.sn),
                },
            };
        }
    };
}

/// Compute the radial cache for a particular case of radial potential roots.
pub fn radialCache(
    comptime T: type,
    case: RadialCase,
    roots: RadialRoots(T),
    r: T,
) RadialCaseCache(T) {
    var ctx = tracy.trace(@src());
    defer ctx.end();
    return switch (case) {
        .case_I => radialCache_case_I(T, roots, r),
        .case_II => radialCache_case_II(T, roots, r),
        .case_III => radialCache_case_III(T, roots, r),
        .case_IV => radialCache_case_IV(T, roots, r),
    };
}

/// Compute the radial cache for a particular case of radial potential roots.
pub fn radialValuesFromCache(
    comptime T: type,
    case: RadialCase,
    cache: RadialCaseCache(T),
    mino_time: T,
    r_sign: T.T,
) PrincipalRadialValues(T) {
    return switch (case) {
        .case_I => fromCache_case_I(T, cache, mino_time, r_sign),
        .case_II => fromCache_case_II(T, cache, mino_time, r_sign),
        .case_III => fromCache_case_III(T, cache, mino_time, r_sign),
        .case_IV => fromCache_case_IV(T, cache, mino_time, r_sign),
    };
}

/// Compute the radial values for a particular case of radial potential roots.
/// The `r` radius is the radius of the null-geodesic at the start of the
/// trajectory. `r` sign should be the sign of the radial component of the
/// four-momentum of the geodesic.
pub fn radialValues(
    comptime T: type,
    case: RadialCase,
    roots: RadialRoots(T),
    r: T,
    mino_time: T,
    r_sign: T.T,
) PrincipalRadialValues(T) {
    return radialValuesFromCache(
        T,
        case,
        radialCache(T, case, roots, r),
        mino_time,
        r_sign,
    );
}

fn clampRoundingErrors(comptime T: type, value: T) T {
    const error_tolerance = @sqrt(std.math.floatEps(T.T));
    var x = value;
    // Check for possible rounding errors if less than zero
    if (x.x < 0 and @abs(x.x) < error_tolerance) {
        x.x = 0;
    }
    // Check for possible rounding errors if in excess of 1
    if (@abs(x.x) > 1 and @abs(x.x) < 1 + error_tolerance) {
        x.x = 1;
    }
    return x;
}

/// Compute the cache for determining the radial position as a function of the
/// mino time for a given set of radial roots.
///
/// For this case, all roots must be real. This is G&L Equation (B26), using
/// (B20), (B15) and (B13).
fn radialCache_case_I(
    comptime T: type,
    roots: RadialRoots(T),
    r: T,
) RadialCaseCache(T) {
    const A = T.Algebra;

    const r1 = roots.r1.real();
    const r2 = roots.r2.real();
    const r3 = roots.r3.real();
    const r4 = roots.r4.real();

    const r31 = A.sub(r3, r1);
    const r32 = A.sub(r3, r2);
    const r41 = A.sub(r4, r1);
    const r42 = A.sub(r4, r2);

    const k = A.div(
        A.mult(r32, r41),
        A.mult(r31, r42),
    );

    const x1_squared = A.mult(
        A.div(r31, r32),
        A.div(A.sub(r, r2), A.sub(r, r1)),
    );

    const x1 = A.sqrt(clampRoundingErrors(T, x1_squared));

    std.debug.assert(@abs(x1.x) <= 1);

    const prefactor = A.div(
        .promote(2.0),
        A.sqrt(A.mult(r31, r42)),
    );

    const anti_derivative = A.mult(
        prefactor,
        elliptic_integrals.Incomplete.firstKind(T, A.asin(x1), k),
    );

    const coeff = A.mult(A.sqrt(A.mult(r31, r42)), .promote(1.0 / 2.0));

    return .{
        .r = r,
        .I_0 = anti_derivative,
        .k = k,
        .x = x1,
        .case_specific = .{
            .case_I = .{
                .coeff = coeff,
                .r1 = r1,
                .r2 = r2,
                .r31 = r31,
                .r32 = r32,
            },
        },
    };
}

fn fromCache_case_I(
    comptime T: type,
    cache: RadialCaseCache(T),
    mino_time: T,
    r_sign: T.T,
) PrincipalRadialValues(T) {
    const A = T.Algebra;

    const sp = cache.case_specific.case_I;

    const X_1 = A.mult(sp.coeff, A.add(
        mino_time,
        A.mult(.promote(r_sign), cache.I_0),
    ));

    const integrals = elliptic_integrals.Jacobi.all(T, X_1, cache.k);
    const sn2_X = A.powi(integrals.sn, 2);

    const numerator = A.sub(A.mult(sp.r2, sp.r31), A.mult(A.mult(sp.r1, sp.r32), sn2_X));
    const denominator = A.sub(sp.r31, A.mult(sp.r32, sn2_X));

    return .{
        .r = A.div(numerator, denominator),
        .X = X_1,
        .integrals = integrals,
    };
}

/// Compute the cache for determining the radial position as a function of the
/// mino time for a given set of radial roots.
///
/// For this case, all roots must be real. This is G&L Equation (B46), using
/// (B45), (B40), (B35).
fn radialCache_case_II(
    comptime T: type,
    roots: RadialRoots(T),
    r: T,
) RadialCaseCache(T) {
    const A = T.Algebra;

    const r1 = roots.r1.real();
    const r2 = roots.r2.real();
    const r3 = roots.r3.real();
    const r4 = roots.r4.real();

    const r31 = A.sub(r3, r1);
    const r32 = A.sub(r3, r2);
    const r41 = A.sub(r4, r1);
    const r42 = A.sub(r4, r2);

    const k = A.div(
        A.mult(r32, r41),
        A.mult(r31, r42),
    );

    const x2_squared = A.mult(
        A.div(r31, r41),
        A.div(A.sub(r, r4), A.sub(r, r3)),
    );

    const x2 = A.sqrt(clampRoundingErrors(T, x2_squared));

    std.debug.assert(@abs(x2.x) <= 1);

    const prefactor = A.div(
        .promote(2.0),
        A.sqrt(A.mult(r31, r42)),
    );
    const anti_derivative = A.mult(
        prefactor,
        elliptic_integrals.Incomplete.firstKind(T, A.asin(x2), k),
    );

    const coeff = A.mult(A.sqrt(A.mult(r31, r42)), .promote(1.0 / 2.0));

    return .{
        .r = r,
        .I_0 = anti_derivative,
        .k = k,
        .x = x2,
        .case_specific = .{
            .case_II = .{
                .coeff = coeff,
                .r3 = r3,
                .r4 = r4,
                .r31 = r31,
                .r41 = r41,
            },
        },
    };
}

fn fromCache_case_II(
    comptime T: type,
    cache: RadialCaseCache(T),
    mino_time: T,
    r_sign: T.T,
) PrincipalRadialValues(T) {
    const A = T.Algebra;

    const sp = cache.case_specific.case_II;

    const X_2 = A.mult(sp.coeff, A.add(
        mino_time,
        A.mult(.promote(r_sign), cache.I_0),
    ));

    const integrals = elliptic_integrals.Jacobi.all(T, X_2, cache.k);
    const sn2_X = A.powi(integrals.sn, 2);

    const numerator = A.sub(A.mult(sp.r4, sp.r31), A.mult(A.mult(sp.r3, sp.r41), sn2_X));
    const denominator = A.sub(sp.r31, A.mult(sp.r41, sn2_X));

    return .{
        .r = A.div(numerator, denominator),
        .X = X_2,
        .integrals = integrals,
    };
}

/// A helper routine for computing
///     sqrt((r3 - r1) * (r4 - r1))
/// when r3 and r4 are complex valued
fn complexSubMultSqrt(comptime T: type, r3: ComplexNumber(T), r4: ComplexNumber(T), r1: T) T {
    const a = r3.addRe(r1.neg());
    const b = r4.addRe(r1.neg());
    return complexSquareProduct(T, a, b);
}

/// A helper routine for computing
///     sqrt(a * b)
/// when a and b are complex valued.
fn complexSquareProduct(comptime T: type, a: ComplexNumber(T), b: ComplexNumber(T)) T {
    const A = T.Algebra;
    const product: ComplexNumber(T) = .{
        .mag = A.mult(a.mag, b.mag),
        .arg = A.add(a.arg, b.arg),
    };
    return product.sqrt().real();
}

/// Equation (B55) or Equation (B58)
fn x_case_III(comptime T: type, r: T, _A: T, _B: T, r1: T, r2: T) T {
    const A = T.Algebra;
    const common = A.div(
        A.mult(_B, A.sub(r, r2)),
        A.mult(_A, A.sub(r, r1)),
    );
    return A.div(
        A.sub(.one, common),
        A.add(.one, common),
    );
}

/// Compute the cache for determining the radial position as a function of the
/// mino time for a given set of radial roots.
///
/// Only two roots are real, and r3 and r4 are complex conjugates. G&L Equation
/// (B75).
fn radialCache_case_III(
    comptime T: type,
    roots: RadialRoots(T),
    r: T,
) RadialCaseCache(T) {
    const A = T.Algebra;

    const r1 = roots.r1.real();
    const r2 = roots.r2.real();

    const r21 = A.sub(r2, r1);

    const _A = complexSubMultSqrt(T, roots.r3, roots.r4, r2);
    const _B = complexSubMultSqrt(T, roots.r3, roots.r4, r1);

    const k = A.div(
        A.sub(A.powi(A.add(_A, _B), 2), A.powi(r21, 2)),
        A.mult(.promote(4.0), A.mult(_A, _B)),
    );

    const x3 = x_case_III(T, r, _A, _B, r1, r2);
    std.debug.assert(@abs(x3.x) < 1);

    const coeff = A.sqrt(A.mult(_A, _B));

    const anti_derivative = A.mult(
        A.div(.one, coeff),
        elliptic_integrals.Incomplete.firstKind(T, A.acos(x3), k),
    );

    const t1 = A.sub(A.mult(_B, r2), A.mult(_A, r1));
    const t2 = A.add(A.mult(_B, r2), A.mult(_A, r1));

    const d1 = A.sub(_B, _A);
    const d2 = A.add(_B, _A);

    return .{
        .r = r,
        .I_0 = anti_derivative,
        .k = k,
        .x = x3,
        .case_specific = .{
            .case_III = .{
                .coeff = coeff,
                .t1 = t1,
                .t2 = t2,
                .d1 = d1,
                .d2 = d2,
            },
        },
    };
}

fn fromCache_case_III(
    comptime T: type,
    cache: RadialCaseCache(T),
    mino_time: T,
    r_sign: T.T,
) PrincipalRadialValues(T) {
    const A = T.Algebra;

    const sp = cache.case_specific.case_III;

    const X_3 = A.mult(sp.coeff, A.add(
        mino_time,
        A.mult(.promote(r_sign), cache.I_0),
    ));

    const integrals = elliptic_integrals.Jacobi.all(T, X_3, cache.k);
    const cn_X = integrals.cn;

    const numerator = A.add(sp.t1, A.mult(sp.t2, cn_X));
    const denominator = A.add(sp.d1, A.mult(sp.d2, cn_X));

    return .{
        .r = A.div(numerator, denominator),
        .X = X_3,
        .integrals = integrals,
    };
}

/// Compute the cache for determining the radial position as a function of the
/// mino time for a given set of radial roots.
/// G&L Equation (B109),
fn radialCache_case_IV(
    comptime T: type,
    roots: RadialRoots(T),
    r: T,
) RadialCaseCache(T) {
    const A = T.Algebra;

    const r31_complex = roots.r3.sub(roots.r1);
    const r32_complex = roots.r3.sub(roots.r2);
    const r41_complex = roots.r4.sub(roots.r1);
    const r42_complex = roots.r4.sub(roots.r2);

    const _C = complexSquareProduct(T, r31_complex, r42_complex);
    const _D = complexSquareProduct(T, r32_complex, r41_complex);

    const k_sqrt = A.sqrt(A.div(
        A.powi(_D, 2),
        A.powi(_C, 2),
    ));

    const k4 = A.div(
        A.mult(.promote(4), k_sqrt),
        A.powi(A.add(.one, k_sqrt), 2),
    );

    // Equation (B11) for the expansion of a2, then (B88)
    const a2 = roots.r2.imag();
    std.debug.assert(a2.x > 0);

    const four_a2_squared = A.mult(.promote(4.0), A.powi(a2, 2));
    const g0_numerator = A.sub(
        four_a2_squared,
        A.powi(A.sub(_C, _D), 2),
    );
    const g0_denominator = A.sub(
        A.powi(A.add(_C, _D), 2),
        four_a2_squared,
    );
    const g0 = A.sqrt(A.div(g0_numerator, g0_denominator));
    std.debug.assert(g0.x > 0 and g0.x < 1);

    // Equation (B83)
    const b1 = roots.r3.real();
    const x4 = A.div(A.add(r, b1), a2);

    const coeff = A.div(.promote(2.0), A.add(_C, _D));
    const phi = A.add(A.atan(x4), A.atan(g0));

    const anti_derivative = A.mult(
        coeff,
        elliptic_integrals.Incomplete.firstKind(
            T,
            A.add(A.atan(x4), A.atan(g0)),
            k4,
        ),
    );

    return .{
        .r = r,
        .I_0 = anti_derivative,
        .k = k4,
        .x = x4,
        .case_specific = .{
            .case_IV = .{
                .coeff = coeff,
                .phi = phi,
                .g0 = g0,
                .a2 = a2,
                .b1 = b1,
            },
        },
    };
}

fn fromCache_case_IV(
    comptime T: type,
    cache: RadialCaseCache(T),
    mino_time: T,
    r_sign: T.T,
) PrincipalRadialValues(T) {
    const A = T.Algebra;

    const sp = cache.case_specific.case_IV;

    // Equation (B104)
    const X_4 = A.mult(
        A.div(.one, sp.coeff),
        A.add(A.mult(.promote(r_sign), mino_time), cache.I_0),
    );

    // Equation (B109)
    const integrals = elliptic_integrals.Jacobi.all(T, X_4, cache.k);
    const sc_X_4 = integrals.sc;
    const brackets = A.div(
        A.sub(sp.g0, sc_X_4),
        A.add(.one, A.mult(sp.g0, sc_X_4)),
    );

    const r = A.mult(.promote(-1), A.add(A.mult(sp.a2, brackets), sp.b1));

    return .{
        .X = X_4,
        .r = r,
        .integrals = integrals,
    };
}

fn RadialAntiderivatives(comptime T: type) type {
    return struct {
        const Self = @This();
        // these are always calculated at the source
        // i_0 == I_0 == mino_time by definition
        I_0: T,
        I_1: T,
        I_2: T,
        I_plus: T,
        I_minus: T,
        am_X: T,
    };
}

fn AngularAntiderivatives(comptime T: type) type {
    return struct {
        /// Azimuthal angular anti-derivative at the start point.
        G_phi_init: T,
        /// Azimuthal angular anti-derivative at the end point.
        G_phi_final: T,
        /// Time angular anti-derivative at the start point.
        G_t_init: T,
        /// Time angular anti-derivative at the end point.
        G_t_final: T,
        /// Half-libration antiderivatives. In the text, these are denoted with a hat.
        G_phi_half: T,
        /// Half-libration antiderivatives. In the text, these are denoted with a hat.
        G_t_half: T,
    };
}

pub fn angularAntiderivatives(
    comptime T: type,
    geometry: KerrMetric(T),
    case: AngularCase,
    cache: AngularCaseCache(T),
    roots: AngularRoots(T),
    theta_init: T,
    theta_final: T,
    sign_theta: T.T,
) AngularAntiderivatives(T) {
    return switch (case) {
        .normal => angularAntiderivatives_normal(
            T,
            geometry,
            cache,
            roots,
            theta_init,
            theta_final,
            sign_theta,
        ),
        .vortical => angularAntiderivatives_vortical(
            T,
            geometry,
            cache,
            roots,
            theta_init,
            theta_final,
            sign_theta,
        ),
    };
}

fn angularAntiderivatives_vortical(
    comptime T: type,
    geometry: KerrMetric(T),
    cache: AngularCaseCache(T),
    roots: AngularRoots(T),
    theta_init: T,
    theta_final: T,
    sign_theta: T.T,
) AngularAntiderivatives(T) {
    const A = T.Algebra;

    const sp = cache.case_specific.vortical;

    const sign: T = .promote(sign_theta);

    // TODO: remove this parameter? Or change the case specific cache?
    _ = roots;

    // Equation (59) or (65)
    const upsilon_init = A.mult(
        A.mult(sp.h, sign.neg()),
        angularAmplitude_vortical(T, theta_init, sp.u_plus, sp.u_minus),
    );
    const upsilon_final = A.mult(
        A.mult(sp.h, sign.neg()),
        angularAmplitude_vortical(T, theta_final, sp.u_plus, sp.u_minus),
    );

    const parameter = A.sub(.one, A.div(sp.u_plus, sp.u_minus));

    // Equation (57), for a moment without the coefficient
    const common_argument_phi = A.div(
        A.sub(sp.u_plus, sp.u_minus),
        A.sub(.one, sp.u_minus),
    );
    const cal_G_phi_init = elliptic_integrals.Incomplete.thirdKind(
        T,
        common_argument_phi,
        upsilon_init,
        parameter,
    );
    const cal_G_phi_final = elliptic_integrals.Incomplete.thirdKind(
        T,
        common_argument_phi,
        upsilon_final,
        parameter,
    );

    // Equation (57), now with the coefficient. The prefactor does not include
    // the `h` term so it can be used in Equation (61) later.
    const prefactor_G_phi = A.div(
        sp.prefactor,
        A.sub(.one, sp.u_minus),
    );
    const _G_phi_final = A.mult(A.mult(sp.h.neg(), prefactor_G_phi), cal_G_phi_final);
    const _G_phi_init = A.mult(A.mult(sp.h.neg(), prefactor_G_phi), cal_G_phi_init);

    // Equation (58), only the elliptic integral
    const cal_G_t_init = elliptic_integrals.Incomplete.secondKind_dk(
        T,
        upsilon_init,
        parameter,
    );
    const cal_G_t_final = elliptic_integrals.Incomplete.secondKind_dk(
        T,
        upsilon_final,
        parameter,
    );

    // Equation (58), with coefficient, also missing the `h` term as above.
    const prefactor_G_t = A.sqrt(A.div(
        sp.u_minus,
        A.powi(geometry.a, 2),
    ));
    const _G_t_final = A.mult(A.mult(sp.h.neg(), prefactor_G_t), cal_G_t_final);
    const _G_t_init = A.mult(A.mult(sp.h.neg(), prefactor_G_t), cal_G_t_init);

    // Equation (61), the half libration value
    const G_phi_half = A.mult(
        prefactor_G_phi,
        elliptic_integrals.Complete.thirdKind(T, common_argument_phi, parameter),
    );

    // Equation (62)
    const G_t_half = A.mult(
        prefactor_G_t,
        elliptic_integrals.Complete.secondKind_dk(T, parameter),
    );

    return .{
        // TODO: I think this is the wrong way round, so probably worth
        // renaming the variables?
        .G_phi_init = _G_phi_init,
        .G_phi_final = _G_phi_final,
        .G_t_init = _G_t_init,
        .G_t_final = _G_t_final,
        .G_phi_half = G_phi_half,
        .G_t_half = G_t_half,
    };
}

fn angularAntiderivatives_normal(
    comptime T: type,
    geometry: KerrMetric(T),
    cache: AngularCaseCache(T),
    roots: AngularRoots(T),
    theta_init: T,
    theta_final: T,
    sign_theta: T.T,
) AngularAntiderivatives(T) {
    const A = T.Algebra;

    const sp = cache.case_specific.normal;
    _ = sign_theta;

    // Equation (45)
    const phi_init = angularAmplitude_normal(T, theta_init, sp.sqrt_u_plus).neg();
    const phi_final = angularAmplitude_normal(T, theta_final, sp.sqrt_u_plus).neg();

    // Common terms
    const sqrt_u_a_sq = A.sqrt(A.mult(A.powi(geometry.a, 2), roots.u_minus).neg());

    // Equation (30), without the coefficient, as that is applied in (47) here.
    const cal_G_phi_init = elliptic_integrals.Incomplete.thirdKind(
        T,
        roots.u_plus,
        phi_init,
        sp.u_ratio,
    );
    const cal_G_phi_final = elliptic_integrals.Incomplete.thirdKind(
        T,
        roots.u_plus,
        phi_final,
        sp.u_ratio,
    );
    const prefactor_G_phi = A.div(.promote(-1), sqrt_u_a_sq);

    // Equation (47)
    const _G_phi_final = A.mult(prefactor_G_phi, cal_G_phi_final);
    const _G_phi_init = A.mult(prefactor_G_phi, cal_G_phi_init);

    // Equation (31)
    const cal_G_t_init = elliptic_integrals.Incomplete.secondKind_dk(T, phi_init, sp.u_ratio);
    const cal_G_t_final = elliptic_integrals.Incomplete.secondKind_dk(T, phi_final, sp.u_ratio);
    const prefactor_G_t = A.div(
        A.mult(.promote(-2), roots.u_plus),
        sqrt_u_a_sq,
    );

    // Equation (48)
    const _G_t_final = A.mult(prefactor_G_t, cal_G_t_final);
    const _G_t_init = A.mult(prefactor_G_t, cal_G_t_init);

    // Equation (34)
    const G_phi_half = A.mult(
        .promote(-2),
        A.mult(prefactor_G_phi, elliptic_integrals.Complete.thirdKind(T, roots.u_plus, sp.u_ratio)),
    );

    // Equation (35)
    const G_t_half = A.mult(
        .promote(-2),
        A.mult(prefactor_G_t, elliptic_integrals.Complete.secondKind_dk(T, sp.u_ratio)),
    );

    return .{
        // TODO: I think this is the wrong way round, so probably worth
        // renaming the variables?
        .G_phi_init = _G_phi_init,
        .G_phi_final = _G_phi_final,
        .G_t_init = _G_t_init,
        .G_t_final = _G_t_final,
        .G_phi_half = G_phi_half,
        .G_t_half = G_t_half,
    };
}

fn radialAntiderivatives(
    comptime T: type,
    geometry: KerrMetric(T),
    case: RadialCase,
    cache: RadialCaseCache(T),
    roots: RadialRoots(T),
    rv: PrincipalRadialValues(T),
    r_init: T,
    mino_time: T,
    r_sign: T.T,
) RadialAntiderivatives(T) {
    return switch (case) {
        .case_I => radialAntiderivatives_case_I_and_II(
            T,
            geometry,
            cache,
            rv,
            r_init,
            .init_I(roots),
            mino_time,
            r_sign,
        ),
        .case_II => radialAntiderivatives_case_I_and_II(
            T,
            geometry,
            cache,
            rv,
            r_init,
            .init_II(roots),
            mino_time,
            r_sign,
        ),
        .case_III => radialAntiderivatives_case_III(
            T,
            geometry,
            cache,
            rv,
            roots,
            mino_time,
            r_sign,
        ),
        .case_IV => radialAntiderivatives_case_IV(
            T,
            geometry,
            cache,
            rv,
            roots,
            mino_time,
            r_sign,
        ),
    };
}

fn Case_I_Case_II_Roots(comptime T: type) type {
    const A = T.Algebra;
    return struct {
        r1: T,
        r2: T,
        r3: T,
        r4: T,

        // The term multiplying the Mino time in (B30a) and (B48).
        // I: r1, II: r3
        r_A: T,
        // The term in the denominator of (B34) or (B54) in the first argument
        // to Pi that is not a difference of two roots.
        // I: r2, II: r4
        r_B: T,
        // Difference between two roots, as multiplyin Pi in (B30a) or (B48).
        // I: r21, II: r43
        r_diff1: T,
        // Difference between two roots, as in the numerator of the first
        // parameter to Pi in (B33) or (B53).
        // I: r32, II: r41
        r_diff2: T,

        fn init_I(roots: RadialRoots(T)) @This() {
            const r1 = roots.r1.real();
            const r2 = roots.r2.real();
            const r3 = roots.r3.real();
            const r4 = roots.r4.real();
            return .{
                .r1 = r1,
                .r2 = r2,
                .r3 = r3,
                .r4 = r4,
                .r_A = r1,
                .r_B = r2,
                .r_diff1 = A.sub(r2, r1),
                .r_diff2 = A.sub(r3, r2),
            };
        }

        fn init_II(roots: RadialRoots(T)) @This() {
            const r1 = roots.r1.real();
            const r2 = roots.r2.real();
            const r3 = roots.r3.real();
            const r4 = roots.r4.real();
            return .{
                .r1 = r1,
                .r2 = r2,
                .r3 = r3,
                .r4 = r4,
                .r_A = r3,
                .r_B = r4,
                .r_diff1 = A.sub(r4, r3),
                .r_diff2 = A.sub(r4, r1),
            };
        }
    };
}

/// The only things that differ between case I and case II are which radial
/// roots are used in the equations, so this common function is to avoid typing
/// out all of the equations twice.
fn radialAntiderivatives_case_I_and_II(
    comptime T: type,
    geometry: KerrMetric(T),
    cache: RadialCaseCache(T),
    rv: PrincipalRadialValues(T),
    r_init: T,
    roots: Case_I_Case_II_Roots(T),
    mino_time: T,
    r_sign: T.T,
) RadialAntiderivatives(T) {
    const A = T.Algebra;

    const r1 = roots.r1;
    const r2 = roots.r2;
    const r3 = roots.r3;
    const r4 = roots.r4;

    const r31 = A.sub(r3, r1);
    const r42 = A.sub(r4, r2);
    // r+/-
    const rp_A = A.sub(geometry.horizon_radius, roots.r_A);
    const rp_B = A.sub(geometry.horizon_radius, roots.r_B);
    const rm_A = A.sub(geometry.horizon_radius_negative, roots.r_A);
    const rm_B = A.sub(geometry.horizon_radius_negative, roots.r_B);

    // Evaluated at the initial location
    const asin_x = A.asin(cache.x);
    const am_X2 = elliptic_integrals.Jacobi.am(T, rv.X, cache.k);

    // Common terms
    const radial_sign: T = .promote(r_sign);
    const potential_init = potentials.radialFromRootsAlt(T, r1, r2, r3, r4, r_init);
    const r31rD2 = A.mult(r31, roots.r_diff2);
    const prefactor_E = A.sqrt(A.mult(r31, r42));
    const prefactor_Pi = A.div(.promote(2.0), prefactor_E);

    // Equation (B117). This equation is a unified inversion formula, and thus
    // will hold for case I, even though it was only derived for case II.
    const numerator_prefactor = A.mult(A.mult(r31rD2, roots.r_diff1), prefactor_E);
    const numerator = A.mult(
        A.mult(
            A.mult(numerator_prefactor, rv.integrals.sn),
            rv.integrals.cn,
        ),
        rv.integrals.dn,
    );

    // There is a typo in the paper. If you do the derivatives by hand, it's a
    // factor sn^2 on the denominator.
    const denominator = A.powi(
        A.sub(r31, A.mult(roots.r_diff2, A.powi(rv.integrals.sn, 2))),
        2,
    );
    const dr_dtau = A.div(numerator, denominator);

    // I: Equation (B31) and II: Equation (B51).
    const _H_term1 = A.div(dr_dtau, A.sub(rv.r, roots.r_A));
    const _H_term2 = A.div(A.sqrt(potential_init), A.sub(r_init, roots.r_A));
    const _H = A.sub(_H_term1, A.mult(radial_sign, _H_term2));

    // I: Equation (32) and II: Equation (B52), but they are identical.
    const _E_init = elliptic_integrals.Incomplete.secondKind(
        T,
        asin_x,
        cache.k,
    );
    const _E_final = elliptic_integrals.Incomplete.secondKind(
        T,
        am_X2,
        cache.k,
    );
    const _E = A.mult(
        prefactor_E,
        A.sub(_E_final, A.mult(radial_sign, _E_init)),
    );

    // I: Equation (B33) and  II: Equation (B42) in Equation (B53)
    const _Pi_init = elliptic_integrals.Incomplete.thirdKind(
        T,
        A.div(roots.r_diff2, r31),
        asin_x,
        cache.k,
    );
    const _Pi_final = elliptic_integrals.Incomplete.thirdKind(
        T,
        A.div(roots.r_diff2, r31),
        am_X2,
        cache.k,
    );
    const _Pi = A.mult(
        prefactor_Pi,
        A.sub(_Pi_final, A.mult(radial_sign, _Pi_init)),
    );

    // I: Equation (B34) and II: Equation (B54)
    const prefactor_Pi_p = A.mult(prefactor_Pi, A.div(roots.r_diff1, A.mult(rp_A, rp_B)));
    const prefactor_Pi_m = A.mult(prefactor_Pi, A.div(roots.r_diff1, A.mult(rm_A, rm_B)));

    const _Pi_p_init = elliptic_integrals.Incomplete.thirdKind(T, A.div(
        A.mult(rp_A, roots.r_diff2),
        A.mult(rp_B, r31),
    ), asin_x, cache.k);
    const _Pi_p_final = elliptic_integrals.Incomplete.thirdKind(T, A.div(
        A.mult(rp_A, roots.r_diff2),
        A.mult(rp_B, r31),
    ), am_X2, cache.k);

    const _Pi_m_init = elliptic_integrals.Incomplete.thirdKind(T, A.div(
        A.mult(rm_A, roots.r_diff2),
        A.mult(rm_B, r31),
    ), asin_x, cache.k);
    const _Pi_m_final = elliptic_integrals.Incomplete.thirdKind(T, A.div(
        A.mult(rm_A, roots.r_diff2),
        A.mult(rm_B, r31),
    ), am_X2, cache.k);

    const _Pi_p = A.mult(prefactor_Pi_p, A.sub(_Pi_p_final, A.mult(radial_sign, _Pi_p_init)));
    const _Pi_m = A.mult(prefactor_Pi_m, A.sub(_Pi_m_final, A.mult(radial_sign, _Pi_m_init)));

    const i_2_term = A.mult(.promote(0.5), A.add(A.mult(r1, r4), A.mult(r2, r3)));
    return .{
        .I_0 = mino_time,
        // I: Equation (B30a) and II: Equation (B48).
        .I_1 = A.add(A.mult(roots.r_A, mino_time), A.mult(roots.r_diff1, _Pi)),
        // Equation(B30b) and Equation (B49), they are identical.
        .I_2 = A.sub(A.sub(_H, A.mult(i_2_term, mino_time)), _E),
        // I: Equation (B30c) and II: Equation (B50).
        .I_plus = A.add(A.div(mino_time, rp_A), _Pi_p).neg(),
        .I_minus = A.add(A.div(mino_time, rm_A), _Pi_m).neg(),
        .am_X = am_X2,
    };
}

pub fn RValues(comptime T: type) type {
    return struct {
        const Self = @This();
        R_1: T,
        R_2: T,
    };
}

fn R_1_2_Auxillary(comptime T: type) type {
    return struct {
        const Self = @This();
        r: T,
        A: T,
        B: T,
        r1: T,
        r2: T,

        pub fn init(comptime D: type, r: D, A: D, B: D, r1: D, r2: D) Self {
            return .{
                .r = castDual(T.T, r),
                .A = castDual(T.T, A),
                .B = castDual(T.T, B),
                .r1 = castDual(T.T, r1),
                .r2 = castDual(T.T, r2),
            };
        }
    };
}

/// The same as `R_1_2` but with an auxillary structure that can be used to
/// refine the divergent part of the calculation to higher precision if needed.
fn R_1_2_aux(
    comptime T: type,
    alpha: T,
    phi: T,
    j: T,
    aux: R_1_2_Auxillary(ad.DualNumber(f128, T.N)),
) RValues(T) {
    return R_1_2_impl(T, alpha, phi, j, aux);
}

fn R_1_2(comptime T: type, alpha: T, phi: T, j: T) RValues(T) {
    return R_1_2_impl(T, alpha, phi, j, null);
}

fn R_1_2_impl(
    comptime T: type,
    alpha: T,
    phi: T,
    j: T,
    aux_term: ?R_1_2_Auxillary(ad.DualNumber(f128, T.N)),
) RValues(T) {
    const A = T.Algebra;
    const alpha_squared = A.powi(alpha, 2);

    const n = A.div(alpha_squared, A.sub(alpha_squared, .one));
    const _Pi = elliptic_integrals.Incomplete.thirdKind(T, n, phi, j);

    const sin_phi = A.sin(phi);
    const cos_phi = A.cos(phi);

    // Equation (B65)
    const p_1 = A.sqrt(A.div(
        A.sub(alpha_squared, .one),
        A.add(j, A.mult(A.sub(.one, j), alpha_squared)),
    ));

    std.debug.assert(p_1.x > 0);

    const j_sqrt_sin_phi = A.sqrt(A.sub(.one, A.mult(j, A.powi(sin_phi, 2))));
    const f_1_common = A.mult(
        p_1,
        j_sqrt_sin_phi,
    );
    const f_1 = A.mult(
        A.mult(p_1, .promote(0.5)),
        A.log(A.abs(A.div(
            A.add(f_1_common, sin_phi),
            A.sub(f_1_common, sin_phi),
        ))),
    );

    // Equation (B62)
    const bracket_1 = A.sub(_Pi, A.mult(alpha, f_1));
    const _R_1 = A.mult(A.div(.one, A.sub(.one, alpha_squared)), bracket_1);

    // Equation (B64)
    const _F = elliptic_integrals.Incomplete.firstKind(T, phi, j);

    const _E = elliptic_integrals.Incomplete.secondKind(T, phi, j);

    const prefactor_E = A.div(
        .one,
        A.add(j, A.mult(alpha_squared, A.sub(.one, j))),
    );

    // The numerical instability in this denominator is avoided by moving the
    // denominator away to the prefactor on the bracket.

    // The `alpha` here is taken in the denominator term.
    const _R_2_term_1_numerator = A.mult(sin_phi, j_sqrt_sin_phi);

    // This is the unstable bit, as it may be very small if:
    //
    //     `alpha * cos_phi ≈ -1`
    //
    // Although this can be still correctly computed, it is (in absolute terms)
    // inaccurate. This can't be avoided, due to double precision losses of `1 +
    // alpha * cos_phi`.
    const _divergent_denom = A.add(.one, A.mult(alpha, cos_phi));
    const _R_2_term_1_inv = b: {
        if (aux_term != null and _divergent_denom.x < 1e-4) {
            @branchHint(.unlikely);
            // Recalculate with more bits to keep the precision high:
            const aux = aux_term.?;

            const T128 = ad.DualNumber(f128, T.N);
            const A128 = T128.Algebra;

            // Need to recalculate x at higher precision because otherwise
            // there is not enough information int he bit representation.
            const _x128 = x_case_III(
                T128,
                aux.r,
                aux.A,
                aux.B,
                aux.r1,
                aux.r2,
            );

            const alpha128 = A128.div(
                A128.add(aux.B, aux.A),
                A128.sub(aux.B, aux.A),
            );

            const _inv_term = A128.div(
                alpha128,
                A128.add(.one, A128.mult(alpha128, _x128)),
            );

            break :b castDual(T.T, _inv_term);
        } else {
            break :b A.div(alpha, _divergent_denom);
        }
    };

    const _R_2_term_1 = A.mult(
        _R_2_term_1_numerator,
        _R_2_term_1_inv,
    );

    const bracket_2 = A.sub(_F, A.mult(
        A.mult(alpha_squared, prefactor_E),
        A.sub(_E, _R_2_term_1),
    ));
    const _R_2_term_2 = A.mult(_R_1, A.mult(
        prefactor_E,
        A.sub(
            A.mult(.promote(2.0), j),
            A.div(alpha_squared, A.sub(alpha_squared, .one)),
        ),
    ));

    const _R_2 = A.add(
        A.mult(A.div(.one, A.sub(alpha_squared, .one)), bracket_2),
        _R_2_term_2,
    );

    return .{
        .R_1 = _R_1,
        .R_2 = _R_2,
    };
}

test R_1_2 {
    const Dual = root.DualNumber(f64, 0);
    {
        const res = R_1_2(Dual, .promote(1.1), .promote(0.3), .promote(0.7));
        try std.testing.expectApproxEqAbs(0.14552777, res.R_1.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.0698569, res.R_2.x, TEST_TOLERANCE);
    }
    {
        const res = R_1_2(Dual, .promote(2.00131606), .promote(2.09375611), .promote(0.99977587));
        try std.testing.expectApproxEqAbs(16.97686880, res.R_1.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(2578.771114, res.R_2.x, TEST_TOLERANCE);
    }
    {
        const res = R_1_2(Dual, .promote(2.00131606), .promote(2.01014216), .promote(0.99977587));
        try std.testing.expectApproxEqAbs(10.153757494168755, res.R_1.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(15.42801619300267, res.R_2.x, TEST_TOLERANCE);
    }
    {
        const res = R_1_2(Dual, .promote(2.00131606436299), .promote(-2.010142159165531), .promote(0.9997758655364375));
        try std.testing.expectApproxEqAbs(-10.153737528360649, res.R_1.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-15.42799607, res.R_2.x, TEST_TOLERANCE);
    }
}

fn radialAntiderivatives_case_III(
    comptime T: type,
    geometry: KerrMetric(T),
    cache: RadialCaseCache(T),
    rv: PrincipalRadialValues(T),
    roots: RadialRoots(T),
    mino_time: T,
    r_sign: T.T,
) RadialAntiderivatives(T) {
    const A = T.Algebra;

    const r1 = roots.r1.real();
    const r2 = roots.r2.real();

    const _A = complexSubMultSqrt(T, roots.r3, roots.r4, r2);
    const _B = complexSubMultSqrt(T, roots.r3, roots.r4, r1);

    const rp1 = A.sub(geometry.horizon_radius, r1);
    const rp2 = A.sub(geometry.horizon_radius, r2);
    const rm1 = A.sub(geometry.horizon_radius_negative, r1);
    const rm2 = A.sub(geometry.horizon_radius_negative, r2);

    const r21 = A.sub(r2, r1);

    const radial_sign: T = .promote(r_sign);

    // Evaluated at the initial location
    const acos_x = A.acos(cache.x);
    const am_X3 = elliptic_integrals.Jacobi.am(T, rv.X, cache.k);

    // Equation (B66)
    const alpha_p = A.div(
        .promote(-1),
        x_case_III(T, geometry.horizon_radius, _A, _B, r1, r2),
    );
    const alpha_m = A.div(
        .promote(-1),
        x_case_III(T, geometry.horizon_radius_negative, _A, _B, r1, r2),
    );

    const _B_minus_A = A.sub(_B, _A);
    const _B_plus_A = A.add(_B, _A);

    // Equation (B58)
    const alpha = A.div(_B_plus_A, _B_minus_A);

    const _R_final = R_1_2(T, alpha, am_X3, cache.k);
    const _R_init = R_1_2_aux(
        T,
        alpha,
        acos_x,
        cache.k,
        .init(T, cache.r, _A, _B, r1, r2),
    );

    // Equation (B81)
    const prefactor_common = A.mult(
        A.mult(.promote(2), r21),
        A.sqrt(A.mult(_A, _B)),
    );
    const prefactor_Pi = A.div(
        prefactor_common,
        A.sub(A.powi(_B, 2), A.powi(_A, 2)),
    );
    const _Pi_1 = A.mult(
        prefactor_Pi,
        A.sub(_R_final.R_1, A.mult(radial_sign, _R_init.R_1)),
    );
    const _Pi_2 = A.mult(
        A.powi(prefactor_Pi, 2),
        A.sub(_R_final.R_2, A.mult(radial_sign, _R_init.R_2)),
    );

    // Equation (B82)
    const prefactor_Pi_p = A.div(
        prefactor_common,
        A.sub(A.mult(_B, rp2), A.mult(_A, rp1)),
    );
    const prefactor_Pi_m = A.div(
        prefactor_common,
        A.sub(A.mult(_B, rm2), A.mult(_A, rm1)),
    );

    const _R_p_final = R_1_2(T, alpha_p, am_X3, cache.k);
    const _R_p_init = R_1_2(T, alpha_p, acos_x, cache.k);
    const _R_m_final = R_1_2(T, alpha_m, am_X3, cache.k);
    const _R_m_init = R_1_2(T, alpha_m, acos_x, cache.k);

    const _Pi_p = A.mult(
        prefactor_Pi_p,
        A.sub(_R_p_final.R_1, A.mult(radial_sign, _R_p_init.R_1)),
    );
    const _Pi_m = A.mult(
        prefactor_Pi_m,
        A.sub(_R_m_final.R_1, A.mult(radial_sign, _R_m_init.R_1)),
    );

    // Equation (B78)
    const common = A.div(
        A.add(A.mult(_B, r2), A.mult(_A, r1)),
        A.add(_B, _A),
    );
    const _I_1 = A.add(A.mult(common, mino_time), _Pi_1);

    // Equation (B79)
    const _I_2 = A.add(
        A.add(
            A.mult(A.powi(common, 2), mino_time),
            A.mult(A.mult(.promote(2), common), _Pi_1),
        ),
        A.mult(A.sqrt(A.mult(_A, _B)), _Pi_2),
    );

    // Equation (B80)
    const denominator_I_p = A.add(A.mult(_B, rp2), A.mult(_A, rp1));
    const denominator_I_m = A.add(A.mult(_B, rm2), A.mult(_A, rm1));

    const _I_p = A.div(
        A.add(A.mult(A.add(_B, _A), mino_time), _Pi_p),
        denominator_I_p,
    );
    const _I_m = A.div(
        A.add(A.mult(A.add(_B, _A), mino_time), _Pi_m),
        denominator_I_m,
    );

    return .{
        .I_0 = mino_time,
        .I_1 = _I_1,
        .I_2 = _I_2,
        .I_plus = _I_p.neg(),
        .I_minus = _I_m.neg(),
        .am_X = am_X3,
    };
}

test "antiderivatives case III" {
    const Dual = ad.DualNumber(f64, 0);
    const metric = KerrMetric(Dual).init(.one, .promote(0.01));
    const mino_time: Dual = .promote(0.03487998131271285);
    const roots: RadialRoots(Dual) = .{
        .r1 = .fromReIm(.promote(-5.989857015584385), .zero),
        .r2 = .fromReIm(.promote(0.00004654463766229355), .zero),
        .r3 = .fromReIm(.promote(2.9949052354733614), .promote(-0.13464367601432883)),
        .r4 = .fromReIm(.promote(2.9949052354733614), .promote(0.13464367601432883)),
    };
    const cache_obs = radialCache_case_III(Dual, roots, .promote(1e4));
    const rv = fromCache_case_III(Dual, cache_obs, mino_time, -1);

    const antiderivs = radialAntiderivatives_case_III(
        Dual,
        metric,
        cache_obs,
        rv,
        roots,
        mino_time,
        -1,
    );

    try std.testing.expectApproxEqAbs(0.03487998, antiderivs.I_0.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(5.86001185, antiderivs.I_1.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(9971.71683053, antiderivs.I_2.x, TEST_TOLERANCE);
}

pub fn SValues(comptime T: type) type {
    return struct {
        const Self = @This();
        S_1: T,
        S_2: T,
    };
}

fn S_1_2(
    comptime T: type,
    alpha: T,
    phi: T,
    j: T,
    comptime what: enum { s1, s1_s2 },
) SValues(T) {
    const A = T.Algebra;
    const alpha_squared = A.powi(alpha, 2);

    const one_alpha_sq = A.add(.one, alpha_squared);
    const alpha_sq_j = A.sub(one_alpha_sq, j);
    const inv_one_alpha_sq = A.div(.one, one_alpha_sq);

    const _F = elliptic_integrals.Incomplete.firstKind(T, phi, j);
    const _Pi = elliptic_integrals.Incomplete.thirdKind(T, one_alpha_sq, phi, j);

    const sin_phi = A.sin(phi);

    // Equation (B95) for p
    const p_2 = A.sqrt(A.div(one_alpha_sq, alpha_sq_j));

    std.debug.assert(p_2.x > 0);

    // Equation (B95) for f
    const j_sqrt_sin_phi = A.sqrt(A.sub(.one, A.mult(j, A.powi(sin_phi, 2))));
    const f_2_common = A.mult(
        p_2,
        j_sqrt_sin_phi,
    );
    const f_2 = A.mult(
        A.mult(p_2, .promote(0.5)),
        A.log(A.abs(A.mult(
            A.div(A.sub(.one, p_2), A.add(.one, p_2)),
            A.div(
                A.add(.one, f_2_common),
                A.sub(.one, f_2_common),
            ),
        ))),
    );

    // Equation (B92)
    const bracket_1 = A.sub(
        A.add(_F, A.mult(alpha_squared, _Pi)),
        A.mult(alpha, f_2),
    );
    const _S_1 = A.mult(inv_one_alpha_sq, bracket_1);

    if (what == .s1) {
        return .{ .S_1 = _S_1, .S_2 = .zero };
    }

    // Required only for S2
    const _E = elliptic_integrals.Incomplete.secondKind(T, phi, j);
    const tan_phi = A.tan(phi);

    // Equation (B94)
    const prefactor_common = A.div(
        .one,
        A.mult(one_alpha_sq, alpha_sq_j),
    );
    const prefactor_S = A.add(
        inv_one_alpha_sq,
        A.div(A.sub(.one, j), alpha_sq_j),
    );

    const bracket_term_numerator = A.mult(
        A.mult(alpha_squared, j_sqrt_sin_phi),
        A.sub(alpha, tan_phi),
    );
    const bracket_term_denominator = A.add(.one, A.mult(alpha, tan_phi));
    const bracket_term = A.sub(
        A.div(bracket_term_numerator, bracket_term_denominator),
        A.mult(alpha_squared, alpha),
    );

    const term_1 = A.mult(
        prefactor_common,
        A.add(
            A.add(
                A.mult(A.sub(.one, j), _F),
                A.mult(alpha_squared, _E),
            ),
            bracket_term,
        ),
    );

    const _S_2 = A.add(term_1.neg(), A.mult(prefactor_S, _S_1));

    return .{
        .S_1 = _S_1,
        .S_2 = _S_2,
    };
}

test S_1_2 {
    const Dual = root.DualNumber(f64, 0);
    {
        const res = S_1_2(Dual, .promote(1.1), .promote(0.3), .promote(0.7), .s1_s2);
        try std.testing.expectApproxEqAbs(0.26132484, res.S_1.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.22685318, res.S_2.x, TEST_TOLERANCE);
    }
    {
        const res = S_1_2(
            Dual,
            .promote(2.00131606),
            .promote(2.09375611),
            .promote(0.99977587),
            .s1_s2,
        );
        try std.testing.expectApproxEqAbs(0.538731, res.S_1.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.43155129, res.S_2.x, TEST_TOLERANCE);
    }
    {
        const res = S_1_2(
            Dual,
            .promote(2.00131606),
            .promote(2.01014216),
            .promote(0.99977587),
            .s1_s2,
        );
        try std.testing.expectApproxEqAbs(0.60255912, res.S_1.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.40890383, res.S_2.x, TEST_TOLERANCE);
    }
    {
        const res = S_1_2(
            Dual,
            .promote(2.00131606436299),
            .promote(-2.010142159165531),
            .promote(0.9997758655364375),
            .s1_s2,
        );
        try std.testing.expectApproxEqAbs(0.22603822421598194, res.S_1.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.66611064, res.S_2.x, TEST_TOLERANCE);
    }
}

/// Calculate G&L Equation (B96).
fn g0_plus_minus(comptime T: type, g0: T, x4: T) T {
    const A = T.Algebra;
    return A.div(
        A.sub(A.mult(x4, g0), .one),
        A.add(g0, x4),
    );
}

fn radialAntiderivatives_case_IV(
    comptime T: type,
    geometry: KerrMetric(T),
    cache: RadialCaseCache(T),
    rv: PrincipalRadialValues(T),
    roots: RadialRoots(T),
    mino_time: T,
    r_sign: T.T,
) RadialAntiderivatives(T) {
    const A = T.Algebra;

    // Equation (B10)
    const b1 = roots.r3.real();
    const a2 = roots.r2.imag();

    const sp = cache.case_specific.case_IV;

    const phi_x = sp.phi;
    const am_X4 = elliptic_integrals.Jacobi.am(T, rv.X, cache.k);

    // Compute all the S terms
    const _S_init = S_1_2(T, sp.g0, phi_x, cache.k, .s1_s2);
    const _S_final = S_1_2(T, sp.g0, am_X4, cache.k, .s1_s2);

    const x4_p = A.div(A.add(geometry.horizon_radius, b1), a2);
    const gp = g0_plus_minus(T, sp.g0, x4_p);
    const _S1_p_init = S_1_2(T, gp, phi_x, cache.k, .s1).S_1;
    const _S1_p_final = S_1_2(T, gp, am_X4, cache.k, .s1).S_1;

    const x4_m = A.div(A.add(geometry.horizon_radius_negative, b1), a2);
    const gm = g0_plus_minus(T, sp.g0, x4_m);
    const _S1_m_init = S_1_2(T, gm, phi_x, cache.k, .s1).S_1;
    const _S1_m_final = S_1_2(T, gm, am_X4, cache.k, .s1).S_1;

    // The commonly reused coefficient with the r sign information
    const coeff = A.mult(.promote(r_sign), sp.coeff);

    const one_g0_sq = A.add(.one, A.powi(sp.g0, 2));

    // Equation (B115)
    const prefactor_Pi = A.mult(A.div(a2, sp.g0), one_g0_sq);
    const _Pi_1 = A.mult(
        A.mult(coeff, prefactor_Pi),
        A.sub(_S_final.S_1, _S_init.S_1),
    );
    const _Pi_2 = A.mult(
        A.mult(coeff, A.powi(prefactor_Pi, 2)),
        A.sub(_S_final.S_2, _S_init.S_2),
    );

    // Equation (B116), +ive horizon root:
    const prefactor_Pi_p = A.mult(
        coeff,
        A.div(
            one_g0_sq,
            A.mult(sp.g0, A.add(sp.g0, x4_p)),
        ),
    );
    const _Pi_p = A.mult(
        prefactor_Pi_p,
        A.sub(_S1_p_final, _S1_p_init),
    );

    // And for the, -ive horizon root:
    const prefactor_Pi_m = A.mult(
        coeff,
        A.div(
            one_g0_sq,
            A.mult(sp.g0, A.add(sp.g0, x4_m)),
        ),
    );
    const _Pi_m = A.mult(
        prefactor_Pi_m,
        A.sub(_S1_m_final, _S1_m_init),
    );

    const common_term = A.sub(A.div(a2, sp.g0), b1);

    // Equation (B112)
    const _I_1 = A.sub(A.mult(common_term, mino_time), _Pi_1);

    // Equation (B113)
    const _I_2 = A.add(
        A.sub(
            A.mult(A.powi(common_term, 2), mino_time),
            A.mult(A.mult(.promote(2), common_term), _Pi_1),
        ),
        _Pi_2,
    );

    // Equation (B114), +ive root
    const common_term_p = A.div(
        sp.g0,
        A.mult(
            a2,
            A.sub(.one, A.mult(sp.g0, x4_p)),
        ),
    );
    const _I_p = A.mult(common_term_p, A.sub(mino_time, _Pi_p));
    // and -ive root
    const common_term_m = A.div(
        sp.g0,
        A.mult(
            a2,
            A.sub(.one, A.mult(sp.g0, x4_m)),
        ),
    );
    const _I_m = A.mult(common_term_m, A.sub(mino_time, _Pi_m));

    return .{
        .I_0 = mino_time,
        .I_1 = _I_1,
        .I_2 = _I_2,
        .I_plus = _I_p,
        .I_minus = _I_m,
        .am_X = am_X4,
    };
}

test "antiderivatives case IV" {
    const Dual = ad.DualNumber(f64, 0);
    const metric = KerrMetric(Dual).init(.one, .promote(0.998));

    {
        const mino_time: Dual = .promote(0.5);
        const roots: RadialRoots(Dual) = .{
            .r1 = .fromReIm(.promote(-0.41365109572302766), .promote(-0.0028344710164291556)),
            .r2 = .fromReIm(.promote(-0.41365109572302766), .promote(0.0028344710164291556)),
            .r3 = .fromReIm(.promote(0.41365109572302755), .promote(-1.2850306998855008)),
            .r4 = .fromReIm(.promote(0.4136510957230279), .promote(1.2850306998855006)),
        };
        const cache_obs = radialCache_case_IV(Dual, roots, .promote(10.0));
        const rv = fromCache_case_IV(Dual, cache_obs, mino_time, -1);

        const antiderivs = radialAntiderivatives_case_IV(
            Dual,
            metric,
            cache_obs,
            rv,
            roots,
            mino_time,
            -1,
        );

        try std.testing.expectApproxEqAbs(1.57030413, antiderivs.am_X.x, TEST_TOLERANCE);

        try std.testing.expectApproxEqAbs(0.5, antiderivs.I_0.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(1.74403899, antiderivs.I_1.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(8.08724554, antiderivs.I_2.x, TEST_TOLERANCE);

        try std.testing.expectApproxEqAbs(0.38528573, antiderivs.I_plus.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.33585335, antiderivs.I_minus.x, TEST_TOLERANCE);
    }

    {
        const mino_time: Dual = .promote(0.1);
        const roots: RadialRoots(Dual) = .{
            .r1 = .fromReIm(.promote(-0.08980018566097332), .promote(-0.6613647911417029)),
            .r2 = .fromReIm(.promote(-0.08980018566097318), .promote(0.661364791141703)),
            .r3 = .fromReIm(.promote(0.08980018566097285), .promote(-0.8920322020410887)),
            .r4 = .fromReIm(.promote(0.08980018566097325), .promote(0.8920322020410887)),
        };
        const cache_obs = radialCache_case_IV(Dual, roots, .promote(1000000.0));
        const rv = fromCache_case_IV(Dual, cache_obs, mino_time, -1);

        const antiderivs = radialAntiderivatives_case_IV(
            Dual,
            metric,
            cache_obs,
            rv,
            roots,
            mino_time,
            -1,
        );

        try std.testing.expectApproxEqAbs(1.7784373616923916, antiderivs.am_X.x, TEST_TOLERANCE);

        try std.testing.expectApproxEqAbs(0.1, antiderivs.I_0.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(11.511918427473164, antiderivs.I_1.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(999989.9595452487, antiderivs.I_2.x, TEST_TOLERANCE);

        try std.testing.expectApproxEqAbs(0.005391509371629459, antiderivs.I_plus.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.005342055256626499, antiderivs.I_minus.x, TEST_TOLERANCE);
    }
}

/// The full set of anti-derivatives needed to calculate physical coordinates.
/// The notation follows G&L.
pub fn TotalAntiderivatives(comptime T: type) type {
    return struct {
        const Self = @This();

        radial: RadialAntiderivatives(T),
        angular: AngularAntiderivatives(T),
        theta_sign: T.T,
        angular_case: AngularCase,
        winding: usize,

        /// Calculate the coordinate time.
        ///
        /// Uses G&L Equation (12) along with Equation (B3).
        pub fn coordinateTime(self: Self, geometry: KerrMetric(T), geod: root.NullGeodesic(T)) T {
            const A = T.Algebra;

            // TODO: this was determined empirically. It seems to correct all
            // of the problems I've been having with tracking the time
            // coordinate, but I cannot necessarily motivate it mathematically.
            // In some sense it must be true, but it is perhaps offsetting
            // another error somewhere else. It is not needed for the azimuthal
            // coordinate.
            const winding: usize = 0;

            const rp = geometry.horizon_radius;
            const rm = geometry.horizon_radius_negative;
            const two_M = A.mult(.promote(2), geometry.M);
            const two_M_squared = A.powi(two_M, 2);

            const common = A.div(A.mult(geometry.a, geod.lambda), two_M);

            const term_plus = A.mult(A.mult(rp, A.sub(rp, common)), self.radial.I_plus);
            const term_minus = A.mult(A.mult(rm, A.sub(rm, common)), self.radial.I_minus);

            const bracket = A.sub(term_plus, term_minus);

            // Equation (B3)
            const I_t = A.add(
                A.mult(A.div(two_M_squared, A.sub(rp, rm)), bracket),
                A.add(
                    A.add(
                        A.mult(two_M_squared, self.radial.I_0),
                        A.mult(two_M, self.radial.I_1),
                    ),
                    self.radial.I_2,
                ),
            );

            // The turning points are tracked as in G&L Equation (51)
            const prefactor_init: T.T = if (@mod(winding, 2) == 0)
                1.0
            else
                -1.0;

            // G&L Equation (52)
            const _G_t = A.add(
                A.mult(
                    self.angular.G_t_half,
                    .promote(@floatFromInt(winding)),
                ),
                A.mult(
                    .promote(self.theta_sign),
                    A.sub(
                        A.mult(.promote(prefactor_init), self.angular.G_t_init),
                        self.angular.G_t_final,
                    ),
                ),
            );

            // TODO: this abs is used to avoid problems related to the signs of
            // theta, but it feels like an unsatisfying hat.
            return A.add(I_t, A.mult(A.powi(geometry.a, 2), A.abs(_G_t).neg()));
        }

        /// Compute the azimuthal coordinate from the total antiderivatives.
        ///
        /// Uses G&L Equation (11) along with Equation (B2).
        pub fn coordinateAzimuth(
            self: Self,
            geometry: KerrMetric(T),
            geod: root.NullGeodesic(T),
        ) T {
            const A = T.Algebra;

            // TODO: this should perhaps be a field for the result type
            // Two checks are needed because the theta_sign splits the upper
            // and lower half, whereas the windings split the left and right.
            // Therefore, to correct all of the quadrants, need to check both
            // possible combinations to find the false image:
            const is_false_image = false and ((self.theta_sign > 0 and self.winding > 0) or
                (self.theta_sign < 0 and self.winding > 1));

            const rp = geometry.horizon_radius;
            const rm = geometry.horizon_radius_negative;
            const two_M: T = A.mult(.promote(2), geometry.M);

            const a = if (is_false_image) geometry.a.neg() else geometry.a;

            const common = A.div(A.mult(a, geod.lambda), two_M);

            const term_plus = A.mult(A.sub(rp, common), self.radial.I_plus);
            const term_minus = A.mult(A.sub(rm, common), self.radial.I_minus);

            const bracket = A.sub(term_plus, term_minus);

            const I_phi = A.mult(
                A.div(A.mult(two_M, a), A.sub(rp, rm)),
                bracket,
            );

            // The turning points are tracked as in G&L Equation (51)
            const prefactor_init: T.T = if (@mod(self.winding, 2) == 0)
                1.0
            else
                -1.0;

            // Equation (52) or Equation (70) to get the prefactor
            var prefactor_half: T.T = @floatFromInt(self.winding);

            if (self.angular_case == .vortical) {
                // Equation (70)
                const additional = self.theta_sign * (1 - prefactor_init) / 2;
                prefactor_half += additional;
            }

            // G&L Equation (52)
            const _G_phi = A.add(
                A.mult(
                    .promote(prefactor_half),
                    self.angular.G_phi_half,
                ),
                A.mult(
                    .promote(self.theta_sign),
                    A.sub(
                        A.mult(.promote(prefactor_init), self.angular.G_phi_init),
                        self.angular.G_phi_final,
                    ),
                ),
            );

            const _phi = A.add(I_phi, A.mult(geod.lambda, _G_phi)).neg();

            // If we're calculating a false image, the sign of the overall
            // coordinate needs to change.
            if (is_false_image) {
                return _phi.neg();
            }
            return _phi;
        }
    };
}

/// Compute the total antiderivatives for a particular case of radial potential
/// roots. Requires `radialValues` first to have been calculated, as it reuses
/// many of their cached values.
///
/// TODO: pack the arguments into a struct.
pub fn totalAntiderivatives(
    comptime T: type,
    geometry: KerrMetric(T),
    case: RadialCase,
    angular_case: AngularCase,
    radial_cache: RadialCaseCache(T),
    radial_roots: RadialRoots(T),
    angular_cache: AngularCaseCache(T),
    angular_roots: AngularRoots(T),
    rv: PrincipalRadialValues(T),
    r_init: T,
    mino_time: T,
    theta_init: T,
    theta_final: T,
    r_sign: T.T,
    sign_theta: T.T,
    final_sign_theta: T.T,
    winding: usize,
) TotalAntiderivatives(T) {
    const radial = radialAntiderivatives(
        T,
        geometry,
        case,
        radial_cache,
        radial_roots,
        rv,
        r_init,
        mino_time,
        r_sign,
    );
    const angular = angularAntiderivatives(
        T,
        geometry,
        angular_case,
        angular_cache,
        angular_roots,
        theta_init,
        theta_final,
        sign_theta,
    );
    return .{
        .angular = angular,
        .radial = radial,
        .theta_sign = final_sign_theta,
        .winding = winding,
        .angular_case = angular_case,
    };
}
