/// References:
/// - A&S: Ambramowitz and Stegum
/// - F11: Fukushima (2011): Precise and fast computation of a general
///   incomplete elliptic integral of third kind by half and doubel argument
///   transformations.
const std = @import("std");
const ad = @import("zad");
const tracy = @import("tracy.zig");

const TEST_TOLERANCE = @import("options").test_numerical_tolerance;

fn maxOf(comptime T: type, values: []const T) T {
    var min = values[0];
    for (values[1..]) |v| {
        if (v.x > min.x) {
            min = v;
        }
    }
    return min;
}

fn minOf(comptime T: type, values: []const T) T {
    var min = values[0];
    for (values[1..]) |v| {
        if (v.x < min.x) {
            min = v;
        }
    }
    return min;
}

const Fukushima = struct {
    /// F11: Equation (25)
    fn T(comptime _T: type, t: _T, h: _T) _T {
        const A = _T.Algebra;
        if (h.x > 0) {
            const h_sqrt = A.sqrt(h);
            return A.div(A.atan(A.mult(t, h_sqrt)), h_sqrt);
        }
        if (h.x == 0) {
            @branchHint(.unlikely);
            return t;
        }
        const h_sqrt = A.sqrt(h.neg());
        const a = A.mult(t, h_sqrt);
        if (@abs(a.x) < 1) {
            return A.div(A.atanh(a), h_sqrt);
        } else {
            // From the Fukushima paper
            const n1 = A.abs(A.add(.one, a));
            const n2 = A.abs(A.sub(.one, a));
            return A.div(
                A.mult(.promote(0.5), A.log(A.div(n1, n2))),
                h_sqrt,
            );
        }
    }
};

test "Fukushima T" {
    const Dual = ad.DualNumber(f64, 1);
    {
        const res = Fukushima.T(Dual, .promote(-3.47), .promote(-0.083));
        try std.testing.expectApproxEqAbs(-15.265295566645671, res.x, TEST_TOLERANCE);
    }
}

/// Reference:
/// Bulrisch 1969, Numerical calculation of elliptic integrals and elliptic functions. III.
/// https://link.springer.com/article/10.1007/BF02165405
const Bulrisch = struct {
    // Adapted from the Algol procedure on p. 312.
    fn completeThirdKind(comptime T: type, _kc: T, _p: T, _a: T, _b: T) T {
        const error_tolerance = 1e-10;
        const A = T.Algebra;

        var kc = A.abs(_kc);
        var p = _p;
        var a = _a;
        var b = _b;

        var e = _kc;
        var m: T = .one;

        var f: T = .zero;
        var g: T = .zero;
        var q: T = .zero;

        if (p.x > 0) {
            p = A.sqrt(p);
            b = A.div(b, p);
        } else {
            f = A.powi(kc, 2);
            q = A.sub(.one, f);
            g = A.sub(.one, p);
            f = A.sub(f, p);
            q = A.mult(A.sub(b, A.mult(a, p)), q);
            p = A.sqrt(A.div(f, g));
            a = A.div(A.sub(a, b), g);
            b = A.sub(A.mult(a, p), A.div(q, A.mult(A.powi(g, 2), p)));
        }

        while (true) {
            f = a;
            const p_inv = A.div(.one, p);
            a = A.add(A.mult(p_inv, b), a);
            g = A.mult(e, p_inv);
            b = A.mult(.promote(2), A.add(A.mult(f, g), b));
            p = A.add(g, p);
            g = m;
            m = A.add(kc, m);
            if (@abs(g.x - kc.x) < (g.x * error_tolerance)) {
                break;
            }
            kc = A.mult(.promote(2), A.sqrt(e));
            e = A.mult(kc, m);
        }

        const numerator = A.mult(
            .promote(std.math.pi / 2.0),
            A.add(A.mult(a, m), b),
        );
        const denomiator = A.mult(m, A.add(m, p));
        return A.div(numerator, denomiator);
    }
};

test "Bulrisch.completeThirdKind" {
    const Dual = ad.DualNumber(f64, 1);
    {
        const result = Bulrisch.completeThirdKind(
            Dual,
            .promote(0.1),
            .derivative(0.2),
            .promote(0.5),
            .promote(0.1),
        );
        try std.testing.expectApproxEqAbs(1.8478186814949367, result.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-5.306722806676336, result.dx[0], TEST_TOLERANCE);
    }
    {
        const result = Bulrisch.completeThirdKind(
            Dual,
            .promote(0.1),
            .derivative(-0.2),
            .promote(0.5),
            .promote(0.1),
        );
        try std.testing.expectApproxEqAbs(-0.3947205295360626, result.x, TEST_TOLERANCE);
    }
}

/// The derivative formulae, assuming constant m.
const JacobiDerivatives = struct {
    pub fn sn(comptime T: type, cn_m: T, dn_m: T) T {
        const A = T.Algebra;
        return A.mult(cn_m, dn_m);
    }

    pub fn cn(comptime T: type, sn_m: T, dn_m: T) T {
        const A = T.Algebra;
        return A.mult(sn_m, dn_m).neg();
    }

    pub fn dn(comptime T: type, sn_m: T, cn_m: T, m: T) T {
        const A = T.Algebra;
        return A.mult(m, A.mult(sn_m, cn_m)).neg();
    }
};

pub const Jacobi = struct {
    const error_tolerance = 1e-10;

    pub fn Integrals(comptime T: type) type {
        return struct {
            sn: T,
            cn: T,
            dn: T,
            sc: T,
        };
    }

    /// The deriative with respect to `u` formulae, assuming `m` is constant.
    pub const Derivatives = JacobiDerivatives;

    /// Compute all Jacobi elliptic integrals in `Integrals`.
    ///
    /// Compare GSL-2.8: specfunc/elljac.c
    /// `specfunc/gsl_sf_elljac.h`
    ///
    /// Note that the parameter m is related to the argument as `m = 1 - kc2`.
    pub fn computeAll(comptime T: type, u: T, m: T) Integrals(T) {
        const ctx = tracy.trace(@src());
        defer ctx.end();
        const A = T.Algebra;

        // Reduction to non-elliptic functions
        if (@abs(m.x) < error_tolerance) {
            @branchHint(.unlikely);
            return .{
                .sn = A.sin(u),
                .cn = A.cos(u),
                .dn = .one,
                .sc = undefined,
            };
        }

        var emc = A.sub(.one, m);
        if (@abs(emc.x) < error_tolerance) {
            @branchHint(.unlikely);
            // A&S: 16.15.4: approximation in terms of hyperbolic functions
            const _cn = A.div(.one, A.cosh(u));
            return .{
                .sn = A.tanh(u),
                .cn = _cn,
                .dn = _cn,
                .sc = undefined,
            };
        }

        var d: T = .one;

        var _u = u;

        if (emc.x < 0) {
            d = A.sub(.one, emc);
            emc = A.mult(emc, d.neg());
            d = A.sqrt(d);
            _u = A.mult(_u, d);
        }

        var a: T = .one;
        var c: T = undefined;

        var em: [13]T = undefined;
        var en: [13]T = undefined;

        var index: usize = 0;
        for (0..em.len) |i| {
            index = i;
            em[i] = a;
            emc = A.sqrt(emc);
            en[i] = emc;
            c = A.mult(.promote(0.5), A.add(a, emc));
            if (@abs(a.x - emc.x) <= error_tolerance * a.x) {
                break;
            }
            emc = A.mult(emc, a);
            a = c;
        }

        _u = A.mult(_u, c);

        var _sn = A.sin(_u);
        var _cn = A.cos(_u);
        var _dn: T = .one;

        if (_sn.x != 0) {
            @branchHint(.likely);
            a = A.div(_cn, _sn);
            c = A.mult(c, a);
            for (0..index + 1) |i| {
                const j = index - i;
                a = A.mult(a, c);
                c = A.mult(c, _dn);
                _dn = A.div(A.add(en[j], a), A.add(em[j], a));
                a = A.div(c, em[j]);
            }
            a = A.div(.one, A.sqrt(A.add(A.powi(c, 2), .one)));
            _sn = A.mult(.promote(std.math.sign(_sn.x)), a);
            _cn = A.mult(c, _sn);
        }

        if (emc.x < 0) {
            a = _dn;
            _dn = _cn;
            _cn = a;
            _sn = A.div(_sn, d);
        }

        return .{
            .sn = _sn,
            .cn = _cn,
            .dn = _dn,
            .sc = undefined,
        };
    }

    fn amImpl(comptime T: type, u: T, m: T) T {
        const A = T.Algebra;
        const m_sub = A.sub(.one, m);
        if (m.x < error_tolerance) {
            // A&S 16.15.4
            const t = A.tanh(u);
            return A.add(
                A.asin(t),
                A.mult(
                    A.mult(
                        m_sub,
                        A.sub(t, A.mult(u, A.sub(.one, A.powi(t, 2)))),
                    ),
                    A.mult(A.cosh(u), .promote(1.0 / 4.0)),
                ),
            );
        }

        var a: T = .one;
        var b: T = A.sqrt(m_sub);
        var c: T = .zero;
        var tmp: T = .zero;

        var en: [10]T = undefined;

        var index: usize = 0;
        for (0..en.len) |i| {
            index = i;
            c = A.mult(A.sub(a, b), .promote(0.5));
            tmp = A.sqrt(A.mult(a, b));
            a = A.mult(A.add(a, b), .promote(0.5));
            b = tmp;
            en[index] = A.div(c, a);

            if (@abs(c.x) < error_tolerance) {
                break;
            }
        }

        a = A.mult(
            A.mult(a, u),
            .promote(@floatFromInt(std.math.powi(usize, 2, index + 1) catch {
                // this should be impossible, as index <= 10
                unreachable;
            })),
        );

        for (0..index + 1) |i| {
            const j = index - i;
            const phi = A.asin(A.mult(A.sin(a), en[j]));
            a = A.div(A.add(a, phi), .promote(2.0));
        }

        return a;
    }

    pub fn am(comptime T: type, u: T, m: T) T {
        const A = T.Algebra;
        if (m.x < 0) {
            const m_inv_sub = A.div(.one, A.sub(.one, m));
            const sqrt_m_sub = A.sqrt(m_inv_sub);

            const u_new = A.div(u, sqrt_m_sub);
            const m_new = A.mult(m, m_inv_sub).neg();

            const phi = amImpl(T, u_new, m_new);

            const sin_phi = A.sin(phi);
            const factor = std.math.divFloor(
                T.T,
                phi.x + std.math.pi / 2.0,
                std.math.pi,
            ) catch {
                unreachable;
            };

            const sign: T.T = if (@mod(factor, 2) == 0) 1.0 else -1.0;

            const a = A.asin(
                A.div(
                    A.mult(sqrt_m_sub, sin_phi),
                    A.sqrt(A.sub(.one, A.mult(m_new, A.powi(sin_phi, 2)))),
                ),
            );

            return A.add(
                .promote(factor * std.math.pi),
                A.mult(.promote(sign), a),
            );
        }
        return amImpl(T, u, m);
    }

    /// Compute all of the Jacobi elliptic integrals.
    pub fn all(comptime T: type, u: T, m: T) Integrals(T) {
        const A = T.Algebra;
        var a = computeAll(T, u, m);
        a.sc = A.div(a.sn, a.cn);
        return a;
    }

    pub fn sn(comptime T: type, u: T, m: T) T {
        return computeAll(T, u, m).sn;
    }

    pub fn cn(comptime T: type, u: T, m: T) T {
        return computeAll(T, u, m).cn;
    }

    pub fn dn(comptime T: type, u: T, m: T) T {
        return computeAll(T, u, m).dn;
    }

    pub fn sc(comptime T: type, u: T, m: T) T {
        return all(T, u, m).sc;
    }
};

test "jacobi integrals" {
    const Dual = ad.DualNumber(f64, 1);
    {
        const result = Jacobi.computeAll(Dual, .derivative(0.5), .promote(0.3));
        try std.testing.expectApproxEqAbs(0.47421562271182066, result.sn.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.8501921971455972, result.sn.dx[0], TEST_TOLERANCE);

        try std.testing.expectApproxEqAbs(0.8804087364264624, result.cn.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.45794005160670753, result.cn.dx[0], TEST_TOLERANCE);

        try std.testing.expectApproxEqAbs(0.9656789647459512, result.dn.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.12525107315562062, result.dn.dx[0], TEST_TOLERANCE);
    }
    {
        const result = Jacobi.computeAll(Dual, .derivative(-7.96), .promote(0.998));
        try std.testing.expectApproxEqAbs(-0.6327364690130791, result.cn.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.4907036872750576, result.cn.dx[0], TEST_TOLERANCE);
    }
    {
        const result = Jacobi.am(Dual, .promote(2.0), .promote(0.5));
        try std.testing.expectApproxEqAbs(1.6741639220482394, result.x, TEST_TOLERANCE);
    }
    {
        const result = Jacobi.am(Dual, .promote(2.0), .promote(-0.5));
        try std.testing.expectApproxEqAbs(2.268093077793418, result.x, TEST_TOLERANCE);
    }
    {
        const result = Jacobi.am(Dual, .promote(-9.6772169739112268), .promote(0.99977586553643749));
        try std.testing.expectApproxEqAbs(-2.010142159165532, result.x, TEST_TOLERANCE);
    }
}

pub const Carlson = struct {
    /// The error tolerance parameter for the Carlson algorithms.
    /// Tweaking this changes the number of loop iterations, and can make a
    /// quite significant change to the performance, since this function is the
    /// majority of the ray-tracing computation.
    ///
    /// GSL uses `error_tolerance` of 1e-3 for double precision, and 3e-2 for
    /// single precision. They quote the following relation:
    ///
    ///    relative error < 16 error_tolerance^6 / (1 - 2 error_tolerance)
    ///
    /// So an `error_tolerance` of 1e-3 has a precision of 1e-17, 0.01 has
    /// 2d-11, whereas 0.03 has 2e-8.
    const error_tolerance = 1e-2;

    /// Compute Carlson's incomplete elliptic integral of the first kind.
    ///
    /// Compare GSL: specfunc/ellint.c:L183
    /// `gsl_sf_ellint_RF_e`
    pub fn firstKind(comptime T: type, x: T, y: T, z: T) T {
        const ctx = tracy.trace(@src());
        defer ctx.end();
        const A = T.Algebra;
        const big = std.math.floatMax(T.T) / 3.0;
        const tiny = std.math.floatMin(T.T) * 3.0;

        // Once the numerical issues have been resolved that lead to nans being
        // passed, this could be removed
        std.debug.assert(!std.math.isNan(x.x));
        std.debug.assert(!std.math.isNan(y.x));
        std.debug.assert(!std.math.isNan(z.x));

        if ((@min(x.x, y.x, z.x) < 0.0) or
            (@min(x.x + y.x, x.x + z.x, y.x + z.x) < tiny) or
            (@max(x.x, y.x, z.x) > big))
        {
            @branchHint(.unlikely);
            // invalid arguments
            unreachable;
        }

        const C1 = 1.0 / 24.0;
        const C2 = 0.1;
        const C3 = 3.0 / 44.0;
        const C4 = 1.0 / 14.0;

        var xt = x;
        var yt = y;
        var zt = z;

        while (true) {
            const sqrt_x = A.sqrt(xt);
            const sqrt_y = A.sqrt(yt);
            const sqrt_z = A.sqrt(zt);

            const _lambda = A.add(
                A.mult(
                    sqrt_x,
                    A.add(sqrt_y, sqrt_z),
                ),
                A.mult(sqrt_y, sqrt_z),
            );

            xt = A.mult(A.add(xt, _lambda), .promote(0.25));
            yt = A.mult(A.add(yt, _lambda), .promote(0.25));
            zt = A.mult(A.add(zt, _lambda), .promote(0.25));

            const inv_average = A.div(
                .one,
                A.mult(A.add(A.add(xt, yt), zt), .promote(1.0 / 3.0)),
            );

            const delx = A.sub(.one, A.mult(xt, inv_average));
            const dely = A.sub(.one, A.mult(yt, inv_average));
            const delz = A.sub(.one, A.mult(zt, inv_average));

            // Typically takes 6-7 iterations at error_tolerance == 3e-4
            // Typically takes 4-5 iterations at error_tolerance == 1e-3
            if (@max(@abs(delx.x), @abs(dely.x), @abs(delz.x)) <= error_tolerance) {
                @branchHint(.unlikely);
                const e2 = A.sub(A.mult(delx, dely), A.powi(delz, 2));
                const e3 = A.mult(A.mult(delx, dely), delz);

                const inner_bracket = A.sub(
                    A.mult(.promote(C1), e2),
                    A.add(.promote(C2), A.mult(.promote(C3), e3)),
                );

                const outer_bracket = A.add(
                    A.add(
                        .one,
                        A.mult(inner_bracket, e2),
                    ),
                    A.mult(.promote(C4), e3),
                );

                return A.mult(outer_bracket, A.sqrt(inv_average));
            }
        }
    }

    /// Compute Carlson's incomplete elliptic integral of the second kind.
    ///
    /// Compare GSL: specfunc/ellint.c:L117
    /// `gsl_sf_ellint_RD_e`
    pub fn secondKind(comptime T: type, x: T, y: T, z: T) T {
        const ctx = tracy.trace(@src());
        defer ctx.end();
        const A = T.Algebra;
        const big = std.math.floatMax(T.T) / 3.0;
        const tiny = std.math.floatMin(T.T) * 3.0;

        if ((@min(x.x, y.x) < 0.0) or
            (@min(x.x + y.x, z.x) < tiny) or
            (@max(x.x, y.x, z.x) > big))
        {
            @branchHint(.unlikely);
            // invalid arguments
            unreachable;
        }

        const C1 = 3.0 / 14.0;
        const C2 = 1.0 / 6.0;
        const C3 = 9.0 / 22.0;
        const C4 = 3.0 / 26.0;
        const C5 = 0.25 * C3;
        const C6 = 1.5 * C4;

        var xt = x;
        var yt = y;
        var zt = z;

        var sum: T = .zero;
        var fac: T = .one;

        while (true) {
            const sqrt_x = A.sqrt(xt);
            const sqrt_y = A.sqrt(yt);
            const sqrt_z = A.sqrt(zt);

            const _lambda = A.add(
                A.mult(
                    sqrt_x,
                    A.add(sqrt_y, sqrt_z),
                ),
                A.mult(sqrt_y, sqrt_z),
            );

            sum = A.add(
                sum,
                A.div(
                    fac,
                    A.mult(sqrt_z, A.add(zt, _lambda)),
                ),
            );

            fac = A.mult(fac, .promote(0.25));

            xt = A.mult(A.add(xt, _lambda), .promote(0.25));
            yt = A.mult(A.add(yt, _lambda), .promote(0.25));
            zt = A.mult(A.add(zt, _lambda), .promote(0.25));

            const average = A.mult(
                .promote(0.2),
                A.add(
                    A.add(xt, yt),
                    A.mult(.promote(3.0), zt),
                ),
            );

            const delx = A.div(A.sub(average, xt), average);
            const dely = A.div(A.sub(average, yt), average);
            const delz = A.div(A.sub(average, zt), average);

            if (@max(@abs(delx.x), @abs(dely.x), @abs(delz.x)) <= error_tolerance) {
                const ea = A.mult(delx, dely);
                const eb = A.mult(delz, delz);
                const ec = A.sub(ea, eb);
                const ed = A.sub(ea, A.mult(.promote(6.0), eb));
                const ee = A.add(A.add(ed, ec), ec);

                const bracket_1 = A.mult(ed, A.sub(
                    A.sub(
                        A.mult(.promote(C5), ed),
                        A.mult(A.mult(.promote(C6), delz), ee),
                    ),
                    .promote(C1),
                ));

                const bracket_2 = A.add(
                    A.mult(.promote(-C3), ec),
                    A.mult(delz, A.mult(.promote(C4), ea)),
                );

                const term_1 = A.mult(
                    delz,
                    A.add(A.mult(.promote(C2), ee), A.mult(delz, bracket_2)),
                );

                const numerator = A.mult(
                    fac,
                    A.add(A.add(.one, bracket_1), term_1),
                );

                const term_2 = A.div(numerator, A.mult(average, A.sqrt(average)));

                return A.add(A.mult(.promote(3.0), sum), term_2);
            }
        }
    }

    /// Compute Carlson's degenerate elliptic integral.
    ///
    /// Compare GSL-2.8: specfunc/ellint.c:L74
    /// `gsl_sf_ellint_RC_e`
    pub fn degenerate(comptime T: type, x: T, y: T) T {
        const ctx = tracy.trace(@src());
        defer ctx.end();
        const A = T.Algebra;
        const sqrtny = 1.3e-19;
        const big = std.math.floatMax(T.T) / 3.0;
        const tiny = std.math.floatMin(T.T) * 3.0;
        const tinybig = big * tiny;
        const comp1 = 2.236 / sqrtny;
        const comp2 = tinybig * tinybig / 25.0;

        if ((x.x < 0.0) or
            (y.x == 0.0) or
            (x.x + @abs(y.x) < tiny) or
            (x.x + @abs(y.x) > big) or
            ((y.x < -comp1) and (x.x > 0) and (x.x < comp2)))
        {
            @branchHint(.unlikely);
            // invalid arguments
            unreachable;
        }

        const C1 = 0.3;
        const C2 = 1.0 / 7.0;
        const C3 = 0.375;
        const C4 = 9.0 / 22.0;

        var xt = x;
        var yt = y;
        var w: T = .one;

        if (y.x <= 0) {
            xt = A.sub(x, y);
            yt = yt.neg();
            w = A.div(A.sqrt(x), A.sqrt(xt));
        }

        while (true) {
            const alamb = A.add(
                A.mult(.promote(2.0), A.mult(A.sqrt(xt), A.sqrt(yt))),
                yt,
            );
            xt = A.mult(.promote(0.25), A.add(xt, alamb));
            yt = A.mult(.promote(0.25), A.add(yt, alamb));
            const average = A.mult(.promote(1.0 / 3.0), A.add(A.add(yt, yt), xt));
            const s = A.div(
                A.sub(yt, average),
                average,
            );

            if (@abs(s.x) < error_tolerance) {
                var term_1 = A.add(.promote(C3), A.mult(s, .promote(C4)));
                term_1 = A.add(.promote(C2), A.mult(s, term_1));
                term_1 = A.add(.promote(C1), A.mult(s, term_1));
                term_1 = A.add(.one, A.mult(A.powi(s, 2), term_1));
                return A.div(A.mult(w, term_1), A.sqrt(average));
            }
        }
    }

    /// Compute Carlson's incomplete elliptic integral of the third kind.
    ///
    /// Compare GSL-2.8: specfunc/ellint.c:L242
    /// `gsl_sf_ellint_RJ_e`
    pub fn thirdKind(comptime T: type, x: T, y: T, z: T, p: T) T {
        const ctx = tracy.trace(@src());
        defer ctx.end();
        const A = T.Algebra;
        const big = std.math.floatMax(T.T) / 3.0;
        const tiny = std.math.floatMin(T.T) * 3.0;

        if ((@min(x.x, y.x, z.x) < 0.0) or
            (@min(x.x + y.x, x.x + z.x, y.x + z.x) < tiny) or
            (@max(x.x, y.x, z.x, @abs(p.x)) > big))
        {
            @branchHint(.unlikely);
            // invalid arguments
            unreachable;
        }

        const C1 = 3.0 / 14.0;
        const C2 = 1.0 / 3.0;
        const C3 = 3.0 / 22.0;
        const C4 = 3.0 / 26.0;
        const C5 = 0.75 * C3;
        const C6 = 1.5 * C4;
        const C7 = 0.5 * C2;
        const C8 = 2.0 * C3;

        var xt = x;
        var yt = y;
        var zt = z;
        var pt = p;

        var a: T = .zero;
        var b: T = .zero;
        var rcx: T = .zero;

        var sum: T = .zero;
        var fac: T = .one;

        if (p.x <= 0) {
            xt = minOf(T, &.{ x, y, z });
            zt = maxOf(T, &.{ x, y, z });
            yt = A.sub(A.add(A.add(x, y), z), A.add(xt, zt));

            a = A.div(.one, A.sub(yt, p));
            b = A.mult(A.mult(a, A.sub(zt, yt)), A.sub(yt, xt));
            pt = A.add(yt, b);
            const rho = A.div(A.mult(xt, zt), yt);
            const tau = A.div(A.mult(p, pt), yt);
            rcx = degenerate(T, rho, tau);
        }

        while (true) {
            const sqrt_x = A.sqrt(xt);
            const sqrt_y = A.sqrt(yt);
            const sqrt_z = A.sqrt(zt);

            const _lambda = A.add(
                A.mult(
                    sqrt_x,
                    A.add(sqrt_y, sqrt_z),
                ),
                A.mult(sqrt_y, sqrt_z),
            );

            const alpha = A.powi(
                A.add(
                    A.mult(pt, A.add(A.add(sqrt_x, sqrt_y), sqrt_z)),
                    A.mult(A.mult(sqrt_x, sqrt_y), sqrt_z),
                ),
                2,
            );

            const beta = A.mult(pt, A.powi(A.add(pt, _lambda), 2));
            sum = A.add(sum, A.mult(fac, degenerate(T, alpha, beta)));

            fac = A.mult(fac, .promote(0.25));

            xt = A.mult(A.add(xt, _lambda), .promote(0.25));
            yt = A.mult(A.add(yt, _lambda), .promote(0.25));
            zt = A.mult(A.add(zt, _lambda), .promote(0.25));
            pt = A.mult(A.add(pt, _lambda), .promote(0.25));

            const average = A.mult(
                .promote(0.2),
                A.add(A.add(A.add(xt, yt), zt), A.add(pt, pt)),
            );

            const delx = A.div(A.sub(average, xt), average);
            const dely = A.div(A.sub(average, yt), average);
            const delz = A.div(A.sub(average, zt), average);
            const delp = A.div(A.sub(average, pt), average);

            if (@max(@abs(delx.x), @abs(dely.x), @abs(delz.x), @abs(delp.x)) <= error_tolerance) {
                const ea = A.add(
                    A.mult(delx, A.add(dely, delz)),
                    A.mult(dely, delz),
                );
                const eb = A.mult(delx, A.mult(dely, delz));
                const ec = A.mult(delp, delp);
                const ed = A.sub(ea, A.mult(.promote(3.0), ec));
                const ee = A.add(eb, A.mult(.promote(2.0), A.mult(delp, A.sub(ea, ec))));

                const term_1 = A.mult(
                    ed,
                    A.sub(
                        A.mult(.promote(C5), ed),
                        A.add(A.mult(.promote(C6), ee), .promote(C1)),
                    ),
                );
                const term_2 = A.mult(
                    eb,
                    A.add(
                        .promote(C7),
                        A.mult(delp, A.sub(A.mult(delp, .promote(C4)), .promote(C8))),
                    ),
                );
                const term_3 = A.mult(
                    delp,
                    A.mult(
                        ea,
                        A.sub(.promote(C2), A.mult(delp, .promote(C3))),
                    ),
                );
                const term_4 = A.mult(ec, A.mult(.promote(C2), delp));

                const bracket = A.add(.one, A.add(term_1, A.add(term_2, A.add(term_3, term_4))));

                const denominator = A.mult(average, A.sqrt(average));

                const ans = A.add(
                    A.mult(.promote(3.0), sum),
                    A.div(A.mult(fac, bracket), denominator),
                );

                if (p.x <= 0) {
                    const intermediate = A.mult(
                        .promote(3.0),
                        A.sub(rcx, firstKind(T, xt, yt, zt)),
                    );
                    return A.mult(a, A.add(A.mult(b, ans), intermediate));
                } else {
                    return ans;
                }
            }
        }
    }

    /// The common shortening of Carlson's incomplete elliptic integral of the
    /// first kind.
    pub const rf = firstKind;
    /// The common shortening of Carlson's incomplete elliptic integral of the
    /// second kind.
    pub const rd = secondKind;
    /// The common shortening of Carlson's incomplete elliptic integral of the
    /// third kind.
    pub const rj = thirdKind;
    /// The common shortening of Carlson's degenerate elliptic integral.
    pub const rc = degenerate;
};

test "Carlson.firstKind" {
    const Dual = ad.DualNumber(f64, 1);
    {
        // From GSL
        const result = Carlson.firstKind(Dual, .promote(5.0e-11), .promote(1.0e-10), .derivative(1.0));
        try std.testing.expectApproxEqAbs(12.36441982979439, result.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-5.682209915322291, result.dx[0], TEST_TOLERANCE);
    }
    {
        // From GSL
        const result = Carlson.firstKind(Dual, .promote(1.0), .promote(2.0), .derivative(3.0));
        try std.testing.expectApproxEqAbs(0.7269459354689082, result.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.04841004683817687, result.dx[0], TEST_TOLERANCE);
    }
    {
        const result = Carlson.firstKind(Dual, .promote(0.3), .promote(0.6), .derivative(0.1));
        try std.testing.expectApproxEqAbs(1.8529130412683374, result.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-1.9786071064514628, result.dx[0], TEST_TOLERANCE);
    }
}

test "Carlson.secondKind" {
    const Dual = ad.DualNumber(f64, 1);
    {
        // From GSL
        const result = Carlson.secondKind(Dual, .promote(5.0e-11), .promote(1.0e-10), .derivative(1.0));
        try std.testing.expectApproxEqAbs(34.0932594919337362, result.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-49.63988924161392, result.dx[0], TEST_TOLERANCE);
    }
    {
        // From GSL
        const result = Carlson.secondKind(Dual, .promote(1.0), .promote(2.0), .derivative(3.0));
        try std.testing.expectApproxEqAbs(0.2904602810289906, result.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.09460511517373632, result.dx[0], TEST_TOLERANCE);
    }
    {
        const result = Carlson.secondKind(Dual, .promote(0.3), .promote(0.6), .derivative(0.1));
        try std.testing.expectApproxEqAbs(11.871642638708785, result.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-90.35092388499453, result.dx[0], TEST_TOLERANCE);
    }
}

test "Carlson.degenerate" {
    const Dual = ad.DualNumber(f64, 1);
    {
        // From GSL
        const result = Carlson.degenerate(Dual, .promote(1.0), .derivative(2.0));
        try std.testing.expectApproxEqAbs(0.7853981633974482, result.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.1426990816987265, result.dx[0], TEST_TOLERANCE);
    }
    {
        const result = Carlson.degenerate(Dual, .promote(0.3), .derivative(0.6));
        try std.testing.expectApproxEqAbs(1.433934302386369, result.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.8684389553518207, result.dx[0], TEST_TOLERANCE);
    }
}

test "Carlson.thirdKind" {
    const Dual = ad.DualNumber(f64, 1);
    {
        // From GSL
        const result = Carlson.thirdKind(Dual, .promote(2.0), .derivative(3.0), .promote(4.0), .promote(5.0));
        try std.testing.expectApproxEqAbs(0.1429757966715675, result.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.013579327237123053, result.dx[0], TEST_TOLERANCE);
    }
    {
        const result = Carlson.thirdKind(Dual, .promote(0.3), .derivative(0.6), .promote(0.4), .promote(0.1));
        try std.testing.expectApproxEqAbs(7.8308384830546895, result.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-4.890649519526418, result.dx[0], TEST_TOLERANCE);
    }
}

pub const Complete = struct {
    /// Compute the Legendre complete elliptic integral of the first kind.
    pub fn firstKind(comptime T: type, parameter: T) T {
        const ctx = tracy.trace(@src());
        defer ctx.end();
        // TODO: some bounds checks on the parameter
        const A = T.Algebra;
        return Carlson.firstKind(T, .zero, A.sub(.one, parameter), .one);
    }

    /// Compute the Legendre complete elliptic integral of the second kind.
    pub fn secondKind(comptime T: type, parameter: T) T {
        const ctx = tracy.trace(@src());
        defer ctx.end();
        const A = T.Algebra;
        // TODO: this could be reused!
        const rd = Carlson.secondKind(T, .zero, A.sub(.one, parameter), .one);
        // TODO: this could be reused!
        const _K = Carlson.firstKind(T, .zero, A.sub(.one, parameter), .one);
        return A.sub(_K, A.mult(.promote(1.0 / 3.0), A.mult(parameter, rd)));
    }

    /// Compute the derivative of the complete second kind, with respect to the
    /// parameter k.
    pub fn secondKind_dk(comptime T: type, parameter: T) T {
        // TODO: check this expression. I guessed it from the incomplete version.
        const A = T.Algebra;
        // TODO: this could be reused!
        const rd = Carlson.secondKind(T, .zero, A.sub(.one, parameter), .one);
        const _K = Carlson.firstKind(T, .zero, A.sub(.one, parameter), .one);
        const _E = A.sub(_K, A.mult(.promote(1.0 / 3.0), A.mult(parameter, rd)));
        return A.sub(A.sub(_E, _K), parameter);
    }

    /// Compute the Legendre complete elliptic integral of the third kind.
    pub fn thirdKind(comptime T: type, n: T, m: T) T {
        const ctx = tracy.trace(@src());
        defer ctx.end();

        // TODO: some bounds checks on the parameter
        const A = T.Algebra;

        // TODO: benchmark these to check which is more performant

        // Alg 1: this one is invalid if n == 1.
        // const m_sub = A.sub(.one, m);
        // const _rf = Carlson.firstKind(T, .zero, m_sub, .one);
        // const _rj = Carlson.thirdKind(T, .zero, m_sub, .one, A.sub(.one, n));
        // return A.add(_rf, A.mult(_rj, A.mult(n, .promote(1.0 / 3.0))));

        // Alg 2:
        if (n.x == 0) {
            return firstKind(T, m);
        }
        if (n.x > 1) {
            return A.sub(
                firstKind(T, m),
                thirdKind(T, A.div(m, n), m),
            );
        }
        if (m.x == 0 or m.x == 1) {
            @branchHint(.unlikely);
            unreachable;
        }
        return Bulrisch.completeThirdKind(
            T,
            A.sqrt(A.sub(.one, m)),
            A.sub(.one, n),
            .one,
            .one,
        );
    }

    /// The common name for the complete elliptic integral of the first kind.
    pub const K = firstKind;
    /// The common name for the complete elliptic integral of the second kind.
    pub const E = secondKind;
    /// The common name for the complete elliptic integral of the third kind.
    pub const Pi = thirdKind;
};

test "complete first-kind" {
    const Dual = ad.DualNumber(f64, 1);
    {
        const result = Complete.firstKind(Dual, .derivative(0.7));
        try std.testing.expectApproxEqAbs(2.0753631352924686, result.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(1.4739562556144803, result.dx[0], TEST_TOLERANCE);
    }
}

test "complete second-kind" {
    const Dual = ad.DualNumber(f64, 1);
    {
        const result = Complete.secondKind(Dual, .derivative(0.7));
        try std.testing.expectApproxEqAbs(1.2416705679458224, result.x, TEST_TOLERANCE);
    }
}

test "complete third-kind" {
    const Dual = ad.DualNumber(f64, 1);
    {
        const result = Complete.thirdKind(Dual, .derivative(0.7), .promote(0.3));
        try std.testing.expectApproxEqAbs(3.21214277914048, result.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(5.574437463416906, result.dx[0], TEST_TOLERANCE);
    }
    {
        const result = Complete.thirdKind(Dual, .derivative(1.0), .promote(0.3));
        try std.testing.expectApproxEqAbs(-0.3509149295535878, result.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.4620744013682388, result.dx[0], TEST_TOLERANCE);
    }
}

pub const Incomplete = struct {
    /// Compute the Legendre incomplete elliptic integral of the first kind.
    pub fn firstKind(comptime T: type, amplitude: T, parameter: T) T {
        const ctx = tracy.trace(@src());
        defer ctx.end();
        const A = T.Algebra;

        if (amplitude.x == 0) return .zero;

        if (parameter.x > 1) {
            // A&S 17.4.15: parameter greater than unity
            const m_sqrt = A.sqrt(parameter);
            var theta = A.asin(A.mult(A.sin(amplitude), m_sqrt));
            const sign_theta = std.math.sign(theta.x);
            theta = A.mult(theta, .promote(sign_theta));
            return A.mult(.promote(sign_theta), A.div(
                firstKind(T, theta, A.div(.one, parameter)),
                m_sqrt,
            ));
        } else if (parameter.x < 0) {
            // A&S 17.4.17: negative parameter
            const m_minus = A.sub(.one, parameter);
            const m_invsqrt = A.div(
                .one,
                A.sqrt(m_minus),
            );
            const ratio = A.div(
                parameter.neg(),
                m_minus,
            );

            const theta = A.sub(.promote(std.math.pi / 2.0), amplitude);
            const complete = Complete.firstKind(T, ratio);
            const incomplete = firstKind(T, theta, ratio);

            return A.mult(m_invsqrt, A.sub(complete, incomplete));
        }

        const sign_amplitude = std.math.sign(amplitude.x);
        var phi = A.abs(amplitude);
        var k: T.T = 0;
        while (phi.x > std.math.pi / 2.0) {
            k += sign_amplitude;
            phi.x -= std.math.pi;
        }

        const sin_phi = A.sin(phi);
        const sin_phi_2 = A.powi(sin_phi, 2);

        const result = Carlson.firstKind(
            T,
            A.sub(.one, sin_phi_2),
            A.sub(.one, A.mult(parameter, sin_phi_2)),
            .one,
        );

        // carry the sign information
        var signed_result = A.mult(.promote(sign_amplitude), A.mult(sin_phi, result));

        if (k != 0) {
            signed_result = A.add(
                signed_result,
                A.mult(
                    .promote(2 * k),
                    Complete.firstKind(T, parameter),
                ),
            );
        }

        return signed_result;
    }

    /// Compute the Legendre incomplete elliptic integral of the second kind.
    pub fn secondKind(comptime T: type, amplitude: T, parameter: T) T {
        const ctx = tracy.trace(@src());
        defer ctx.end();
        const A = T.Algebra;

        if (amplitude.x == 0) return .zero;

        if (parameter.x > 1) {
            // A&S 17.4.16: parameter greater than unity
            unreachable;
            // const m_sqrt = A.sqrt(parameter);
            // const m_inv = A.div(.one, parameter);
            // const new_amplitude = A.mult(amplitude, m_sqrt);

            // return A.sub(
            //     A.mult(m_sqrt, secondKind(T, new_amplitude, m_inv)),
            //     A.mult(A.sub(parameter, .one), amplitude),
            // );
            // } else if (parameter.x < 0) {
            // A&S 17.4.18: negative parameter
            // unreachable;
        }

        const sign_amplitude = std.math.sign(amplitude.x);
        var phi = A.abs(amplitude);
        var k: T.T = 0;
        while (phi.x > std.math.pi / 2.0) {
            k += sign_amplitude;
            phi.x -= std.math.pi;
        }

        const sin_phi = A.sin(phi);
        const cos_phi_squared = A.powi(A.cos(phi), 2);
        const sin_phi_squared = A.sub(.one, cos_phi_squared);
        const a = A.sub(.one, A.mult(sin_phi_squared, parameter));

        const _rf = Carlson.firstKind(T, cos_phi_squared, a, .one);
        const _rd = Carlson.secondKind(T, cos_phi_squared, a, .one);

        const value = A.sub(_rf, A.mult(
            A.mult(A.mult(sin_phi_squared, parameter), _rd),
            .promote(1.0 / 3.0),
        ));

        // carry the sign information
        var signed_result = A.mult(.promote(sign_amplitude), A.mult(sin_phi, value));

        if (k != 0) {
            signed_result = A.add(
                signed_result,
                A.mult(
                    .promote(2 * k),
                    Complete.secondKind(T, parameter),
                ),
            );
        }

        return signed_result;
    }

    fn thirdKindImpl(comptime T: type, n: T, sin_phi: T, parameter: T) T {
        const A = T.Algebra;
        const sin_phi_squared = A.powi(sin_phi, 2);
        const cos_phi_squred = A.sub(.one, sin_phi_squared);

        const a = A.sub(.one, A.mult(parameter, sin_phi_squared));

        const _rf = Carlson.firstKind(T, cos_phi_squred, a, .one);
        const _rj = Carlson.thirdKind(T, cos_phi_squred, a, .one, A.sub(
            .one,
            A.mult(n, sin_phi_squared),
        ));

        return A.mult(
            sin_phi,
            A.add(
                _rf,
                A.mult(
                    .promote(1.0 / 3.0),
                    A.mult(_rj, A.mult(n, sin_phi_squared)),
                ),
            ),
        );
    }

    /// Compute the Legendre incomplete elliptic integral of the third kind.
    ///
    /// Compare `gsl_sf_ellint_P_e`. Note that the order of the arguments
    /// differs, and what is here `parameter` is actually `m = k * k`.
    pub fn thirdKind(comptime T: type, n: T, amplitude: T, parameter: T) T {
        const ctx = tracy.trace(@src());
        defer ctx.end();
        const A = T.Algebra;

        const sin_phi = A.sin(amplitude);

        if (parameter.x < 0) {
            const m_sub = A.sub(.one, parameter);
            const new_amplitude = A.asin(A.mult(
                A.sqrt(
                    A.div(
                        m_sub,
                        A.sub(.one, A.mult(parameter, A.powi(sin_phi, 2))),
                    ),
                ),
                sin_phi,
            ));

            const m_sub_inv = A.div(.one, m_sub);
            const new_n = A.mult(A.sub(n, parameter), m_sub_inv);
            const new_parameter = A.mult(parameter.neg(), m_sub_inv);

            const _Pi = thirdKind(T, new_n, new_amplitude, new_parameter);
            const _F = firstKind(T, new_amplitude, new_parameter);

            return A.mult(
                A.div(A.sqrt(m_sub_inv), new_n),
                A.add(A.mult(new_parameter, _F), A.mult(A.mult(m_sub_inv, n), _Pi)),
            );
        }

        if (n.x > 1) {
            const n_sub = A.sub(.one, n);
            const t = A.div(
                A.tan(amplitude),
                A.sqrt(A.sub(.one, A.mult(parameter, A.powi(sin_phi, 2)))),
            );
            const h = A.div(A.mult(n_sub, A.sub(n, parameter)), n);

            const _Pi = thirdKind(T, A.div(parameter, n), amplitude, parameter);
            const _F = firstKind(T, amplitude, parameter);

            return A.add(
                A.sub(
                    Fukushima.T(T, t, h),
                    _Pi,
                ),
                _F,
            );
        }

        if (@abs(amplitude.x) > (std.math.pi / 2.0)) {
            const sign_phi: T = .promote(std.math.sign(amplitude.x));
            var phi = amplitude;
            var factor: usize = 0;
            while (@abs(phi.x) > (std.math.pi / 2.0)) {
                factor += 1;
                phi = A.sub(phi, A.mult(sign_phi, .promote(std.math.pi)));
            }
            const _Pi_complete = A.mult(
                .promote(2.0 * sign_phi.x),
                A.mult(.promote(@floatFromInt(factor)), Complete.thirdKind(T, n, parameter)),
            );
            return A.sub(
                _Pi_complete,
                thirdKindImpl(
                    T,
                    n,
                    A.sin(phi.neg()),
                    parameter,
                ),
            );
        }

        return thirdKindImpl(T, n, sin_phi, parameter);
    }

    /// Derivative of the incomplete elliptic integral of the second kind with
    /// respect to k.
    pub fn secondKind_dk(comptime T: type, amplitude: T, parameter: T) T {
        const ctx = tracy.trace(@src());
        defer ctx.end();
        const A = T.Algebra;
        const _E = secondKind(T, amplitude, parameter);
        const _F = firstKind(T, amplitude, parameter);
        return A.div(A.sub(_E, _F), A.mult(.promote(2), parameter));
    }

    /// Common shortening for the elliptic integral of the first kind.
    pub const F = firstKind;
    /// Common shortening for the elliptic integral of the second kind.
    pub const E = secondKind;
    /// Common shortening for the elliptic integral of the third kind.
    pub const Pi = thirdKind;
};

test "incomplete first-kind" {
    const Dual = ad.DualNumber(f64, 1);
    {
        // From GSL
        const result = Incomplete.firstKind(Dual, .promote(std.math.pi / 3.0), .promote(0.99 * 0.99));
        try std.testing.expectApproxEqAbs(1.3065333392738766762, result.x, TEST_TOLERANCE);
    }
    {
        // From GSL
        const result = Incomplete.firstKind(Dual, .promote(std.math.pi / 3.0), .promote(0.5 * 0.5));
        try std.testing.expectApproxEqAbs(1.0895506700518854093, result.x, TEST_TOLERANCE);
    }
    {
        // From GSL
        const result = Incomplete.firstKind(Dual, .promote(std.math.pi / 3.0), .promote(0.01 * 0.01));
        try std.testing.expectApproxEqAbs(1.0472129063770918952, result.x, TEST_TOLERANCE);
    }
    {
        const result = Incomplete.firstKind(Dual, .derivative(0.5), .promote(0.7));
        try std.testing.expectApproxEqAbs(0.5150014689374526, result.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(1.0916706560156484, result.dx[0], TEST_TOLERANCE);
    }
    {
        const result = Incomplete.firstKind(Dual, .promote(0.5), .promote(2.7));
        try std.testing.expectApproxEqAbs(0.5789104025367218, result.x, TEST_TOLERANCE);
    }
    {
        const result = Incomplete.firstKind(Dual, .promote(3.1), .promote(0.7));
        try std.testing.expectApproxEqAbs(4.1091252207894655, result.x, TEST_TOLERANCE);
    }
    {
        const result = Incomplete.firstKind(Dual, .promote(5.1), .promote(0.7));
        try std.testing.expectApproxEqAbs(6.898726928887488, result.x, TEST_TOLERANCE);
    }
    {
        const result = Incomplete.firstKind(Dual, .promote(1.1), .promote(-0.7));
        try std.testing.expectApproxEqAbs(1.0030688571399822, result.x, TEST_TOLERANCE);
    }
    {
        const result = Incomplete.firstKind(Dual, .promote(2.09375611), .promote(0.99977587));
        try std.testing.expectApproxEqAbs(9.858271365970047, result.x, TEST_TOLERANCE);
    }
    {
        const result = Incomplete.firstKind(Dual, .promote(2.0937561097665736), .promote(0.9997758655364375));
        try std.testing.expectApproxEqAbs(9.858251462169509, result.x, TEST_TOLERANCE);
    }
    {
        // Check antisymmetry property
        const result = Incomplete.firstKind(Dual, .promote(-2.0937561097665736), .promote(0.9997758655364375));
        try std.testing.expectApproxEqAbs(-9.858251462169509, result.x, TEST_TOLERANCE);
    }
    {
        const result = Incomplete.firstKind(Dual, .promote(0.17454113151938627), .promote(0.09843479323515428));
        try std.testing.expectApproxEqAbs(0.17462795276790985, result.x, TEST_TOLERANCE);
    }
    {
        const result = Incomplete.firstKind(Dual, .derivative(2.16960), .promote(0.99929));
        try std.testing.expectApproxEqAbs(8.849064593413068, result.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(1.7727848102458736, result.dx[0], TEST_TOLERANCE);
    }
}

test "incomplete second-kind" {
    const Dual = ad.DualNumber(f64, 1);
    {
        // From GSL
        const result = Incomplete.secondKind(Dual, .promote(std.math.pi / 3.0), .promote(0.99 * 0.99));
        try std.testing.expectApproxEqAbs(0.8704819220377943536, result.x, TEST_TOLERANCE);
    }
    {
        // From GSL
        const result = Incomplete.secondKind(Dual, .promote(std.math.pi / 3.0), .promote(0.5 * 0.5));
        try std.testing.expectApproxEqAbs(1.0075555551444720293, result.x, TEST_TOLERANCE);
    }
    {
        // From GSL
        const result = Incomplete.secondKind(Dual, .promote(std.math.pi / 3.0), .promote(0.01 * 0.01));
        try std.testing.expectApproxEqAbs(1.0471821963889481104, result.x, TEST_TOLERANCE);
    }
    {
        const result = Incomplete.secondKind(Dual, .derivative(0.5), .promote(0.3));
        try std.testing.expectApproxEqAbs(0.4939911447289684, result.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.9649069104738658, result.dx[0], TEST_TOLERANCE);
    }
    {
        const result = Incomplete.secondKind(Dual, .promote(2.09375611), .promote(0.99977587));
        try std.testing.expectApproxEqAbs(1.1347451, result.x, TEST_TOLERANCE);
    }
    {
        const result = Incomplete.secondKind(Dual, .promote(2.010142159165531), .promote(0.9997758655364375));
        try std.testing.expectApproxEqAbs(1.0960437630027704, result.x, TEST_TOLERANCE);
    }
    {
        const result = Incomplete.secondKind(Dual, .promote(-2.010142159165531), .promote(0.9997758655364375));
        try std.testing.expectApproxEqAbs(-1.0960437630027704, result.x, TEST_TOLERANCE);
    }
}

test "incomplete third-kind" {
    const Dual = ad.DualNumber(f64, 1);
    {
        const result = Incomplete.thirdKind(Dual, .derivative(0.7), .promote(0.3), .promote(0.1));
        try std.testing.expectApproxEqAbs(0.3068870557197901, result.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.009562543919461083, result.dx[0], TEST_TOLERANCE);
    }
    {
        const result = Incomplete.thirdKind(
            Dual,
            .derivative(1.3327492505862528),
            .promote(2.09375611),
            .promote(0.99977587),
        );
        try std.testing.expectApproxEqAbs(-36.950167910312736, result.x, TEST_TOLERANCE);
    }
    {
        const result = Incomplete.thirdKind(
            Dual,
            .derivative(1.3327492505862528),
            .promote(-2.09375611),
            .promote(0.99977587),
        );
        try std.testing.expectApproxEqAbs(36.950167910312736, result.x, TEST_TOLERANCE);
    }
    {
        const result = Incomplete.thirdKind(
            Dual,
            .derivative(5.761904761904757),
            .promote(0.3),
            .promote(0.7),
        );
        try std.testing.expectApproxEqAbs(0.380478334936009, result.x, TEST_TOLERANCE);
    }
    {
        const result = Incomplete.thirdKind(
            Dual,
            .derivative(1.4074536408094236),
            .promote(-1.002572032700238),
            .promote(0.45950796050720055),
        );
        try std.testing.expectApproxEqAbs(-8.422822080313027, result.x, TEST_TOLERANCE);
    }
    {
        const result = Incomplete.thirdKind(
            Dual,
            .derivative(10.4074536408094236),
            .promote(-1.002572032700238),
            .promote(0.45950796050720055),
        );
        try std.testing.expectApproxEqAbs(-0.04489936124242888, result.x, TEST_TOLERANCE);
    }
    {
        const result = Incomplete.thirdKind(
            Dual,
            .derivative(10.4074536408094236),
            .promote(1.002572032700238),
            .promote(0.45950796050720055),
        );
        try std.testing.expectApproxEqAbs(0.04489936124242888, result.x, TEST_TOLERANCE);
    }
}
