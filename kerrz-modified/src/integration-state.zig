const std = @import("std");
const root = @import("root.zig");
const elliptic_integrals = @import("elliptic-integrals.zig");
const potentials = @import("potentials.zig");
const geodesic = @import("geodesic.zig");
const antiderivatives = @import("antiderivatives.zig");
const tracy = @import("tracy.zig");

const KerrMetric = root.KerrMetric;

const NullGeodesic = geodesic.NullGeodesic;

const RadialRoots = potentials.RadialRoots;
const AngularRoots = potentials.AngularRoots;
const RadialCase = potentials.RadialCase;
const AngularCase = potentials.AngularCase;

const RadialCaseCache = antiderivatives.RadialCaseCache;
const AngularCaseCache = antiderivatives.AngularCaseCache;

const TEST_TOLERANCE = @import("options").test_numerical_tolerance;

/// All of the anti-derivatives, potentials and roots needed to calculate
/// geodesic motion. These are only the things that depend on the initial
/// conditions, and are not e.g. functions of the Mino time.
///
/// This structure is lazy and calculates new potentials only as they are
/// needed. As such, avoid field access and use the accessor methods if
/// accessing values for the first time.
pub fn IntegrationState(comptime T: type) type {
    return struct {
        const Self = @This();
        /// Tracks the radian winding sign.
        ///
        /// TODO: this is currently a hack on top of the
        /// ring-corona-to-observer traces for the case where the radial sign
        /// is positive, and a radial turning point would be encountered. It
        /// had previously been giving negative Mino time, so adjusting this
        /// sign fixes it to give the correct coordinates (modulo time and
        /// azimuth, as those are handled elsewhere).
        r_winding_sign: T.T = 1,

        /// The radial roots of the geodesic.
        radial_roots: RadialRoots(T),

        /// The angular roots of the geodesic.
        angular_roots: AngularRoots(T),

        /// The radial case, as determined by the radial roots.
        radial_case: RadialCase,

        /// The angular case, as determined by the sign of the Carter's
        /// constant.
        angular_case: AngularCase,

        /// Cached values related to calcuating radial motion for the various
        /// cases.
        radial_cache: RadialCaseCache(T),

        /// Cached values related to calcuating angular motion for the various
        /// cases.
        angular_cache: AngularCaseCache(T),

        /// Adapt the dual number type.
        pub fn adapt(self: Self, comptime NewT: type) IntegrationState(NewT) {
            return .{
                .radial_roots = self.radial_roots.adapt(NewT),
                .angular_roots = self.angular_roots.adapt(NewT),
                .radial_case = self.radial_case,
                .angular_case = self.angular_case,
                .radial_cache = self.radial_cache.adapt(NewT),
                .angular_cache = self.angular_cache.adapt(NewT),
            };
        }

        /// Initialise the integrator state. Calculate all quantities related
        /// to the start point of the geodesic and are needed for all possible
        /// motion calculations.
        pub fn init(metric: KerrMetric(T), geod: NullGeodesic(T)) Self {
            var ctx = tracy.trace(@src());
            defer ctx.end();

            const eta_sign: i32 = @intFromFloat(std.math.sign(geod.eta.x));

            const angular_roots = potentials.rootsOfAngularPotential(
                T,
                metric,
                geod.eta,
                geod.lambda,
            );

            const angular_case = AngularCase.fromSign(eta_sign);

            const angular_cache = antiderivatives.angularCache(
                T,
                metric,
                angular_roots,
                angular_case,
                geod.x_init.th,
            );

            const radial_roots = potentials.rootsOfRadialPotential(
                T,
                metric,
                geod.eta,
                geod.lambda,
            );

            const radial_case = radial_roots.determineCase(
                metric.horizon_radius,
                geod.x_init.r,
            );

            const radial_cache = antiderivatives.radialCache(
                T,
                radial_case,
                radial_roots,
                geod.x_init.r,
            );

            return .{
                .radial_roots = radial_roots,
                .angular_roots = angular_roots,
                .radial_case = radial_case,
                .angular_case = angular_case,
                .radial_cache = radial_cache,
                .angular_cache = angular_cache,
            };
        }

        /// Calculate the Mino time until both (real) angular turning points
        /// have been reached.
        pub fn minoTimeToAngularTurns(
            self: Self,
            winding: usize,
            theta_sign: T.T,
        ) struct { tau_0: T, tau_1: T } {
            const A = T.Algebra;

            var theta_0 = T.zero;
            var theta_1 = T.zero;

            switch (self.angular_case) {
                .normal => {
                    // For normal motion, turning points are theta 1 and theta
                    // 4 Equation (28)
                    theta_0 = A.acos(A.sqrt(self.angular_roots.u_plus));
                    theta_1 = A.acos(A.sqrt(self.angular_roots.u_plus).neg());
                },
                .vortical => {
                    // For vortical motion, turning points depend on the
                    // hemisphere Equation (55)
                    theta_0 = A.acos(A.sqrt(self.angular_roots.u_minus));
                    theta_1 = A.acos(A.sqrt(self.angular_roots.u_plus));
                },
            }

            var mino_time_0 = self.minoTimeToAngle(
                theta_0,
                winding,
                theta_sign,
            );

            var mino_time_1 = self.minoTimeToAngle(
                theta_1,
                winding,
                theta_sign,
            );

            if (self.angular_cache.G_theta_half_libration.x < mino_time_0.x) {
                mino_time_0.x = @mod(mino_time_0.x, self.angular_cache.G_theta_half_libration.x);
            }
            if (self.angular_cache.G_theta_half_libration.x < mino_time_1.x) {
                mino_time_1.x = @mod(mino_time_1.x, self.angular_cache.G_theta_half_libration.x);
            }

            return .{
                .tau_0 = mino_time_0,
                .tau_1 = mino_time_1,
            };
        }

        /// Calculate the Mino time until the radial turning points
        /// have been reached.
        ///
        /// This method will only calculate the mino times for the real roots
        /// exterior of the horizon, which depends on the particular radial
        /// case.
        ///
        /// If no time could be calculated, returns null for that root.
        pub fn minoTimeToRadialTurns(
            self: Self,
            r_sign: T.T,
        ) struct { r_0: ?T, r_1: ?T } {
            switch (self.radial_case) {
                .case_I => { // r1 < r2 <= r <= r3 < r4
                    const r2 = self.radial_roots.r2.real();
                    const r3 = self.radial_roots.r3.real();

                    return .{
                        .r_0 = self.minoTimeToRadius(r2, r_sign),
                        .r_1 = self.minoTimeToRadius(r3, r_sign),
                    };
                },
                .case_II => { // r1 < r2 < r3 < r4 <= r
                    const r4 = self.radial_roots.r4.real();
                    return .{
                        .r_0 = self.minoTimeToRadius(r4, r_sign),
                        .r_1 = null,
                    };
                },
                .case_III, .case_IV => { // All roots interior
                    return .{ .r_0 = null, .r_1 = null };
                },
            }
        }

        /// Calculate the Mino time to a particular angle for `n` windings.
        pub fn minoTimeToAngle(
            self: Self,
            angle: T,
            n: usize,
            theta_sign: T.T,
        ) T {
            const A = T.Algebra;
            const angular_values = antiderivatives.angularFromCache(
                T,
                self.angular_case,
                self.angular_cache,
                angle,
            );

            // G&L Equations (A9) or (A10)
            const n_correction = if (theta_sign < 0) n + 1 else n;
            const total_sign: T.T = if (@rem(n_correction, 2) == 0) 1 else -1;

            // G&L Equation (52)
            return A.add(
                A.mult(self.angular_cache.G_theta_half_libration, .promote(@floatFromInt(n_correction))),
                A.mult(
                    .promote(theta_sign),
                    A.sub(
                        A.mult(.promote(total_sign), angular_values.G_theta_final),
                        self.angular_cache.G_theta_init,
                    ),
                ),
            );
        }

        /// Calculate the Mino time to a particular radius for `n` windings.
        pub fn minoTimeToRadius(
            self: Self,
            radius: T,
            r_sign: T.T,
        ) T {
            const A = T.Algebra;
            const radial_values = antiderivatives.radialCache(
                T,
                self.radial_case,
                self.radial_roots,
                radius,
            );

            // Equation (B7) and also Equation (B25)
            return A.mult(
                .promote(r_sign),
                A.sub(
                    A.mult(.promote(self.r_winding_sign), radial_values.I_0),
                    self.radial_cache.I_0,
                ),
            );
        }

        /// Calculate the radial coordinate at a particular Mino time.
        pub fn radiusAtMinoTime(
            self: Self,
            mino_time: T,
            r_sign: T.T,
        ) antiderivatives.PrincipalRadialValues(T) {
            return antiderivatives.radialValuesFromCache(
                T,
                self.radial_case,
                self.radial_cache,
                mino_time,
                r_sign,
            );
        }

        /// Calculate the poloidal anglular coordinate at a particular Mino
        /// time.
        pub fn angleAtMinoTime(
            self: Self,
            mino_time: T,
            theta_sign: T.T,
        ) antiderivatives.PrincipalAngularValues(T) {
            return antiderivatives.angularValues(
                T,
                self.angular_case,
                self.angular_cache,
                mino_time,
                theta_sign,
            );
        }
    };
}
