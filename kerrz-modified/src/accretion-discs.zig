const std = @import("std");
const ad = @import("zad");
const rootsolve = @import("rootsolve");
const options = @import("options");
const geometry = @import("geometry.zig");
const geodesic = @import("geodesic.zig");
const orbits = @import("orbits.zig");
const iterators = @import("iterators.zig");

const FourVector = geometry.FourVector;
const NullGeodesic = geodesic.NullGeodesic;
const TraceResult = geodesic.TraceResult;
const TracingConfig = geodesic.TracingConfig;
const KerrMetric = geometry.KerrMetric;

const TEST_TOLERANCE = options.test_numerical_tolerance;

/// A utilty disc that can be used to represent a sphere of constant radius.
pub fn SphereDisc(comptime T: type) type {
    return struct {
        const Self = @This();
        /// The radius of the sphere in units of `rg`.
        radius: T = .promote(10.0),

        /// Trace the geodesic. In this case, the winding parameter is an upper
        /// limit.
        pub fn trace(
            self: Self,
            metric: KerrMetric(T),
            geod: NullGeodesic(T),
            config: TracingConfig,
        ) TraceResult(T) {
            const A = T.Algebra;

            var modified_config = config;
            modified_config.winding = 0;

            const builder = geod.traceBuilder(metric, modified_config);

            switch (builder.state.radial_case) {
                .case_II => {
                    const r1 = builder.state.radial_roots.r1.real();
                    const r3 = builder.state.radial_roots.r3.real();
                    const r4 = builder.state.radial_roots.r4.real();

                    const r31 = A.sub(r3, r1);
                    const r41 = A.sub(r4, r1);

                    const x2_squared = A.mult(
                        A.div(r31, r41),
                        A.div(A.sub(self.radius, r4), A.sub(self.radius, r3)),
                    );

                    if (x2_squared.x < 0 or x2_squared.x >= 1) {
                        // No solution for the disc exists.
                        return builder.atMinoTime(.zero);
                    }
                },
                else => {},
            }

            const mino_time = builder.state.minoTimeToRadius(
                self.radius,
                @floatFromInt(geod.radial_sign),
            );

            var result = builder.atMinoTime(mino_time);

            switch (result.status) {
                .no_status, .intersected_disc => {
                    if (std.math.approxEqAbs(T.T, result.r.x, self.radius.x, 1e-4)) {
                        result.status = .intersected_disc;
                    }
                },
                .event_horizon, .infinity => {},
            }
            return result;
        }
    };
}

pub fn ShakuraSunyaev(comptime T: type) type {
    return struct {
        const Self = @This();
        /// The accretion rate in units of Eddington.
        eddington_fraction: T = .promote(0.2),

        state: struct {
            /// The radiative efficiency. Determined by the ISCO.
            radiative_efficiency: T,
            /// Should be set to the ISCO, where the radiative efficiency is also set.
            inner_radius: T,

            pub fn init(metric: KerrMetric(T)) @This() {
                const A = T.Algebra;
                const e_am = orbits.keplerianEnergyAngularMomentum(T, metric, metric.isco);
                return .{
                    .radiative_efficiency = A.sub(.one, e_am.energy),
                    .inner_radius = metric.isco,
                };
            }
        },

        const Solver = ThickDiscSolver(T, ShakuraSunyaev);

        pub fn trace(
            self: Self,
            metric: KerrMetric(T),
            geod: NullGeodesic(T),
            config: TracingConfig,
        ) TraceResult(T) {
            return Solver.solveIntersection(self, metric, geod, config);
        }

        pub fn heightAtRadius(self: Self, radius: T) ?T {
            const A = T.Algebra;
            if (radius.x < self.state.inner_radius.x) {
                return null;
            }
            const factor = A.sub(.one, A.div(self.state.inner_radius, radius));
            return A.mult(
                .promote(3),
                A.div(
                    A.mult(self.eddington_fraction, factor),
                    self.state.radiative_efficiency,
                ),
            );
        }

        pub fn adapt(self: Self, comptime NewT: type) ShakuraSunyaev(NewT) {
            return .{
                .eddington_fraction = .adaptFrom(self.eddington_fraction),
                .state = .{
                    .radiative_efficiency = .adaptFrom(self.state.radiative_efficiency),
                    .inner_radius = .adaptFrom(self.state.inner_radius),
                },
            };
        }
    };
}

/// The standard, equatorial, geometrically (infinitely) thin accretion discs.
pub fn ThinDisc(comptime T: type) type {
    return struct {
        const Self = @This();
        /// The inner radius of the accretion disc in units of `rg`.
        inner_radius: T = .promote(0),
        /// The outer radius of the accretion disc in units of `rg`.
        outer_radius: T = .promote(30.0),

        /// Trace the geodesic. In this case, the winding parameter is an upper
        /// limit.
        pub fn trace(
            self: Self,
            metric: KerrMetric(T),
            geod: NullGeodesic(T),
            config: TracingConfig,
        ) TraceResult(T) {
            const angle: T = .promote(std.math.pi / 2.0);
            var conf = config;

            const max_winding: usize = @min(config.winding, 1);

            for (0..max_winding + 1) |winding| {
                conf.winding = winding;
                var result = geod.traceToAngle(metric, angle, conf);

                switch (result.status) {
                    .event_horizon, .infinity => return result,
                    else => {},
                }

                if (result.r.x <= self.outer_radius.x and result.r.x >= self.inner_radius.x) {
                    result.status = .intersected_disc;
                    return result;
                }

                if (winding == max_winding) {
                    if (result.status == .no_status) {
                        result.status = .infinity;
                    }
                    return result;
                }
            }

            unreachable;
        }
    };
}

test "thin disc" {
    const Dual = ad.DualNumber(f64, 0);
    const FV = FourVector(Dual);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    const x: FV = .{
        .t = .zero,
        .r = .promote(1e8),
        .th = .promote(std.math.degreesToRadians(20)),
        .ph = .zero,
    };
    const geod = NullGeodesic(Dual).fromImpactParameters(
        metric,
        x,
        .promote(2.0),
        .promote(2.0),
    );

    const result = geod.traceDisc(metric, .{ .thin_disc = .{} }, .{ .winding = 0 });
    try std.testing.expectEqual(geodesic.Status.intersected_disc, result.status);
}

/// A pseudo-geometry, representing an infinite plane that is parallel to the
/// equatorial plane but lifted out to some `height`.
pub fn DatumPlane(comptime T: type) type {
    return struct {
        const Self = @This();

        // Internally used for the root solving. Requires one extra slot for
        // the root-solving.
        const Dual1 = ad.DualNumber(T.T, T.N + 1);

        /// The height (in `rg`) that the plane is above the equatorial plane.
        height: Dual1 = .zero,
        /// The inner radius of the datum (disc) plane.
        inner_radius: Dual1 = .zero,
        /// The outer radius of the datum (disc) plane.
        outer_radius: Dual1 = .promote(std.math.floatMax(T.T)),

        /// Trace the geodesic until it intersects with the datum plane.
        ///
        /// The algorithm for finding intersections with the datum plane is as
        /// follows:
        ///
        ///                                       b
        ///                                       / Geodesic
        ///                                      /
        ///         Datum plane -> -------------/c-------------------
        ///                                   /
        ///    Black hole ->  O ------------/-------------- <- equatorial plane
        ///                                 a
        ///
        /// The geodesic is traced to the mid plane intersection at `a`. At
        /// that point, the geodesic height as a function of Mino time `h(tau)
        /// - height` is bracketed by the Mino time at the equatorial plane and
        /// 0. A root solving algorithm is then used to find the `tau` that
        /// gives an intersection with the plane, up to some tolerance.
        ///
        /// The exception is for geodesics that would hit the event horizon.
        /// These will sometimes not intersect the equatorial plane before
        /// falling into the black hole, in which case the maximal Mino time is
        /// found for the photon to reach some `r + ε` outside of the event
        /// horizon radius. The intersection point is then a bracketed root as
        /// before.
        ///
        /// For the case where the origin of the geodesic is below the datum
        /// plane (i.e. the datum plane is a ceiling), an `unreachable` occurs,
        /// as this case has not yet been handled. In theory, the bracketing
        /// tau could be found by increasing the Mino time monotonically a few
        /// times until some sufficiently large value is found (in which case
        /// the trace is terminated) or a bracket is determined.
        pub fn trace(
            self: Self,
            metric: KerrMetric(T),
            geod: NullGeodesic(T),
            config: TracingConfig,
        ) TraceResult(T) {
            const A = Dual1.Algebra;
            const angle: Dual1 = .promote(std.math.pi / 2.0);

            const metric_1 = metric.adapt(Dual1);
            const geod_1 = geod.adapt(Dual1);

            // Overwrite the maximum winding for datum plane traces, as the
            // winding number is determined afterwards.
            var modified_config = config;
            modified_config.winding = 0;

            // The winding number is here fixed to zero, since the algorithm
            // currently implemented for datum intersection does not include
            // any geodesic winding.
            // TODO: allow for non-zero winding.
            var result = geod_1.traceToAngle(metric_1, angle, modified_config);
            // TODO: avoid retracing and remove these kinds of things
            const reg_builder = geod.traceBuilder(metric, modified_config);
            var reg_result = reg_builder.atMinoTime(
                result.mino_time.popSlot(),
            );

            const ctx: RootSolveContext = .{
                .trace_builder = .fromResult(geod_1, metric_1, result),
                .height = self.height,
            };

            switch (result.status) {
                .event_horizon => {
                    // TODO: find Mino-time using another means and avoid
                    // retracing
                    result = ctx.trace_builder.atMinoTime(.promote(1.0));
                },
                else => {},
            }

            var init_tau: Dual1 = .promote(result.mino_time.x * 0.01);
            const pivot = ctx.trace_builder.atMinoTime(init_tau);

            var m_pivot = A.sub(measureHeight(
                Dual1,
                pivot.r,
                pivot.theta,
            ), self.height);

            const m_init = A.sub(measureHeight(
                Dual1,
                geod_1.x_init.r,
                geod_1.x_init.th,
            ), self.height);

            const m_final = A.sub(measureHeight(
                Dual1,
                result.r,
                result.theta,
            ), self.height);

            // If not bracketting, then there can be no intersection.
            if (std.math.sign(m_pivot.x) == std.math.sign(m_final.x)) {
                // Try the origin
                m_pivot = m_init;
                init_tau = .zero;
                if (std.math.sign(m_pivot.x) == std.math.sign(m_final.x)) {
                    if (reg_result.status == .no_status) {
                        reg_result.status = .infinity;
                    }
                    return reg_result;
                }
            }

            const sol = rootsolve.univariateSolveDual(
                Dual1,
                ctx,
                intersectionMeasure,
                A.mult(.promote(0.2), result.mino_time),
                .{
                    .lower_bound = init_tau,
                    .upper_bound = result.mino_time,
                    .value_at_lower_bound = m_pivot,
                    .value_at_upper_bound = m_final,
                },
            ) catch return reg_result;

            // TODO: this strips derivative information... I would ideally need
            // a derivative on the root, i.e. the univariateSolve should return
            // a dual number with one slot less, so that the derivatives in all
            // the other slots are maintained at that particular mino time.
            var new_result = reg_builder.atMinoTime(sol.x);

            // Check if within the datum plane boundary
            const proj_r = new_result.r.x * @sin(new_result.theta.x);
            if (proj_r > self.outer_radius.x or proj_r < self.inner_radius.x) {
                // Restore the original status
                new_result.status = result.status;
                if (new_result.status == .no_status) {
                    new_result.status = .infinity;
                }
                return new_result;
            }

            // Determine the status of the geodesic
            if (self.height.x < metric.horizon_radius.x) {
                // Calculate the `x` intersection of the horizon with the datum
                // plane
                const rs = metric.horizon_radius.x;
                const h = self.height.x;
                const rs_x = @sqrt(rs * rs - h * h);
                if (new_result.r.x * @sin(new_result.theta.x) < rs_x) {
                    new_result.status = .event_horizon;
                    return new_result;
                }
            }

            // Update the status codes
            switch (new_result.status) {
                .no_status, .event_horizon => {
                    new_result.status = .intersected_disc;
                },
                else => {},
            }

            return new_result;
        }

        const RootSolveContext = struct {
            trace_builder: NullGeodesic(Dual1).PathBuilder,
            height: Dual1,
        };

        fn intersectionMeasure(x: Dual1, ctx: RootSolveContext) Dual1 {
            const A = Dual1.Algebra;
            const result = ctx.trace_builder.atMinoTime(x);
            const h = measureHeight(Dual1, result.r, result.theta);
            return A.sub(ctx.height, h);
        }
    };
}

fn measureHeight(comptime T: type, r: T, theta: T) T {
    const A = T.Algebra;
    return A.abs(A.mult(r, A.cos(theta)));
}

fn testDatumPlane(alpha: f64, beta: f64, height: f64) !FourVector(ad.DualNumber(f64, 0)) {
    const Dual = ad.DualNumber(f64, 0);
    const FV = FourVector(Dual);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    const x: FV = .{
        .t = .zero,
        .r = .promote(1e6),
        .th = .promote(std.math.degreesToRadians(70)),
        .ph = .zero,
    };

    const geod = NullGeodesic(Dual).fromImpactParameters(
        metric,
        x,
        .promote(alpha),
        .promote(beta),
    );
    const result = geod.traceDisc(
        metric,
        .{ .datum_plane = .{ .height = .promote(height) } },
        .{ .winding = 0 },
    );

    const total = result.totalAntiderivatives(metric, geod);

    return .{
        .t = total.coordinateTime(metric, geod),
        .r = result.r,
        .th = result.theta,
        .ph = total.coordinateAzimuth(metric, geod),
    };
}

test "datum plane coordinates" {
    const r1 = try testDatumPlane(2.0, -5.0, 0.001);
    try std.testing.expectApproxEqAbs(14.682487364709566, r1.r.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(1.5707282184439986, r1.th.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(0.13067546292730106, r1.ph.x, TEST_TOLERANCE);

    const r2 = try testDatumPlane(2.0, 5.0, 0.001);
    try std.testing.expectApproxEqAbs(3.173939183622281, r2.r.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(1.5704812575539306, r2.th.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(2.6011044333982234, r2.ph.x, TEST_TOLERANCE);
}

/// A general union of all accretion disc types.
pub fn AccretionDisc(comptime T: type) type {
    return union(enum) {
        const Self = @This();

        datum_plane: DatumPlane(T),
        thin_disc: ThinDisc(T),
        equatorial_plane: void,
        sphere: SphereDisc(T),
        shakura_sunyaev: ShakuraSunyaev(T),

        /// Trace to the accretion disc.
        pub fn trace(
            self: Self,
            metric: KerrMetric(T),
            geod: NullGeodesic(T),
            config: TracingConfig,
        ) TraceResult(T) {
            switch (self) {
                .equatorial_plane => {
                    var result = geod.traceToAngle(
                        metric,
                        .promote(std.math.pi / 2.0),
                        .{},
                    );
                    switch (result.status) {
                        .event_horizon, .infinity => {},
                        // Everything else intersected with the infinite plane.
                        else => result.status = .intersected_disc,
                    }
                    return result;
                },
                inline else => |disc| return disc.trace(metric, geod, config),
            }
        }

        /// Adapt to a different dual number type.
        pub fn adapt(self: Self, comptime NewT: type) AccretionDisc(NewT) {
            return switch (self) {
                .datum_plane => |d| .{ .datum_plane = .{
                    .height = .adaptFrom(d.height),
                    .inner_radius = .adaptFrom(d.inner_radius),
                    .outer_radius = .adaptFrom(d.outer_radius),
                } },
                .thin_disc => |d| .{ .thin_disc = .{
                    .inner_radius = .adaptFrom(d.inner_radius),
                    .outer_radius = .adaptFrom(d.outer_radius),
                } },
                .equatorial_plane => .{
                    .equatorial_plane = {},
                },
                .sphere => |d| .{
                    .sphere = .{ .radius = .adaptFrom(d.radius) },
                },
                .shakura_sunyaev => |d| .{
                    .shakura_sunyaev = d.adapt(NewT),
                },
            };
        }
    };
}

/// Solves intersection points with thick disc geometries.
///
/// The `DiscType` must be the function that returns the disc type, so that it
/// may be instantiated for different dual number types needed during the
/// solve.
///
/// The `DiscType` must additionally define the following functions:
///
///      /// Calculate the disc height at a given projected radius. Should
///      /// return `null` only if the disc does not exist at the given
///      /// projection radius.
///      fn heightAtRadius(self: DiscType(T), proj_radius: T) ?T
///
fn ThickDiscSolver(comptime T: type, comptime DiscType: fn (comptime type) type) type {
    return struct {
        const Self = @This();
        const Dual1 = T.PushSlot();
        const A = Dual1.Algebra;

        builder: NullGeodesic(Dual1).PathBuilder,
        disc: DiscType(Dual1),
        m: KerrMetric(Dual1),
        finished: bool = false,
        last_trace: ?TraceResult(Dual1) = null,

        /// Calculate the difference between the projected height of
        /// the current geodesic onto the spin axis and the scale
        /// height of the accretion disc.
        ///
        /// This function will return null if the geodesic falls into
        /// the black hole (it will also set the `finished` flag of the
        /// context). It will also return null if the `heightAtRadius`
        /// function return null, indicating that no disc is present at
        /// that radius.
        fn objectiveNull(x: Dual1, this: *@This()) ?Dual1 {
            const trace = this.builder.atMinoTime(x);
            this.last_trace = trace;

            // A safety catch:
            if (trace.r.x <= this.m.horizon_radius.x) {
                this.finished = true;
                return null;
            }

            const r_projected = A.mult(trace.r, A.sin(trace.theta));
            const h_projected = A.mult(trace.r, A.cos(trace.theta));

            const h = this.disc.adapt(Dual1).heightAtRadius(r_projected) orelse
                return null;

            // Sanity check:
            std.debug.assert(std.math.isFinite(h.x));

            return A.sub(h, h_projected);
        }

        /// An alternative wrapper that will never return `null`.
        fn objective(x: Dual1, this: *@This()) Dual1 {
            return objectiveNull(x, this) orelse Dual1.one.neg();
        }

        fn zeroMinoTime(builder: NullGeodesic(T).PathBuilder) TraceResult(T) {
            var result = builder.atMinoTime(.zero);
            result.status = .event_horizon;
            return result;
        }

        /// The algorithm is a simple root solve after a bracketing interval
        /// has been found. Finding the bracketing interval is tricky,
        /// especially given some of the quirks of the problem. Consider the
        /// following cross-section comic:
        ///
        ///                  |    | x B
        ///                  |.C  |
        ///                  |    |
        ///     -------.     |    |
        ///       x A   ------ BH -----
        ///
        /// The vertical bars denote the innermost radius of the disc, between
        /// which `heightAtRadius` will return null. `BH` is the black hole,
        /// and the dashed line shows a segment of the disc scale height above
        /// the equatorial plane. The two labelled `x` are two points along
        /// some geodesic that bound the intersection root. The problem is if
        /// the `univariateSolve` call then takes a step C at any point, it
        /// will get a null height and convergence will fail. The the easiest
        /// way to avoid that is to take very small Mino time steps (which is
        /// currently what happens).
        ///
        /// Small Mino time steps are taken until the difference measure has a
        /// sign flip, and which point the bracket is given to the rootsolver
        /// to finish.
        fn solveIntersectionImpl(
            disc: DiscType(T),
            metric: KerrMetric(T),
            geod: NullGeodesic(T),
            config: TracingConfig,
        ) !TraceResult(T) {
            // Type promotions
            const ddisc = disc.adapt(Dual1);
            const dgeod = geod.adapt(Dual1);

            const builder = dgeod.traceBuilder(metric.adapt(Dual1), config);

            var self: Self = .{
                .builder = builder,
                .disc = ddisc,
                .m = metric.adapt(Dual1),
            };

            // Trace to the equatorial plane as the starting point
            const midplane_result = builder.traceToAngle(.promote(std.math.pi / 2.0));

            var max_mino_time: Dual1 = .promote(5);

            if (midplane_result.mino_time.x > 0 and midplane_result.r.x > 20.0) {
                max_mino_time = midplane_result.mino_time;
            }

            // Init diff should never be null. If the height function returns
            // null, that is treated as some positive distance outside of the
            // disc.
            var curr_diff: Dual1 = objectiveNull(.zero, &self) orelse
                Dual1.one.neg();

            var left_mino: Dual1 = .zero;
            var right_mino: Dual1 = .zero;
            var right_diff: Dual1 = .zero;

            // TODO: make this maximum step configurable.
            for (1..200) |step| {
                var mino_time = Dual1.promote(config.mino_step * @as(Dual1.T, @floatFromInt(step)));

                if (mino_time.x > max_mino_time.x) {
                    mino_time = max_mino_time;
                }

                const diff = objectiveNull(mino_time, &self) orelse {
                    if (self.finished) {
                        // Geodesic went into the black hole
                        var result = self.last_trace.?.adapt(T);
                        result.status = .event_horizon;
                        return result;
                    }
                    continue;
                };

                // There's always a chance:
                if (@abs(diff.x) < config.error_tolerance) {
                    var result = self.last_trace.?;
                    result.status = .intersected_disc;
                    return result.adapt(T);
                }

                if (std.math.sign(diff.x) != std.math.sign(curr_diff.x)) {
                    // Bracketing interval found.
                    right_mino = mino_time;
                    right_diff = diff;
                    break;
                }

                left_mino = mino_time;
                curr_diff = diff;

                if (mino_time.x == max_mino_time.x) break;
            }

            if (left_mino.x >= right_mino.x) {
                // No solution found across all of the steps.
                var result = self.last_trace.?.adapt(T);
                result.status = .infinity;
                return result;
            }

            // The guess is the midpoint between the left and right bracket,
            // i.e. an uninformed guess.
            const guess = A.mult(.promote(0.5), A.add(left_mino, right_mino));
            const sol = try rootsolve.univariateSolveDual(
                Dual1,
                &self,
                objective,
                guess,
                .{
                    .error_tolerance = config.error_tolerance,
                    .lower_bound = left_mino,
                    .upper_bound = right_mino,
                    .value_at_lower_bound = curr_diff,
                    .value_at_upper_bound = right_diff,
                },
            );

            // TODO: do something with this information.
            _ = sol;

            var result = self.last_trace.?.adapt(T);
            result.status = .intersected_disc;
            return result;
        }

        /// Find the first (as in, earliest Mino time) intersection of a null
        /// geodesic with a thick accretion disc.
        pub fn solveIntersection(
            disc: DiscType(T),
            metric: KerrMetric(T),
            geod: NullGeodesic(T),
            config: TracingConfig,
        ) TraceResult(T) {
            var modified_config = config;
            modified_config.winding = 0;
            return solveIntersectionImpl(disc, metric, geod, modified_config) catch {
                std.debug.print(" >>> Failed to converge\n", .{});
                const builder = geod.traceBuilder(metric, modified_config);
                return zeroMinoTime(builder);
            };
        }
    };
}
