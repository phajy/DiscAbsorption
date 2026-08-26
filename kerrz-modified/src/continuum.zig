const std = @import("std");
const ad = @import("zad");
const zfits = @import("zfits");
const rootsolve = @import("rootsolve");
const options = @import("options");

const utils = @import("utils.zig");

const TEST_TOLERANCE = options.test_numerical_tolerance;

const interpolations = @import("interpolations.zig");
const emissivity = @import("emissivity.zig");
const redshift = @import("redshift.zig");
const orbits = @import("orbits.zig");
const iterators = @import("iterators.zig");
const geometry = @import("geometry.zig");
const geodesic = @import("geodesic.zig");
const accretion_discs = @import("accretion-discs.zig");
const transfer_tables = @import("transfer-tables.zig");
const TransferFunction = @import("tools.zig").TransferFunction;

const FourVector = geometry.FourVector;
const NullGeodesic = geodesic.NullGeodesic;
const TraceResult = geodesic.TraceResult;
const KerrMetric = geometry.KerrMetric;
const DatumPlane = accretion_discs.DatumPlane;
const CoronalModel = emissivity.CoronalModel;

pub fn ContinuumGeodesic(comptime T: type) type {
    return struct {
        const Self = @This();
        /// The geodesic itself.
        geod: NullGeodesic(T),
        /// The computed result at the endpoint of the geodesic, which is here
        /// at the observer.
        result: TraceResult(T),
        /// How many function calls were required to evaluate the result.
        f_calls: usize,
        /// The absolute error on the position (in rg).
        err: T.T,
        /// The coordinate time, with the observer distance subtracted. This is
        /// pre-calculated since it is often needed.
        delta_t: T,
        /// The local angle on the sky from the perspective of the corona. This
        /// is the initial angle as measured by the local point of the geodesic
        /// travelling to the observer.
        local_theta: T,
        /// The Jacobian term |∂θ / ∂Y| of the solution, where Y is the
        /// angle on the local sky.
        jac: T.T,
        /// The impact parameters of this geodesic.
        alpha: T,
        beta: T,

        /// Calculate |∂ cosθ / ∂ cosY| from the known Jacobian value.
        pub fn jacobianAsCosine(self: Self) T.T {
            const sin_t = @sin(self.result.theta.x);
            const sin_d = @sin(self.local_theta.x);
            return (sin_t / sin_d) * self.jac;
        }

        /// Calculate the energyshift of this geodesic as measured by the
        /// observer.
        pub fn calculateEnergyshift(
            self: Self,
            metric: KerrMetric(T),
            corona: CoronalModel(T),
        ) T {
            std.debug.assert(corona == .lamppost);
            // A safety check:
            std.debug.assert(std.math.approxEqAbs(
                T.T,
                corona.lamppost.x.r.x,
                self.geod.x_init.r.x,
                TEST_TOLERANCE,
            ));

            const ts_obs = metric.tangentSpaceAlt(self.result.r, self.result.theta);
            const v_obs = orbits.stationary(T, ts_obs);

            // Need the energyshift at corona, sensitive to its velocity:
            return redshift.redshift(
                T,
                ts_obs,
                v_obs,
                self.result.velocity(metric, self.geod),
                corona.lamppost.ts,
                corona.lamppost.v,
                self.geod.initialVelocity(metric),
            );
        }
    };
}

pub const Error = error{InvalidForModel};

/// Calculate the corona-to-observer time, weighted by the flux to obtain the
/// 'centroid' time. That is, it returns a single time that can be used to
/// estimate when the majority of the emission would arrive at the observer.
///
/// The allocator is needed for temporary allocations, and all memory will be
/// freed before this function returns.
pub fn continuumGeodesic(
    comptime T: type,
    metric: KerrMetric(T),
    x_obs: FourVector(T),
    model: CoronalModel(T),
) !ContinuumGeodesic(T) {
    return try switch (model) {
        .lamppost => |corona| continuumGeodesic_lamppost(T, metric, x_obs, corona),
        .disc, .ring, .umbrella => return Error.InvalidForModel,
    };
}

/// Reverse the direction of the velocity vector to point oppositely.
fn reverseVelocity(comptime T: type, vel: FourVector(T)) FourVector(T) {
    return .{
        .t = vel.t,
        .r = vel.r.neg(),
        .th = vel.th.neg(),
        .ph = vel.ph,
    };
}

/// For the lamppost, there is complete axis-symmetry. So the easiest way to
/// compute a photon from the observer to the corona is to do so in reverse,
/// and to find the local angle `theta` that results in a photon at infinity
/// with the desired observer inclination.
fn continuumGeodesic_lamppost(
    comptime T: type,
    metric: KerrMetric(T),
    x_obs: FourVector(T),
    model: emissivity.Lamppost(T),
) !ContinuumGeodesic(T) {
    const A = T.Algebra;
    const Dual1 = T.PushSlot();

    const metric1 = metric.adapt(Dual1);
    const x_obs1 = x_obs.adapt(Dual1);

    const sol = try traceFromOnAxisPoint(T, metric, x_obs, model);
    const total = sol.result.totalAntiderivatives(metric1, sol.geod);
    const ip = sol.toImpactParameters(Dual1, metric1, x_obs1);

    const jac = @abs(sol.result.theta.dx[0]);

    return .{
        .err = sol.err,
        // Plus one for the above re-evaluation.
        .f_calls = sol.f_calls + 1,
        .delta_t = A.sub(total.coordinateTime(metric1, sol.geod).popSlot(), x_obs.r),
        .result = sol.result.adapt(T),
        .geod = sol.geod.adapt(T),
        .jac = jac,
        .local_theta = .promote(sol.angle),
        .alpha = .adaptFrom(ip.alpha),
        .beta = .adaptFrom(ip.beta),
    };
}

test "lensing factors" {
    const T = ad.DualNumber(f64, 0);
    const metric = KerrMetric(T).init(.one, .promote(0.998));
    {
        const x_obs: FourVector(T) = .{
            .t = .zero,
            .r = .promote(1e8),
            .th = .promote(std.math.degreesToRadians(45)),
            .ph = .zero,
        };

        const lamppost = emissivity.Lamppost(T).init(
            metric,
            .{ .height = .promote(10.0) },
        );

        const sol = try continuumGeodesic_lamppost(T, metric, x_obs, lamppost);
        try std.testing.expectApproxEqAbs(
            0,
            sol.err,
            TEST_TOLERANCE,
        );
        // Check the Jacobian values:
        try std.testing.expectApproxEqAbs(
            1.135347170876976,
            sol.jac,
            TEST_TOLERANCE,
        );
        try std.testing.expectApproxEqAbs(
            0.802072375641173,
            1.0 / sol.jacobianAsCosine(),
            TEST_TOLERANCE,
        );

        try std.testing.expectApproxEqAbs(
            0.6996020272321579,
            sol.local_theta.x,
            TEST_TOLERANCE,
        );
    }
}

fn testLensingTime(dist: f64, incl: f64) !f64 {
    const T = ad.DualNumber(f64, 0);
    const metric = KerrMetric(T).init(.one, .promote(0.998));
    const lamppost = emissivity.Lamppost(T).init(
        metric,
        .{ .height = .promote(10.0) },
    );
    const x_obs: FourVector(T) = .{
        .t = .zero,
        .r = .promote(dist),
        .th = .promote(std.math.degreesToRadians(incl)),
        .ph = .zero,
    };

    const sol = try continuumGeodesic_lamppost(T, metric, x_obs, lamppost);
    try std.testing.expectApproxEqAbs(
        0,
        sol.err,
        TEST_TOLERANCE,
    );

    const total = sol.result.totalAntiderivatives(metric, sol.geod);

    return total.coordinateTime(metric, sol.geod).x;
}

test "lensing times" {
    try std.testing.expectApproxEqAbs(
        10000020.585304461,
        try testLensingTime(1e7, 45),
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        10000015.843675144,
        try testLensingTime(1e7, 2),
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        100000026.80693872,
        try testLensingTime(1e8, 88),
        TEST_TOLERANCE,
    );

    try std.testing.expectApproxEqAbs(
        100011.62907058328,
        try testLensingTime(1e5, 45),
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        100010.0096396974,
        try testLensingTime(1e5, 30.0),
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        100018.70602544083,
        try testLensingTime(1e5, 88),
        TEST_TOLERANCE,
    );
}

test "lamppost there and back again" {
    const T = ad.DualNumber(f64, 0);
    const metric = KerrMetric(T).init(.one, .promote(0.998));
    {
        const x_obs: FourVector(T) = .{
            .t = .zero,
            .r = .promote(1e3),
            .th = .promote(std.math.degreesToRadians(66)),
            .ph = .zero,
        };

        const lamppost = emissivity.Lamppost(T).init(
            metric,
            .{ .height = .promote(3.3) },
        );

        const sol = try continuumGeodesic_lamppost(T, metric, x_obs, lamppost);

        // Check that the solution is valid.
        try std.testing.expectApproxEqAbs(
            0,
            sol.err,
            TEST_TOLERANCE,
        );

        // Check that the coordinate time is correct, c.f. Gradus.jl: Technically
        // Gradus.jl gives 1011.733... but that's good enough for me.
        try std.testing.expectApproxEqAbs(
            1011.3505072828794 - x_obs.r.x,
            sol.delta_t.x,
            TEST_TOLERANCE,
        );

        // And that integrating backwards gives the right result:
        const geod = NullGeodesic(T).fromConstantsOfMotion(
            x_obs,
            sol.geod.E,
            sol.geod.L,
            sol.geod.Q,
        );
        const builder = geod.traceBuilder(metric, .{});
        const res = builder.atMinoTime(sol.result.mino_time);
        try std.testing.expectApproxEqAbs(
            lamppost.x.r.x,
            res.r.x,
            TEST_TOLERANCE,
        );

        // Also test the same can be achieved via the impact parameters:
        const geod_impacts = NullGeodesic(T).fromImpactParameters(
            metric,
            x_obs,
            sol.alpha,
            sol.beta,
        );
        const builder_impacts = geod_impacts.traceBuilder(metric, .{});
        const res_impacts = builder_impacts.atMinoTime(sol.result.mino_time);
        try std.testing.expectApproxEqAbs(
            lamppost.x.r.x,
            res_impacts.r.x,
            TEST_TOLERANCE,
        );

        // And that the local angles match up correctly:
        const angles = res_impacts.localAngles(metric, geod_impacts, lamppost.v);
        try std.testing.expectApproxEqAbs(
            @abs(sol.local_theta.x - std.math.pi),
            angles.theta.x,
            TEST_TOLERANCE,
        );
    }

    {
        const x_obs: FourVector(T) = .{
            .t = .zero,
            .r = .promote(1e8),
            .th = .promote(std.math.degreesToRadians(80)),
            .ph = .zero,
        };

        const lamppost = emissivity.Lamppost(T).init(
            metric,
            .{ .height = .promote(2.2) },
        );

        const sol = try traceFromOnAxisPoint(T, metric, x_obs, lamppost);

        // Check that the solution is valid.
        try std.testing.expectApproxEqAbs(
            0,
            sol.err,
            TEST_TOLERANCE,
        );

        // And by impact parameter:
        const ip = sol.toImpactParameters(T, metric, x_obs);
        const geod2 = NullGeodesic(T).fromImpactParameters(metric, x_obs, ip.alpha, ip.beta);
        const builder2 = geod2.traceBuilder(metric, .{});
        const res2 = builder2.atMinoTime(sol.result.mino_time.popSlot());
        try std.testing.expectApproxEqAbs(
            lamppost.x.r.x,
            res2.r.x,
            TEST_TOLERANCE,
        );
    }
}

fn GeodesicAndResult(comptime Dual: type) type {
    return struct {
        const Self = @This();
        geod: NullGeodesic(Dual),
        result: TraceResult(Dual),
        err: Dual.T,
        f_calls: usize,
        /// The polar angle on the local sky of the emitting point.
        angle: Dual.T,

        /// Calculate the impact parameters on the image plane corresponding to
        /// this solution.
        pub fn toImpactParameters(
            self: Self,
            comptime T: type,
            metric: KerrMetric(T),
            x_obs: FourVector(T),
        ) geodesic.ImpactParameters(T) {
            // And that integrating backwards gives the right result:
            const v_final = self.result.velocity(metric.adapt(Dual), self.geod).adapt(T);
            const geod = NullGeodesic(T).fromVelocity(
                metric,
                x_obs,
                reverseVelocity(T, v_final),
            );
            return geod.calculateImpactParameters(metric, x_obs.th);
        }

        pub fn adapt(self: *const Self, comptime T: type) GeodesicAndResult(T) {
            return .{
                .geod = self.geod.adapt(T),
                .result = self.result.adapt(T),
                .err = self.err,
                .f_calls = self.f_calls,
                .angle = self.angle,
            };
        }
    };
}

/// Find a geodesic that travels from a point on the spin axis to the
/// observer's inclination.
///
/// For convenience, this currently uses the lamppost structure, as it
/// calculates the tangent space and so on.
///
/// TODO: make it not use the lamppost structure.
fn traceFromOnAxisPoint(
    comptime T: type,
    metric: KerrMetric(T),
    x_obs: FourVector(T),
    model: emissivity.Lamppost(T),
) !GeodesicAndResult(T.PushSlot()) {
    const Dual1 = T.PushSlot();
    const A = Dual1.Algebra;

    const SolverContext = struct {
        m: KerrMetric(Dual1),
        x: FourVector(Dual1),
        lp: emissivity.Lamppost(Dual1),

        last_result: ?TraceResult(Dual1) = null,
        last_geod: ?NullGeodesic(Dual1) = null,

        fn findInclination(theta: Dual1, ctx: *@This()) Dual1 {
            const geod = NullGeodesic(Dual1).fromSkyAnglesTangentSpace(
                ctx.m,
                ctx.lp.ts,
                ctx.lp.v,
                theta,
                .zero,
            );

            const result = geod.traceToRadius(
                ctx.m,
                ctx.x.r,
                .{},
            );

            ctx.last_geod = geod;
            ctx.last_result = result;

            // TODO: mod 2 pi or some such?
            return A.sub(result.theta, ctx.x.th);
        }
    };

    var ctx: SolverContext = .{
        .m = metric.adapt(Dual1),
        .x = x_obs.adapt(Dual1),
        .lp = model.adapt(Dual1),
    };

    const sol = try rootsolve.univariateSolveDual(
        Dual1,
        &ctx,
        SolverContext.findInclination,
        .promote(std.math.degreesToRadians(1.0)),
        .{
            .lower_bound = .zero,
            .fallback = .bisect,
            .error_tolerance = 1e-9,
            .dx_tolerance = 1e-12,
        },
    );

    return .{
        .geod = ctx.last_geod.?,
        .result = ctx.last_result.?,
        .err = sol.err.x,
        .f_calls = sol.f_calls,
        .angle = sol.x.x,
    };
}

pub const ContinuumTransferOptions = struct {
    /// The maximum number of continuum traces to calculate.
    max_points: usize = 100,
};

pub fn InterpolatedContinuumPoint(comptime T: type) type {
    return struct {
        const Self = @This();

        /// The observer-to-corona time with the observer distance subtracted.
        delta_t: T,
        /// The observer-to-corona Mino time.
        mino_time: T,
        /// The observed redshift.
        g: T,
        /// The azimuthal coordinate along the ring.
        phi: T,

        /// The local polar angle on the sky of the corona.
        local_theta: T,
        /// The local azimuthal angle on the sky of the corona.
        local_phi: T,

        /// The impact parameters themselves.
        alpha: T,
        beta: T,

        // The lensing factor. See the docstring of `ContinuumTrace`.
        jac: T,

        fn interpolateBetween(self: Self, other: Self, weight: T) Self {
            var out: Self = undefined;
            inline for (@typeInfo(Self).@"struct".fields) |field| {
                @field(out, field.name) = interpolations.lerpValue(
                    T,
                    weight,
                    @field(self, field.name),
                    @field(other, field.name),
                );
            }
            return out;
        }
    };
}

pub fn ContinuumTrace(comptime T: type) type {
    return struct {
        const Self = @This();
        /// The observer-to-corona time with the observer distance subtracted.
        delta_t: T,
        /// The observer-to-corona Mino time.
        mino_time: T,
        /// The observed redshift.
        g: T,
        /// The azimuthal coordinate along the ring.
        phi: T,

        /// The local polar angle on the sky of the corona.
        local_theta: T,
        /// The local azimuthal angle on the sky of the corona.
        local_phi: T,

        /// The impact parameters themselves.
        alpha: T,
        beta: T,

        /// The solid angle Jacobian determinant, closely related to the
        /// lensing factor. That is:
        ///
        ///     J = ∂θ / ∂Y
        ///
        /// which can be used for mapping flux bundles as `dA / dΩ`, after
        /// multiplying by relevant `sin` factors.
        jac: T.T = 0,

        /// The geodesic itself.
        geod: NullGeodesic(T),
        /// The computed result at the endpoint of the geodesic, which is here
        /// at the observer.
        result: TraceResult(T),

        /// The total number of geodesic calls needed to solve and evaluate
        /// this continuum trace.
        f_calls: usize,

        fn toPoint(self: Self) InterpolatedContinuumPoint(T.T) {
            return .{
                .delta_t = self.delta_t.x,
                .mino_time = self.mino_time.x,
                .g = self.g.x,
                .phi = self.phi.x,
                .local_theta = self.local_theta.x,
                .local_phi = self.local_phi.x,
                .alpha = self.alpha.x,
                .beta = self.beta.x,
                .jac = self.jac,
            };
        }

        /// Interpolate between two traced with a given weight.
        pub fn interpolateBetween(
            self: Self,
            other: Self,
            weight: T.T,
        ) InterpolatedContinuumPoint(T.T) {
            return self.toPoint().interpolateBetween(other.toPoint(), weight);
        }
    };
}

fn ContinuumTransferTracer(comptime T: type) type {
    return struct {
        const Self = @This();
        x_obs: FourVector(T),
        ts_obs: geometry.KerrMetric(T).TangentSpace,
        v_obs: FourVector(T),
        model: CoronalModel(T),
        metric: KerrMetric(T),

        fn init(
            metric: KerrMetric(T),
            x_obs: FourVector(T),
            model: CoronalModel(T),
        ) Self {
            const ts_obs = metric.tangentSpace(x_obs);
            const v_obs = orbits.stationary(T, ts_obs);

            return .{
                .ts_obs = ts_obs,
                .v_obs = v_obs,
                .model = model,
                .metric = metric,
                .x_obs = x_obs,
            };
        }

        fn retraceForRing(
            self: *const Self,
            lp_geod: NullGeodesic(T),
            lp_result: TraceResult(T),
            alpha: T,
            beta: T,
            f_calls: usize,
        ) ContinuumTrace(T) {
            const Dual = T.PushSlot();
            const A = T.Algebra;

            const v_medium = self.model.ring.v;
            const v_medium2 = v_medium.adapt(Dual);

            // The redshift can be calculated from the original lamppost-style
            // trace using the velocity profile of the ring.
            const g = redshift.redshift(
                T,
                self.ts_obs,
                self.v_obs,
                lp_geod.initialVelocity(self.metric),
                self.model.ring.ts,
                v_medium,
                lp_result.velocity(self.metric, lp_geod),
            );

            const metric2 = self.metric.adapt(Dual);

            // The strategy here is to find the local angles in the frame of
            // the co-rotating ring, and to then use them to calculate the
            // lensing Jacobian back to the observer.
            const local_angles = lp_result.localAnglesReverse(self.metric, lp_geod, v_medium);

            const lp_total = lp_result.totalAntiderivatives(self.metric, lp_geod);
            const lp_phi = lp_total.coordinateAzimuth(self.metric, lp_geod);

            var d_theta = Dual.adaptFrom(local_angles.theta);
            const d_phi = Dual.adaptFrom(local_angles.phi);

            // Setup the derivative slots:
            d_theta.dx[0] = 1;

            // Setup the new geodesic to trace.
            var geod = NullGeodesic(Dual).fromSkyAnglesTangentSpace(
                metric2,
                // Ignore the actual azimuthal coordinate since everything is
                // basically axisymmetric.
                // TODO: is this alright? It might introduce an error R /
                // R_obs?
                self.model.ring.ts.adapt(Dual),
                v_medium2,
                d_theta,
                d_phi,
            );

            var config: geodesic.TracingConfig = .{};

            if (geod.radial_sign < 0) {
                config.r_winding = 1;
            }

            // Trace to the observer radius:
            const result = geod.traceToRadius(
                metric2,
                .adaptFrom(self.x_obs.r),
                config,
            );

            // Sanity check, we should have recovered the original observer
            // position:
            const error_tolerance = 1e-3;
            std.debug.assert(std.math.approxEqAbs(
                T.T,
                result.theta.x,
                self.x_obs.th.x,
                error_tolerance,
            ));

            // TODO: it would seem the full Jacobian term isn't needed, and
            // blows up at the poles if the azimuthal (`obs_phi`) part is
            // included?
            const jac = result.theta.dx[0];

            const delta_t = A.sub(lp_total.coordinateTime(self.metric, lp_geod), self.x_obs.r);

            return .{
                .g = g,
                .delta_t = delta_t,
                .mino_time = .adaptFrom(result.mino_time),
                .phi = lp_phi,
                .local_theta = local_angles.theta,
                .local_phi = local_angles.phi,
                .alpha = alpha,
                .beta = beta,
                .jac = jac,
                .geod = geod.adapt(T),
                .result = result.adapt(T),
                // +1 for the retrace
                .f_calls = f_calls + 1,
            };
        }
    };
}

pub fn ContinuumTransfer(comptime T: type) type {
    return struct {
        const Self = @This();

        theta: []const T.T,
        traces: []ContinuumTrace(T),
        // Centroid impact parameters
        alpha: T,
        beta: T,
        /// The centroid time difference after subtracting off the observer's radial
        /// position
        delta_t: T,
        r_obs: T,
        theta_obs: T,

        /// Free the allocated resources.
        pub fn deinit(self: *const Self, allocator: std.mem.Allocator) void {
            allocator.free(self.theta);
            allocator.free(self.traces);
        }

        fn tracePhi(t: ContinuumTrace(T)) T.T {
            return t.phi.x;
        }

        /// Interpolate a trace at a given azimuthal coordinate along the ring.
        pub fn interpolateAzimuth(self: *const Self, phi: T.T) InterpolatedContinuumPoint(T.T) {
            const info = interpolations.periodicInterpolate(
                ContinuumTrace(T),
                self.traces,
                phi,
                tracePhi,
            );
            var t1 = self.traces[info.left];
            if (info.left_mod) {
                t1.phi.x -= 2 * std.math.pi;
            }
            return t1.interpolateBetween(self.traces[info.right], info.weight);
        }

        /// Serialise to a FITS table.
        pub fn toFITS(self: *const Self, allocator: std.mem.Allocator) !zfits.Hdu {
            var hdu = zfits.Hdu.init(allocator, .{ .binary_table = .empty });
            errdefer hdu.deinit();

            try hdu.setName(
                "CONTINUUM",
                "This is a continuum transfer function table",
            );
            try utils.addKerrzFITSInfo(&hdu);

            // Write the centroid information
            try hdu.appendHeaderRecord("R_OBS", .{
                .comment = "The distance to the observer",
                .value = .{ .float = @floatCast(self.r_obs.x) },
            });
            try hdu.appendHeaderRecord("TH_OBS", .{
                .comment = "The observer inclination",
                .value = .{ .float = @floatCast(self.theta_obs.x) },
            });
            try hdu.appendHeaderRecord("C_DTIME", .{
                .comment = "The centroid time difference from r_obs",
                .value = .{ .float = @floatCast(self.delta_t.x - self.r_obs.x) },
            });
            try hdu.appendHeaderRecord("C_ALPHA", .{
                .comment = "The centroid alpha impact parameter",
                .value = .{ .float = @floatCast(self.alpha.x) },
            });
            try hdu.appendHeaderRecord("C_BETA", .{
                .comment = "The centroid beta impact parameter",
                .value = .{ .float = @floatCast(self.beta.x) },
            });

            // Setup the columns
            try hdu.data.binary_table.appendColumn(.{
                .label = "theta",
                .comment = "Angle on the observer's image plane",
                .units = "radians",
            });
            try hdu.data.binary_table.appendColumn(.{
                .label = "alpha",
                .comment = "The x-axis impact parameter",
                .units = "rg",
            });
            try hdu.data.binary_table.appendColumn(.{
                .label = "beta",
                .comment = "The y-axis impact parameter",
                .units = "rg",
            });
            try hdu.data.binary_table.appendColumn(.{
                .label = "mino_time",
                .comment = "Total elapsed Mino time",
            });
            try hdu.data.binary_table.appendColumn(.{
                .label = "deltime",
                .comment = "Time difference from `r_obs`",
                .units = "tg",
            });
            try hdu.data.binary_table.appendColumn(.{
                .label = "g",
                .comment = "Energyshift of this photon",
            });
            try hdu.data.binary_table.appendColumn(.{
                .label = "phi",
                .comment = "Azimuthal coordinate on the ring",
                .units = "radians",
            });
            try hdu.data.binary_table.appendColumn(.{
                .label = "localth",
                .comment = "Local polar angle on the sky",
                .units = "radians",
            });
            try hdu.data.binary_table.appendColumn(.{
                .label = "localph",
                .comment = "Local azimuthal angle on the sky",
                .units = "radians",
            });
            try hdu.data.binary_table.appendColumn(.{
                .label = "lensing",
                .comment = "Lensing factor |dth/dY|",
            });

            for (self.theta, self.traces) |th, t| {
                const row = try hdu.data.binary_table.addRow();
                row.cols[0].one.float_32 = @floatCast(th);
                row.cols[1].one.float_32 = @floatCast(t.alpha.x);
                row.cols[2].one.float_32 = @floatCast(t.beta.x);
                row.cols[3].one.float_32 = @floatCast(t.mino_time.x);
                row.cols[4].one.float_32 = @floatCast(t.delta_t.x);
                row.cols[5].one.float_32 = @floatCast(t.g.x);
                row.cols[6].one.float_32 = @floatCast(t.phi.x);
                row.cols[7].one.float_32 = @floatCast(t.local_theta.x);
                row.cols[8].one.float_32 = @floatCast(t.local_phi.x);
                row.cols[9].one.float_32 = @floatCast(t.jac);
            }

            return hdu;
        }
    };
}

/// Compute a transfer function for a given coronal geometry.
pub fn transferFunction(
    comptime T: type,
    allocator: std.mem.Allocator,
    metric: KerrMetric(T),
    x_obs: FourVector(T),
    model: CoronalModel(T),
    opts: ContinuumTransferOptions,
) !ContinuumTransfer(T) {
    const Dual1 = T.PushSlot();
    const error_tolerance = 1e-3;
    if (model != .ring) unreachable;

    const metric1 = metric.adapt(Dual1);
    const x_obs1 = x_obs.adapt(Dual1);

    const r_target = model.ring.getRadius();
    const h_target = model.ring.getHeight();
    const disc: accretion_discs.AccretionDisc(T) = .{
        .datum_plane = .{
            .height = h_target.pushSlot(),
        },
    };

    // Calculate where to put the zero point of the image projection for the
    // root solver. This is needed because if the height of the disc is
    // significant, then the (alpha = 0, beta = 0) point on the image plane is
    // no longer contained within the projected ring and the solver will fail.
    //
    // The easiest way to determine this is to calculate the impact paramters
    // for the lamppost and then center on that:
    const lp = try traceFromOnAxisPoint(
        T,
        metric,
        x_obs,
        .init(metric, .{ .height = h_target }),
    );

    // Now convert the solution to a set of impact parameters.
    const lp_impact_parameters = lp.toImpactParameters(
        Dual1,
        metric1,
        x_obs1,
    );

    const tracer: ContinuumTransferTracer(T) = .init(metric, x_obs, model);

    const traces = try allocator.alloc(ContinuumTrace(T), opts.max_points);
    errdefer allocator.free(traces);
    var theta_itt = iterators.RangeIterator(T.T).init(0.0, 2 * std.math.pi, traces.len);
    const thetas = try theta_itt.drain(allocator);
    errdefer allocator.free(thetas);

    for (traces, thetas) |*trace, theta| {
        const image_theta: Dual1 = .promote(theta + error_tolerance);

        const sol = try transfer_tables.impactOffsetForRadius(
            Dual1,
            metric1,
            x_obs1,
            image_theta,
            r_target.pushSlot(),
            .{
                .disc = disc.adapt(Dual1),
                .beta_offset = lp_impact_parameters.beta.x,
            },
        );

        var impact_params = transfer_tables.toImpactParameters(
            Dual1,
            .promote(sol.r_offset),
            .promote(sol.theta_image_plane),
        );
        // Apply the zero offset.
        impact_params.beta = Dual1.Algebra.add(
            impact_params.beta,
            lp_impact_parameters.beta,
        );

        trace.* = tracer.retraceForRing(
            sol.geod.adapt(T),
            sol.result.adapt(T),
            impact_params.alpha.popSlot(),
            impact_params.beta.popSlot(),
            sol.f_calls,
        );
    }

    // For the centroid geodesic
    const total = lp.result.totalAntiderivatives(metric1, lp.geod);
    const ip = lp.toImpactParameters(Dual1, metric1, x_obs1);

    return .{
        .traces = traces,
        .theta = thetas,
        .alpha = .adaptFrom(ip.alpha),
        .beta = .adaptFrom(ip.beta),
        .delta_t = .adaptFrom(total.coordinateTime(metric1, lp.geod)),
        .r_obs = x_obs.r,
        .theta_obs = x_obs.th,
    };
}

test "edgecases" {
    const T = ad.DualNumber(f64, 0);
    const metric = KerrMetric(T).init(.one, .promote(0.998));
    {
        const x_obs: FourVector(T) = .{
            .t = .zero,
            .r = .promote(1e5),
            .th = .promote(0.001),
            .ph = .zero,
        };

        const lamppost = emissivity.Lamppost(T).init(
            metric,
            .{ .height = .promote(6) },
        );

        const sol = try continuumGeodesic_lamppost(T, metric, x_obs, lamppost);
        try std.testing.expectApproxEqAbs(
            0,
            sol.err,
            TEST_TOLERANCE,
        );
        // Check the Jacobian values:
        try std.testing.expectApproxEqAbs(
            1.2165087833577988,
            sol.jac,
            TEST_TOLERANCE,
        );
        try std.testing.expectApproxEqAbs(
            0.6758048148393498,
            1.0 / sol.jacobianAsCosine(),
            TEST_TOLERANCE,
        );

        try std.testing.expectApproxEqAbs(
            0.0008218598093330716,
            sol.local_theta.x,
            TEST_TOLERANCE,
        );
    }
}
