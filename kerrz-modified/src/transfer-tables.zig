/// This file is for computing results towards Cunningham transfer function
/// tables and the like.
const std = @import("std");
const ad = @import("zad");
const zfits = @import("zfits");
const rootsolve = @import("rootsolve");
const dinterp = @import("dinterp");
const options = @import("options");

const utils = @import("utils.zig");

const redshift = @import("redshift.zig");
const orbits = @import("orbits.zig");
const interpolations = @import("interpolations.zig");
const spectra = @import("spectra.zig");
const emissivity = @import("emissivity.zig");
const accretion_discs = @import("accretion-discs.zig");
const polarisation = @import("polarisation.zig");

const geometry = @import("geometry.zig");
const geodesic = @import("geodesic.zig");
const iterators = @import("iterators.zig");

const Matrix = @import("matrix.zig").Matrix;
const EmissivityProfile = emissivity.EmissivityProfile;
const RangeIterator = iterators.RangeIterator;
const Grid = iterators.Grid;
const DualNumber = ad.DualNumber;
const KerrMetric = geometry.KerrMetric;
const FourVector = geometry.FourVector;
const NullGeodesic = geodesic.NullGeodesic;
const PolarisationConstant = polarisation.PolarisationConstant;

const TEST_TOLERANCE = options.test_numerical_tolerance;

/// This controls how close to the edges in g_star space to treat as stable.
/// Very close to `gstar = 0` or `1`, the Jacobians diverge, which can skew the
/// rest of the transfer function integration.
///
/// This value has been somewhat empirically determined to balance accurate
/// integration against the refinement of the `gstar` grid.
///
/// Tweaking
/// it can significantly clean up the line profiles.
const EDGE_TOLERANCE = 1.2e-4;

const logger = std.log.scoped(.transfer_tables);

fn gstarIterator(comptime T: type, size: usize) RangeIterator(T) {
    return gstarIteratorTol(T, size, EDGE_TOLERANCE);
}

fn gstarIteratorTol(comptime T: type, size: usize, h: T) RangeIterator(T) {
    return .init(h, 1 - h, size);
}

fn angularDelta(comptime T: type, a1: T, a2: T) T {
    const two_pi = std.math.pi * 2.0;
    // Normalize inputs to `[0, 2π)`.
    const _a1 = @mod(a1, two_pi);
    const _a2 = @mod(a2, two_pi);
    // Map to `[-π,π]` relative difference.
    return @mod(_a2 - _a1 + std.math.pi, two_pi) - std.math.pi;
}

fn angularMidpoint(comptime T: type, a1: T, a2: T) T {
    const two_pi = std.math.pi * 2.0;
    const _a1 = @mod(a1, two_pi);
    const delta = angularDelta(T, a1, a2);
    // Normalize to `[0, 2π)`.
    return @mod(_a1 + delta / 2.0, two_pi);
}

fn rotateAngle(angle: anytype) @TypeOf(angle) {
    return @mod(angle - std.math.pi / 2.0, std.math.pi * 2.0);
}

fn GoldenSectionOptions(comptime T: type) type {
    return struct {
        optimise: enum { min, max },
        lower_bound: T,
        upper_bound: T,
        max_steps: usize = 20,
    };
}

fn optimiseGoldenSection(
    comptime T: type,
    ctx: anytype,
    comptime score_func: fn (@TypeOf(ctx), T) T,
    opts: GoldenSectionOptions(T),
) T {
    const inv_phi = (@sqrt(5.0) - 1.0) / 2.0;

    var a = opts.lower_bound;
    var c = opts.upper_bound;
    var b = a + inv_phi * (c - a);

    var f_calls: usize = 1;
    var saved_y = score_func(ctx, b);

    while (a != c and f_calls < opts.max_steps) {
        const x = a + inv_phi * (b - a);
        const y = score_func(ctx, x);
        f_calls += 1;

        if (switch (opts.optimise) {
            .min => y < saved_y,
            .max => y > saved_y,
        }) {
            saved_y = y;
            c = b;
            b = x;
        } else {
            a = c;
            c = x;
        }
    }

    return b;
}

pub fn ImpactOffsetResult(comptime T: type) type {
    return struct {
        theta_image_plane: T.T,
        r_offset: T.T,
        /// The mino time to intersection.
        mino_time: T,
        err: T.T,
        f_calls: usize,

        /// The geodesic result that solved the impact offset problem.
        result: geodesic.TraceResult(T),
        /// The geodesic itself that solved the impact offset problem.
        geod: geodesic.NullGeodesic(T),
    };
}

pub fn ImpactOffsetOptions(comptime T: type) type {
    return struct {
        initial_guess: T.T = 4.0,
        /// How much to offset the zero-point on the beta axis.
        /// For the algorithm to work, the zero point must be within the
        /// projection, since it does not allow negative radii.
        beta_offset: T.T = 0,
        disc: accretion_discs.AccretionDisc(T) = .equatorial_plane,
    };
}

/// Translate `(r, theta)` on the image plane to a set of impact parameters.
pub fn toImpactParameters(
    comptime T: type,
    r: T,
    theta: T,
) geodesic.ImpactParameters(T) {
    const A = T.Algebra;
    const alpha = A.mult(A.cos(theta), r);
    const beta = A.mult(A.sin(theta), r);
    return .{
        .alpha = alpha,
        .beta = beta,
    };
}

/// Calculate the impact parameter offset for a particular angle on the image
/// plane that gives a set of impact parameters that hit a particular radius on
/// the accretion disc. That is, it finds `R` given `θ` that can be used to
/// recover impact parameters α and β that map a geodesic to a particular
/// annulus on the disc. See the schematic below
///
///     +-------------------------+
///     |                         | Here, the `O` marks the origin of the image
///     |                (α, β)   | plane.
///     |               x         |
///     |              /.         | The algorithm root-solves along the line from
///     |            R/ .         | the origin of the image plane inclined at θ
///     |            /  .         | off the α = 0 line (x-axis) in the clockwise
///     |           /   .         | direction. The angle marked is therefore the
///     |          / 2π - θ       | reciprocal, 2π - θ.
///     |         O .....         |
///     |         origin          |
///     |                         |
///     +-------------------------+
///
/// This function currently cannot compute derivatives, and the geodesic must
/// be re-traced to obtain any e.g. Jacobian terms.
///
/// This function also only considers the upper half of the disc, and does not
/// include the contributions from false images.
pub fn impactOffsetForRadius(
    comptime T: type,
    metric: KerrMetric(T),
    x_obs: FourVector(T),
    theta_image_plane: T,
    target_radius: T,
    opts: ImpactOffsetOptions(T),
) !ImpactOffsetResult(T) {
    const Ctx = TargetRadiusContext(T);
    return try impactOffsetForRadiusImpl(
        T,
        metric,
        x_obs,
        theta_image_plane,
        target_radius,
        opts,
        Ctx.impactOffsetForRadius_target,
    );
}

fn impactOffsetForRadiusImpl(
    comptime T: type,
    metric: KerrMetric(T),
    x_obs: FourVector(T),
    theta_image_plane: T,
    target_radius: T,
    opts: ImpactOffsetOptions(T),
    comptime optimFunc: anytype,
) !ImpactOffsetResult(T) {
    const error_tolerance = 1e-6;
    const Ctx = TargetRadiusContext(T);
    if (T.N < 1) {
        @compileError("Must call impactOffsetForRadius using dual numbers with at least one slot");
    }
    var ctx: Ctx = .{
        .metric = metric,
        .x_obs = x_obs,
        .theta_image_plane = theta_image_plane,
        .target_radius = target_radius,
        .disc = opts.disc,
        .beta_offset = .promote(opts.beta_offset),
    };

    const sol = try rootsolve.univariateSolveDual(
        T,
        &ctx,
        optimFunc,
        .promote(opts.initial_guess),
        .{
            .lower_bound = .promote(@sqrt(error_tolerance)),
            .fallback = .bisect,
            .error_tolerance = error_tolerance,
            .dx_tolerance = error_tolerance,
        },
    );

    return .{
        .r_offset = sol.x.x,
        .err = sol.err.x,
        .theta_image_plane = theta_image_plane.x,
        .f_calls = sol.f_calls,
        .mino_time = ctx.last_result.?.mino_time,
        .result = ctx.last_result.?,
        .geod = ctx.last_geodesic.?,
    };
}

fn TargetRadiusContext(comptime T: type) type {
    const A = T.Algebra;
    return struct {
        const Self = @This();

        metric: KerrMetric(T),
        x_obs: FourVector(T),
        theta_image_plane: T,
        target_radius: T,
        disc: accretion_discs.AccretionDisc(T),
        beta_offset: T = .zero,

        last_result: ?geodesic.TraceResult(T) = null,
        last_geodesic: ?geodesic.NullGeodesic(T) = null,

        /// The target function used by the rootsolve call in `impactOffsetForRadius`.
        fn impactOffsetForRadius_target(r: T, ctx: *Self) T {
            const p = toImpactParameters(T, r, ctx.theta_image_plane);

            const geod = NullGeodesic(T).fromImpactParameters(
                ctx.metric,
                ctx.x_obs,
                p.alpha,
                A.add(p.beta, ctx.beta_offset),
            );

            const result = ctx.disc.trace(ctx.metric, geod, .{});
            // Project the result into the equatorial plane.
            const r_projected = A.mult(result.r, A.sin(result.theta));

            // Store this trace:
            ctx.last_result = result;
            ctx.last_geodesic = geod;

            return A.sub(r_projected, ctx.target_radius);
        }
    };
}

const TestOffsetParameters = struct {
    spin: f64,
    theta_obs: f64,
    theta_image: f64,
    r_target: f64,
};

fn testOffset(params: TestOffsetParameters) !ImpactOffsetResult(DualNumber(f64, 1)) {
    const Dual = DualNumber(f64, 1);
    const metric: KerrMetric(Dual) = .init(.promote(1.0), .promote(params.spin));
    const x_obs: FourVector(Dual) = .{
        .t = .zero,
        .r = .promote(1e7),
        .th = .promote(std.math.degreesToRadians(params.theta_obs)),
        .ph = .zero,
    };
    const result = try impactOffsetForRadius(
        Dual,
        metric,
        x_obs,
        .promote(std.math.degreesToRadians(params.theta_image)),
        .promote(params.r_target),
        .{},
    );
    return result;
}

test "impact offset" {
    const result = try testOffset(.{
        .spin = 0.998,
        .r_target = 3.0,
        .theta_image = -30,
        .theta_obs = 78,
    });
    try std.testing.expectApproxEqAbs(
        0.0,
        result.err,
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        1.2830834052434665,
        result.r_offset,
        TEST_TOLERANCE,
    );
}

const IndexSet = struct {
    low1: usize,
    low2: usize,
    up1: usize,
    up2: usize,
};

/// A single instance of Cunningham's transfer function.
///
/// The lower branch is always between g_max_index and g_min_index going from
/// maximum to minimum energyshift (i.e. reversed), whereas the upper branch is
/// from g_min index to the end, and from 0 to g_max_index.
pub fn CunninghamTransferFunction(comptime T: type) type {
    return struct {
        const Self = @This();

        /// Convert g to gstar.
        pub fn gstarFrom(g: T.T, g_min: T.T, g_max: T.T) T.T {
            return (g - g_min) / (g_max - g_min);
        }

        /// Convert gstar to g.
        pub fn gFrom(gstar: T.T, g_min: T.T, g_max: T.T) T.T {
            return (g_max - g_min) * gstar + g_min;
        }

        /// Convert g to gstar. See also `gstarFrom` and `gFrom`.
        pub fn toGStar(self: Self, g: T.T) T.T {
            return gstarFrom(g, self.g_min, self.g_max);
        }

        /// A single point in the transfer function, and all associated values.
        pub const Trace = struct {
            /// The offset angle on the image plane.
            image_angle: T.T,
            /// The offset radius on the image plane. Together with `image_angle` they
            /// describe the impact parameters.
            image_radius: T.T,
            /// The alpha impact parameter.
            alpha: T.T,
            /// The beta impact parameter.
            beta: T.T,
            /// The redshift value of this trace.
            g: T.T,
            /// The normalised redshift value of this trace.
            g_star: T.T,
            /// The azimuthal coordinate of this trace.
            phi: T.T,
            /// The light-travel time of this trace, with the observer distance
            /// subtracted. That is, `t - r_obs`.
            delta_t: T.T,
            /// The |d(r, g)/d(alpha, beta)| Jacobian term, specifically the
            /// absolute value of the determinant.
            jacobian: T.T,
            /// Cunningham's actual transfer function, or the value thereof.
            f: T.T,
            /// The root solver error.
            r_err: T.T,
            /// The polarisation normalised Stokes X = Q/I parameter.
            stokes_x: T.T,
            /// The polarisation normalised Stokes Y = U/I parameter.
            stokes_y: T.T,
            /// The local theta (polar) angle between the geodesic and the disc
            /// normal.
            local_theta: T.T,
            /// The local phi (azimuthal) angle between the geodesic and the disc
            /// normal.
            local_phi: T.T,

            /// The zero trace, where every field is set to zero. This is used
            /// in de serialising a transfer function table.
            pub const zero = std.mem.zeroes(Trace);

            fn ascendingAngle(_: void, lhs: Trace, rhs: Trace) bool {
                // Here's a trick to make the branch separation easier later in
                // the Cunningham transfer functions.  When angle = 0, it
                // traces the right hand side of the image plane, and for angle
                // = π the left. But that's also likely where the extrema in
                // redshift will end up, so, reorient the plane so that angle =
                // 0 is offset somewhere else.
                return rotateAngle(lhs.image_angle) < rotateAngle(rhs.image_angle);
            }

            fn gStarMinPredicate(g_star: T.T, t: Trace) bool {
                return t.g_star > g_star;
            }

            fn gStarMaxPredicate(g_star: T.T, t: Trace) bool {
                return t.g_star < g_star;
            }

            fn interpolateBetween(self: Trace, other: Trace, weight: T.T) Trace {
                var out: Trace = undefined;
                inline for (@typeInfo(Trace).@"struct".fields) |field| {
                    @field(out, field.name) = interpolations.lerpValue(
                        T.T,
                        weight,
                        @field(self, field.name),
                        @field(other, field.name),
                    );
                }
                return out;
            }
        };

        /// The target radius on the accretion disc.
        target_radius: T.T,
        /// Traces calculated for this radius.
        traces: []Trace,
        /// The minimum energyshift for this annulus.
        g_min: T.T,
        /// The index of the minimum energy shift.
        g_min_index: usize,
        /// The maximum energyshift for this annulus.
        g_max: T.T,
        /// The index of the maximum energy shift.
        g_max_index: usize,

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.traces);
            self.* = undefined;
        }

        pub const InterpolatedBranches = struct {
            lower: T.T,
            upper: T.T,

            /// The lower branch light change in travel time.
            delta_t_lower: T.T,
            /// The upper branch light change in travel time.
            delta_t_upper: T.T,

            /// The local polar angle on the upper branch.
            theta_upper: T.T,
            /// The local polar angle on the lower branch.
            theta_lower: T.T,

            const zero: InterpolatedBranches = .{
                .lower = 0,
                .upper = 0,
                .delta_t_lower = 0,
                .delta_t_upper = 0,
                .theta_lower = 0,
                .theta_upper = 0,
            };
        };

        const WeightsAndIndices = struct {
            weight_lower: T.T,
            weight_upper: T.T,
            inds: IndexSet,
        };

        fn interpolationWeightsAndIndices(self: Self, g_star: T.T) ?WeightsAndIndices {
            // Calculations should be an order of magnitude better than what is
            // serialised. The tolerance here is used to avoid the common
            // pathologies at the extrema that are not easy to detect.
            const tolerance = EDGE_TOLERANCE / 10.0;
            if (g_star < tolerance or g_star > 1 - tolerance) return null;
            const inds = self.interpolationIndices(g_star);

            const weight_lower = interpolations.lerpWeight(
                T.T,
                g_star,
                self.traces[inds.low1].g_star,
                self.traces[inds.low2].g_star,
            );
            const weight_upper = interpolations.lerpWeight(
                T.T,
                g_star,
                self.traces[inds.up1].g_star,
                self.traces[inds.up2].g_star,
            );

            return .{
                .inds = inds,
                .weight_lower = weight_lower,
                .weight_upper = weight_upper,
            };
        }

        /// An interpolated Trace at a particular `g_star`. Gives both the
        /// trace on the upper and lower branch.
        pub const InterpolatedTrace = struct {
            upper: Trace,
            lower: Trace,
        };

        /// Interpolate every field of Trace at some `g`. Returns null if `g`
        /// is out of bounds.
        pub fn interpolateTrace(self: Self, g: T.T) ?InterpolatedTrace {
            if (g < self.g_min or g > self.g_max) return .zero;
            return self.interpolateTraceGstar(self.toGStar(g));
        }

        /// Same as `interpolateTrace` but the argument is `g_star` instead of
        /// `g`.
        ///
        /// Returns `null` if `g_star` is out of bounds.
        pub fn interpolateTraceGstar(self: Self, g_star: T.T) ?InterpolatedTrace {
            const wi = self.interpolationWeightsAndIndices(g_star) orelse
                return null;
            const lower = self.traces[wi.inds.low1].interpolateBetween(
                self.traces[wi.inds.low2],
                wi.weight_lower,
            );
            const upper = self.traces[wi.inds.up1].interpolateBetween(
                self.traces[wi.inds.up2],
                wi.weight_upper,
            );
            return .{ .upper = upper, .lower = lower };
        }

        /// Interpolate the value of `f` at some `g`. Returns zero if `g` is
        /// out of bounds.
        pub fn interpolateBranches(self: Self, g: T.T) InterpolatedBranches {
            if (g < self.g_min or g > self.g_max) return .zero;
            return self.interpolateBranchesGstar(self.toGStar(g));
        }

        /// Same as `interpolateBranches` but the argument is `g_star` instead
        /// of `g`.
        pub fn interpolateBranchesGstar(self: Self, g_star: T.T) InterpolatedBranches {
            const wi = self.interpolationWeightsAndIndices(g_star) orelse
                return .zero;

            return .{
                .lower = interpolations.lerpValue(
                    T.T,
                    wi.weight_lower,
                    self.traces[wi.inds.low1].f,
                    self.traces[wi.inds.low2].f,
                ),
                .upper = interpolations.lerpValue(
                    T.T,
                    wi.weight_upper,
                    self.traces[wi.inds.up1].f,
                    self.traces[wi.inds.up2].f,
                ),
                .delta_t_lower = interpolations.lerpValue(
                    T.T,
                    wi.weight_lower,
                    self.traces[wi.inds.low1].delta_t,
                    self.traces[wi.inds.low2].delta_t,
                ),
                .delta_t_upper = interpolations.lerpValue(
                    T.T,
                    wi.weight_upper,
                    self.traces[wi.inds.up1].delta_t,
                    self.traces[wi.inds.up2].delta_t,
                ),
                .theta_lower = interpolations.lerpValue(
                    T.T,
                    wi.weight_lower,
                    self.traces[wi.inds.low1].local_theta,
                    self.traces[wi.inds.low2].local_theta,
                ),
                .theta_upper = interpolations.lerpValue(
                    T.T,
                    wi.weight_upper,
                    self.traces[wi.inds.up1].local_theta,
                    self.traces[wi.inds.up2].local_theta,
                ),
            };
        }

        fn interpolationIndices(self: Self, g_star: T.T) IndexSet {
            std.debug.assert(self.g_max_index < self.g_min_index);

            const upper_index = std.sort.partitionPoint(
                Trace,
                self.traces[self.g_max_index..self.g_min_index],
                g_star,
                Trace.gStarMinPredicate,
            ) + self.g_max_index;

            const i = std.sort.partitionPoint(
                Trace,
                self.traces[0..self.g_max_index],
                g_star,
                Trace.gStarMaxPredicate,
            );

            if (i == 0) {
                // Check the other half
                const j = std.sort.partitionPoint(
                    Trace,
                    self.traces[self.g_min_index..self.traces.len],
                    g_star,
                    Trace.gStarMaxPredicate,
                ) + self.g_min_index;

                if (j == self.traces.len) {
                    return .{
                        .up1 = upper_index,
                        .up2 = upper_index - 1,
                        .low1 = j - 1,
                        .low2 = i,
                    };
                } else {
                    return .{
                        .up1 = upper_index,
                        .up2 = upper_index - 1,
                        .low1 = j - 1,
                        .low2 = j,
                    };
                }
            }

            return .{
                .up1 = upper_index,
                .up2 = upper_index -| 1,
                .low1 = i - 1,
                .low2 = i,
            };
        }

        /// Given a target radius and a list of angles on the image plane,
        /// returns a `CunninghamTransferFunction` with all values populated
        /// for each angle.
        /// The lifetime must then be managed by the caller, and deinit called.
        pub fn traceAllAlloc(
            allocator: std.mem.Allocator,
            tracer: Tracer,
            angles: []const T.T,
            target_radius: T.T,
            opts: TracerOptions,
        ) !Self {
            var traces: std.ArrayList(Trace) = try .initCapacity(allocator, angles.len);
            defer traces.deinit(allocator);

            var initial_guess: T.T = @max(opts.minimum_guess, opts.initial_guess);

            for (angles) |angle| {
                const trace = try tracer.single(angle, target_radius, initial_guess);
                traces.appendAssumeCapacity(trace);
                initial_guess = @max(opts.minimum_guess, trace.image_radius);
            }

            var bad_indices: [512]usize = undefined;
            var bad_count: usize = 0;
            for (traces.items, 0..) |trace, i| {
                if (isPathological(trace)) {
                    bad_indices[bad_count] = i;
                    bad_count += 1;
                }
            }

            // Remove all pathological traces
            traces.orderedRemoveMany(bad_indices[0..bad_count]);

            return .fromTraces(
                target_radius,
                try traces.toOwnedSlice(allocator),
            );
        }

        /// Options for controlling the heuristic-based transfer function
        /// calculator.
        pub const HeuristicOptions = struct {
            /// The maximum number of points to trace.
            max_points: usize = 200,
            /// The number of points to trace before starting to trace based
            /// off of the heuristic function.
            initial_trace: usize = 20,
            /// What the initial guess should be.
            initial_guess: T.T = 10.0,
            /// The minimum value for the initial guess.
            minimum_guess: T.T = 4.0,
            /// The error tolerance.
            error_tolerance: T.T = 1e-6,

            fn toTracerOptions(self: HeuristicOptions) TracerOptions {
                return .{
                    .initial_guess = self.initial_guess,
                    .minimum_guess = self.minimum_guess,
                };
            }
        };

        const HeuristicState = struct {
            trace: Trace,
            score: T.T,
            /// The neighbour to compare the score against
            neighbour: *HeuristicState,
        };

        /// Similar to `traceAll` but with a user-supplied heuristic function
        /// that compares neighbouring traces and should return a score. Those
        /// neighbours with the highest scores will be refined.
        ///
        /// Traces at most `traces.len` traces.
        pub fn traceAllHeuristicAlloc(
            allocator: std.mem.Allocator,
            tracer: Tracer,
            target_radius: T.T,
            comptime heuristic: fn (Trace, Trace) T.T,
            opts: HeuristicOptions,
        ) !Self {
            const states = try allocator.alloc(HeuristicState, opts.max_points);
            defer allocator.free(states);

            var count: usize = 0;

            var initial_guess = @max(opts.initial_guess, opts.minimum_guess);

            // Trace an initial values.
            const diff: T.T = std.math.pi * 2.0 / @as(T.T, @floatFromInt(opts.initial_trace));
            for (0..opts.initial_trace) |i| {
                const angle = @as(T.T, @floatFromInt(i)) * diff;
                const trace = try tracer.single(
                    angle + opts.error_tolerance,
                    target_radius,
                    initial_guess,
                );
                states[count] = .{
                    .trace = trace,
                    .score = if (count > 0) heuristic(states[count - 1].trace, trace) else 0,
                    .neighbour = if (count > 0) &states[count - 1] else undefined,
                };
                count += 1;
                initial_guess = @max(opts.minimum_guess, trace.image_radius);
            }
            // And make it cyclic.
            states[0].score = heuristic(states[count - 1].trace, states[0].trace);
            states[0].neighbour = &states[count - 1];

            while (count < opts.max_points) {
                // Find the current maximum score.
                var max_score: T.T = 0;
                var max_index: usize = 0;
                for (0..count) |i| {
                    if (states[i].score > max_score) {
                        max_score = states[i].score;
                        max_index = i;
                    }
                }

                const this = &states[max_index];
                const neighbour = states[max_index].neighbour;

                const angle_low = neighbour.trace.image_angle;
                const angle_high = this.trace.image_angle;
                // Find the new angle as the midpoint between the existing.
                const angle = angularMidpoint(T.T, angle_low, angle_high);

                const guess = @max(opts.minimum_guess, states[max_index].trace.image_radius);

                const trace = tracer.single(angle, target_radius, guess) catch |e| {
                    logger.err(
                        "Failed to calculate root for theta = {d:.6} ({d:.4}, initial guess: {d:.6})",
                        .{ angle, std.math.radiansToDegrees(angle), guess },
                    );
                    return e;
                };

                states[count] = .{
                    .trace = trace,
                    .score = heuristic(neighbour.trace, trace),
                    .neighbour = neighbour,
                };

                // Update the current maximum.
                states[max_index].neighbour = &states[count];
                states[max_index].score = heuristic(trace, this.trace);

                count += 1;
            }

            // Count how many are pathological.
            var bad_indices: [512]usize = undefined;
            var bad_count: usize = 0;
            for (states, 0..) |state, i| {
                if (isPathological(state.trace)) {
                    bad_indices[bad_count] = i;
                    bad_count += 1;
                }
            }

            // Unpack and remove all pathological traces.
            const traces = try allocator.alloc(Trace, count - bad_count);
            errdefer allocator.free(traces);
            var t_index: usize = 0;
            var s_index: usize = 0;
            var b_index: usize = 0;
            while (t_index < traces.len) {
                if (b_index < bad_count and s_index == bad_indices[b_index]) {
                    // Skip the bad one.
                    b_index += 1;
                    s_index += 1;
                    continue;
                }
                traces[t_index] = states[s_index].trace;
                t_index += 1;
                s_index += 1;
            }

            return .fromTraces(
                target_radius,
                traces,
            );
        }

        /// Uses an optimising algorithm to try to find better estimates of
        /// g_min and g_max, saving each point along the way. Uses a maximum of
        /// `N` evaluations to try to improve the result per extrema, for a
        /// total of `2N` points.
        pub fn optimiseExtremaAlloc(
            self: *Self,
            allocator: std.mem.Allocator,
            tracer: Tracer,
            N: usize,
        ) !void {
            var list: std.ArrayList(Trace) = .fromOwnedSlice(self.traces);
            defer list.deinit(allocator);

            // Preallocate
            try list.ensureUnusedCapacity(allocator, 2 * N);

            const Ctx = struct {
                self_: *Self,
                list_: *std.ArrayList(Trace),
                tracer_: Tracer,
                guess: T.T,

                fn scoreFunction(ctx: @This(), angle: T.T) T.T {
                    const trace = ctx.tracer_.single(angle, ctx.self_.target_radius, ctx.guess) catch {
                        unreachable;
                    };
                    if (!isPathological(trace)) {
                        ctx.list_.appendAssumeCapacity(trace);
                    }
                    return trace.g;
                }
            };

            const offset_angle = 0.3;
            const current_min = list.items[self.g_min_index];
            const current_max = list.items[self.g_max_index];

            // Optimise g_min
            var ctx: Ctx = .{
                .self_ = self,
                .list_ = &list,
                .tracer_ = tracer,
                .guess = current_min.image_radius,
            };

            _ = optimiseGoldenSection(T.T, ctx, Ctx.scoreFunction, .{
                .optimise = .min,
                .max_steps = N,
                .lower_bound = current_min.image_angle - offset_angle,
                .upper_bound = current_min.image_angle + offset_angle,
            });

            // Optimise g_max
            ctx.guess = current_max.image_radius;
            _ = optimiseGoldenSection(T.T, ctx, Ctx.scoreFunction, .{
                .optimise = .max,
                .max_steps = N,
                .lower_bound = current_max.image_angle - offset_angle,
                .upper_bound = current_max.image_angle + offset_angle,
            });

            self.* = .fromTraces(
                self.target_radius,
                try list.toOwnedSlice(allocator),
            );
        }

        /// Given a set of CTF, trace `M` refinemenets of the surrounding `N`
        /// points around g_min and g_max, for a total of `4 * N * M`
        /// additional traces.
        ///
        /// This works e.g. by taking the index `g_min_index - N, g_min_index -
        /// N + 1` and tracing another `M` evenly spaced angles between them,
        /// and doing so all the way up to `g_min + N`.
        ///
        /// After calling this function, the lifetime must be managed by the
        /// caller, and `deinit` must be called.
        pub fn refineExtremaAlloc(
            self: *Self,
            allocator: std.mem.Allocator,
            tracer: Tracer,
            N: usize,
            M: usize,
            opts: TracerOptions,
        ) !void {
            var list: std.ArrayList(Trace) = .fromOwnedSlice(self.traces);
            defer list.deinit(allocator);

            // Preallocate
            try list.ensureUnusedCapacity(allocator, 4 * N * M);

            for (&[_]usize{ self.g_min_index, self.g_max_index }) |start_index| {
                for (0..2 * N) |n| {
                    // TODO: allow wrapparounds
                    const i = start_index + n - N;

                    const angle_i = list.items[i].image_angle;
                    const delta_theta = angularDelta(
                        T.T,
                        angle_i,
                        list.items[i + 1].image_angle,
                    ) / @as(T.T, @floatFromInt(M + 1));

                    for (1..M) |m| {
                        const m_f: T.T = @floatFromInt(m);
                        const angle = @mod(angle_i + delta_theta * m_f, std.math.pi * 2.0);

                        const guess = @max(opts.minimum_guess, list.items[i].image_radius);

                        const trace = try tracer.single(angle, self.target_radius, guess);
                        list.appendAssumeCapacity(trace);
                    }
                }
            }

            var bad_indices: [512]usize = undefined;
            var bad_count: usize = 0;
            for (list.items, 0..) |trace, i| {
                if (isPathological(trace)) {
                    bad_indices[bad_count] = i;
                    bad_count += 1;
                }
            }

            // Remove all pathological traces
            list.orderedRemoveMany(bad_indices[0..bad_count]);

            self.* = .fromTraces(
                self.target_radius,
                try list.toOwnedSlice(allocator),
            );
        }

        pub const TracerOptions = struct {
            /// What the initial guess should be.
            initial_guess: T.T = 10.0,
            /// The minimum value for the initial guess.
            minimum_guess: T.T = 4.0,
        };

        /// The tracer state for calculating Cunningham's transfer function.
        pub const Tracer = struct {
            metric: KerrMetric(T),
            x_obs: FourVector(T),
            disc: accretion_discs.AccretionDisc(T),

            /// Initialise a tracer.
            pub fn init(
                metric: KerrMetric(T),
                x_obs: FourVector(T),
                disc: accretion_discs.AccretionDisc(T),
            ) Tracer {
                return .{
                    .metric = metric,
                    .x_obs = x_obs,
                    .disc = disc,
                };
            }

            /// Trace a single point on the image plane and calculate its CTF
            /// values. The angle should be specified 0 to 2π.
            fn single(
                self: Tracer,
                angle: T.T,
                target_radius: T.T,
                initial_guess: T.T,
            ) !Trace {
                // For root solving:
                const Dual1 = T.PushSlot();
                // For determining the Jacobian:
                const Dual2 = Dual1.PushSlot();

                // TODO: here use datum plane for thick discs.

                const res = try impactOffsetForRadius(
                    Dual1,
                    self.metric.adapt(Dual1),
                    self.x_obs.adapt(Dual1),
                    .promote(angle),
                    .promote(target_radius),
                    .{
                        .initial_guess = initial_guess,
                        .disc = self.disc.adapt(Dual1),
                    },
                );

                const p_2 = toImpactParameters(
                    Dual2,
                    .promote(res.r_offset),
                    .promote(res.theta_image_plane),
                );

                var v_medium = orbits.keplerianPlungingAlt(
                    T,
                    self.metric,
                    .promote(target_radius),
                );
                // Apply time-reversal trick, since the photon is actually
                // originating at the disc:
                v_medium.r = v_medium.r.neg();

                // retrace to calculate the Jacobian
                const dual_results = traceDual(
                    Dual2,
                    self.metric.adapt(Dual2),
                    self.x_obs.adapt(Dual2),
                    p_2.alpha,
                    p_2.beta,
                    self.disc.adapt(Dual2),
                    v_medium.adapt(Dual2),
                );

                return .{
                    .g = dual_results.g,
                    .alpha = p_2.alpha.x,
                    .beta = p_2.beta.x,
                    .jacobian = dual_results.jacobian,
                    .phi = dual_results.phi,
                    .image_radius = res.r_offset,
                    .delta_t = dual_results.t - self.x_obs.r.x,
                    .image_angle = angle,
                    .r_err = res.err,
                    .stokes_x = dual_results.stokes_x,
                    .stokes_y = dual_results.stokes_y,
                    .local_theta = dual_results.local_theta,
                    .local_phi = dual_results.local_phi,
                    // This has to be set later, since it depends knowing all
                    // of the `g` values.
                    .g_star = 0,
                    .f = 0,
                };
            }
        };

        const DualTrace = struct {
            jacobian: T.T,
            g: T.T,
            phi: T.T,
            t: T.T,
            stokes_x: T.T,
            stokes_y: T.T,
            local_theta: T.T,
            local_phi: T.T,
        };

        /// Trace a geodesic with the derivative slots for the alpha and beta
        /// impact parameters set.
        fn traceDual(
            comptime Dual: type,
            metric: KerrMetric(Dual),
            x_obs: FourVector(Dual),
            alpha: Dual,
            beta: Dual,
            disc: accretion_discs.AccretionDisc(Dual),
            v_medium: FourVector(Dual),
        ) DualTrace {
            if (Dual.N < 2) {
                @compileError("Requires at least two slots for traceDual");
            }

            const geod = NullGeodesic(Dual).fromImpactParameters(
                metric,
                x_obs,
                alpha.diff(0),
                beta.diff(1),
            );

            const result = geod.traceDisc(
                metric,
                disc,
                .{},
            );
            const eshift = redshift.redshiftFromResult(
                Dual,
                metric,
                geod,
                result,
                v_medium,
            );
            const total = result.totalAntiderivatives(metric, geod);

            // The Jacobian will be (d r g) / (d a b), but we need |(d a b) /
            // (d r g)|, so invert it:
            const jacobian = 1.0 / @abs(ad.jacobianDeterminant(Dual, result.r, eshift));

            // Sanity check:
            std.debug.assert(std.math.isFinite(jacobian));

            // TODO: reuse the tangent space calculation here
            const pc = result.polarisationConstant(metric, geod, v_medium);
            const xy = pc.polarisationXYAtInfinityAlt(
                metric,
                geod.x_init.th,
                alpha,
                beta,
            );

            const local_angles = result.localAngles(metric, geod, v_medium);

            return .{
                .g = eshift.x,
                .jacobian = jacobian,
                .phi = total.coordinateAzimuth(metric, geod).x,
                .t = total.coordinateTime(metric, geod).x,
                .stokes_x = xy.x.x,
                .stokes_y = xy.y.x,
                .local_phi = local_angles.phi.x,
                .local_theta = local_angles.theta.x,
            };
        }

        /// Returns true if the trace has a pathology, e.g. has nearly a zero
        /// Jacobian.
        ///
        /// TODO: Note this is a heuristic. It is a little bit tuned, and
        /// really a better way of discovering problematic transfer functions
        /// should be implemented.
        fn isPathological(trace: Trace) bool {
            return std.math.approxEqAbs(T.T, 0.0, 1.0 / trace.jacobian, 1e-6);
        }

        /// Calculate all of the remaining fields for the traces that depend on
        /// having a full set of trace calculations.
        fn fromTraces(target_radius: T.T, traces: []Trace) Self {
            std.sort.heap(Trace, traces, {}, Trace.ascendingAngle);

            const g_extrema = findEnergyShiftExtrema(traces);
            std.debug.assert(g_extrema.g_max_index != g_extrema.g_min_index);

            const delta_g = g_extrema.g_max - g_extrema.g_min;

            const common_prefactor = 1.0 / (std.math.pi * target_radius);

            // Now loop over every trace and calculate Cunningham's transfer function.
            for (traces) |*trace| {
                trace.g_star = gstarFrom(trace.g, g_extrema.g_min, g_extrema.g_max);
                // Convert the Jacobian to be (partial g^star)
                const jacobian = trace.jacobian * delta_g;

                // Then calculate Cunningham's actual transfer function.
                const sqrt_g_star = @sqrt(trace.g_star * (1 - trace.g_star));
                trace.f = common_prefactor * trace.g * sqrt_g_star * jacobian;
            }

            return .{
                .target_radius = target_radius,
                .traces = traces,
                .g_min = g_extrema.g_min,
                .g_min_index = g_extrema.g_min_index,
                .g_max = g_extrema.g_max,
                .g_max_index = g_extrema.g_max_index,
            };
        }

        fn findEnergyShiftExtrema(traces: []const Trace) struct {
            g_min: T.T,
            g_min_index: usize,
            g_max: T.T,
            g_max_index: usize,
        } {
            // It's impossible to double the energy, let alone make it 10 fold,
            // so this is a fine initial value.
            var g_min: T.T = 10.0;
            var g_max: T.T = 0;
            var g_min_index: usize = 0;
            var g_max_index: usize = 0;

            for (traces, 0..) |trace, i| {
                if (trace.g > g_max) {
                    g_max = trace.g;
                    g_max_index = i;
                }
                if (trace.g < g_min) {
                    g_min = trace.g;
                    g_min_index = i;
                }
            }

            return .{
                .g_max = g_max,
                .g_min = g_min,
                .g_min_index = g_min_index,
                .g_max_index = g_max_index,
            };
        }

        fn tracePhi(t: Trace) T.T {
            return t.phi;
        }

        /// Interpolate a `Trace` based on the azimuthal coordinate along the
        /// disc.
        pub fn interpolateDiscAzimuth(self: *const Self, phi: T.T) Trace {
            const info = interpolations.periodicInterpolate(Trace, self.traces, phi, tracePhi);
            var t1 = self.traces[info.left];
            if (info.left_mod) {
                t1.phi -= 2 * std.math.pi;
            }
            return t1.interpolateBetween(self.traces[info.right], info.weight);
        }

        /// Normalise all of the azimuthal angles between 0 and 2π.
        pub fn normaliseAngles(self: *Self) void {
            for (self.traces) |*trace| {
                trace.phi = @mod(std.math.pi + trace.phi, std.math.pi * 2.0);
            }
        }

        /// Used to determine the partition point in the radial axis of a set
        /// of CTFs. Use as:
        ///
        ///     std.sort.partitionPoint(CTF, ctfs, r.x, CTF.partitionRadial)
        ///
        /// to get first index with target radius greater than `r`.
        pub fn partitionRadial(r: T.T, self: Self) bool {
            return self.target_radius < r;
        }

        /// Used as the functional argument to `std.sort.` to sort a set of
        /// CTFs by target radius.
        pub fn sortByRadius(_: void, lhs: Self, rhs: Self) bool {
            return lhs.target_radius < rhs.target_radius;
        }
    };
}

test "transfer branch interpolation" {
    const Dual = DualNumber(f64, 0);
    const CTF = CunninghamTransferFunction(Dual);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    const x_obs: FourVector(Dual) = .{
        .t = .zero,
        .r = .promote(1e7),
        .th = .promote(std.math.degreesToRadians(60)),
        .ph = .zero,
    };

    const tracer = CTF.Tracer.init(metric, x_obs, .equatorial_plane);

    const angles: []const f64 = &.{
        0.001, 0.2, 0.4, 0.6, 0.8,
        1.0,   1.2, 1.4, 1.6, 1.8,
        2.0,   2.2, 2.4, 2.6, 2.8,
        3.0,   3.2, 3.4, 3.6, 3.8,
        4.0,   4.2, 4.4, 4.6, 4.8,
        5.0,   5.2, 5.4, 5.6, 5.8,
        6.0,   6.2,
    };

    var ctf = try CTF.traceAllAlloc(std.testing.allocator, tracer, angles, 5.0, .{});
    defer ctf.deinit(std.testing.allocator);

    // Keeping this around as it will likely need debugging again.
    // for (ctf.traces, 0..) |trace, i| {
    //     std.debug.print("{d: >3}: angle:{d:.4} g:{d:.4} gstar:{d:.4}\n", .{ i, trace.image_angle, trace.g, trace.g_star });
    //     if (i == ctf.g_max_index) {
    //         std.debug.print(" -- max index: {d}\n", .{i});
    //     }
    //     if (i == ctf.g_min_index) {
    //         std.debug.print(" -- min index: {d}\n", .{i});
    //     }
    // }

    {
        const interpolated = ctf.interpolationIndices(0.8604);
        try std.testing.expectEqualDeep(
            IndexSet{ .low1 = 5, .low2 = 6, .up1 = 9, .up2 = 8 },
            interpolated,
        );
    }
    {
        const interpolated = ctf.interpolationIndices(0.2826);
        try std.testing.expectEqualDeep(
            IndexSet{ .low1 = 31, .low2 = 0, .up1 = 17, .up2 = 16 },
            interpolated,
        );
    }
    {
        const interpolated = ctf.interpolationIndices(0.002047);
        try std.testing.expectEqualDeep(
            IndexSet{ .low1 = 25, .low2 = 26, .up1 = 25, .up2 = 24 },
            interpolated,
        );
    }

    const values = ctf.interpolateBranches(1.1);
    try std.testing.expectApproxEqAbs(0.1650734805778647, values.upper, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(0.2337572255338412, values.lower, TEST_TOLERANCE);
}

const TestCTFOptions = struct {
    spin: f64,
    incl: f64,
    r: f64,
};

fn calculateTestTransferFunction(
    allocator: std.mem.Allocator,
    opts: TestCTFOptions,
) !CunninghamTransferFunction(DualNumber(f64, 0)) {
    const Dual = DualNumber(f64, 0);
    const CTF = CunninghamTransferFunction(Dual);

    const metric: KerrMetric(Dual) = .init(.one, .promote(opts.spin));
    const x_obs: FourVector(Dual) = .{
        .t = .zero,
        .r = .promote(1e7),
        .th = .promote(std.math.degreesToRadians(opts.incl)),
        .ph = .zero,
    };

    const angles = try allocator.alloc(f64, 80);
    defer allocator.free(angles);
    for (angles, 0..) |*a, i| {
        const f = @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(angles.len + 1));
        const angle = f * std.math.pi * 2.0 + 1e-4;
        a.* = angle;
    }

    const tracer = CTF.Tracer.init(metric, x_obs, .equatorial_plane);
    var ctf = try CTF.traceAllAlloc(allocator, tracer, angles, opts.r, .{});
    try ctf.optimiseExtremaAlloc(allocator, tracer, 17);
    return ctf;
}

test "transfer function regression" {
    {
        var ctf = try calculateTestTransferFunction(
            std.testing.allocator,
            .{ .spin = 0.998, .incl = 80, .r = 1.2469706551751847 },
        );
        defer ctf.deinit(std.testing.allocator);
        // These have been calculated with Gradus.jl for comparison. kerrz
        // should at best match but preferably always outperform.
        try std.testing.expect(0.035052 >= ctf.g_min);
        try std.testing.expect(0.757320 <= ctf.g_max);
    }
}

test "disc interpolation" {
    {
        var ctf = try calculateTestTransferFunction(
            std.testing.allocator,
            .{ .spin = 0.998, .incl = 40, .r = 5.0 },
        );
        defer ctf.deinit(std.testing.allocator);

        const p = ctf.interpolateDiscAzimuth(0.5);
        try std.testing.expectApproxEqAbs(0.6052619993682032, p.g, TEST_TOLERANCE);
    }
}

/// Columns that can be selected for exporting.
pub const SelectedField = enum {
    /// The Cunningham transfer function.
    f,
    /// The disc-to-observer coordinate time.
    delta_t,
    /// The alpha impact parameter.
    alpha,
    /// The beta impact parameter.
    beta,
    /// The local emission angle.
    local_theta,
    /// Stokes X = Q/I.
    stokes_x,
    /// Stokes Y = U/I.
    stokes_y,
    /// The azimuthal coordinate on the disc.
    phi,

    /// Get a string that describes the units of this field.
    pub fn getUnits(self: SelectedField) ?[]const u8 {
        return switch (self) {
            .alpha => "rg",
            .beta => "rg",
            .local_theta => "radians",
            .delta_t => "tg",
            else => null,
        };
    }

    /// Get a short string that described what this field is.
    pub fn getShortDescriptor(self: SelectedField) ?[]const u8 {
        return switch (self) {
            .f => "Cunningham's transfer function",
            .delta_t => "Disc-to-observer time, minus r_obs",
            .alpha => "The alpha impact parameter",
            .beta => "The alpha impact parameter",
            .local_theta => "Angle between photon and disc",
            else => null,
        };
    }

    /// Work out the index of the field in `Trace` that this selected field
    /// corresponds to.
    pub fn toIndex(self: SelectedField, comptime T: type) usize {
        switch (self) {
            inline else => |name| {
                return std.meta.fieldIndex(
                    CunninghamTransferFunction(T).Trace,
                    @tagName(name),
                ).?;
            },
        }
    }

    fn getField(self: SelectedField, comptime T: type, trace: *const CunninghamTransferFunction(T).Trace) T.T {
        return self.getFieldPtr(T, @constCast(trace)).*;
    }

    fn getFieldPtr(self: SelectedField, comptime T: type, trace: *CunninghamTransferFunction(T).Trace) *T.T {
        return switch (self) {
            inline else => |name| &@field(trace, @tagName(name)),
        };
    }

    fn fromName(name: []const u8) ?SelectedField {
        const trimmed_name = name[0 .. name.len - ("_upper").len];
        // Handle exceptional names separately.
        if (std.mem.eql(u8, trimmed_name, "theta")) {
            return .local_theta;
        }
        return std.meta.stringToEnum(SelectedField, trimmed_name);
    }

    fn nameUpper(self: SelectedField) []const u8 {
        return switch (self) {
            .local_theta => "theta_upper",
            inline else => |name| @tagName(name) ++ "_upper",
        };
    }

    fn nameLower(self: SelectedField) []const u8 {
        return switch (self) {
            .local_theta => "theta_lower",
            inline else => |name| @tagName(name) ++ "_lower",
        };
    }
};

/// A table of transfer functions that provides some utility functions for
/// seralising and integrating.
///
/// This is the principle way that a caller will use the transfer function
/// calculators.
pub fn CunninghamTransferFunctionTable(comptime T: type) type {
    return struct {
        pub const Context = TransferFunctionContext(T);
        pub const TransferFunction = CunninghamTransferFunction(T);
        pub const CTF = CunninghamTransferFunction(T);
        pub const Options = Context.Options;

        const Self = @This();

        tracer: CTF.Tracer,
        transfer_functions: std.ArrayList(CTF),

        /// This records which fields of the Trace are active / known. If a
        /// table is read in from file, not all of the fields may be present,
        /// and this can be used to check if the necessary information for a
        /// calculation is present.
        active_fields: std.bit_set.StaticBitSet(@typeInfo(CTF.Trace).@"struct".fields.len) = .initFull(),

        /// Initialise a Cunningham transfer function table. Caller must call
        /// `deinit()`.
        pub fn init(
            metric: KerrMetric(T),
            x_obs: FourVector(T),
            disc: accretion_discs.AccretionDisc(T),
        ) Self {
            const tracer: CTF.Tracer = .init(metric, x_obs, disc);
            return .{
                .tracer = tracer,
                .transfer_functions = .empty,
            };
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            for (self.transfer_functions.items) |*tf| {
                tf.deinit(allocator);
            }
            self.transfer_functions.deinit(allocator);
            self.* = undefined;
        }

        fn isSorted(self: *const Self) bool {
            var r_prev = self.transfer_functions.items[0].target_radius;
            for (self.transfer_functions.items[1..]) |next| {
                if (next.target_radius < r_prev) {
                    return false;
                }
                r_prev = next.target_radius;
            }
            return true;
        }

        /// Interpolate a `Trace` at a set of coordinates on the accretion
        /// disc.
        pub fn interpolateDiscCoordinates(
            self: *const Self,
            r: T.T,
            phi: T.T,
        ) CTF.Trace {
            std.debug.assert(phi >= 0 and phi < std.math.pi * 2);

            const ctfs = self.transfer_functions.items;

            const index = std.sort.partitionPoint(
                CTF,
                ctfs,
                r,
                CTF.partitionRadial,
            );

            // Check the edge cases:
            if (index == ctfs.len) {
                const tf = ctfs[index - 1];
                return tf.interpolateDiscAzimuth(phi);
            }
            if (index == 0) {
                const tf = ctfs[index];
                return tf.interpolateDiscAzimuth(phi);
            }

            const tf1 = ctfs[index - 1];
            const tf2 = ctfs[index];

            const p1 = tf1.interpolateDiscAzimuth(phi);
            const p2 = tf2.interpolateDiscAzimuth(phi);

            const w = interpolations.lerpWeight(
                T.T,
                r,
                tf1.target_radius,
                tf2.target_radius,
            );

            return p1.interpolateBetween(p2, w);
        }

        /// Calculate Cunningham's tranfer function for a particular target
        /// radius.
        ///
        /// In appending the CTF, the table takes ownership of the table, and
        /// will call `deinit()` when cleaning up.
        pub fn append(
            self: *Self,
            allocator: std.mem.Allocator,
            transfer_function: CTF,
        ) !void {
            try self.transfer_functions.append(allocator, transfer_function);
        }

        /// Calculate a transfer function for a particular radius.
        ///
        /// See also `calculateRadiusAndAppend`.
        pub fn calculateRadius(
            self: *const Self,
            allocator: std.mem.Allocator,
            target_radius: T.T,
            opts: Context.Options,
        ) !CTF {
            var ctf = try Context.calculateRadius(
                allocator,
                self.tracer,
                target_radius,
                opts,
            );
            errdefer ctf.deinit(allocator);

            if (opts.refine.N > 0) {
                const ref = opts.refine;
                try ctf.refineExtremaAlloc(
                    allocator,
                    self.tracer,
                    ref.N,
                    ref.M,
                    .{
                        .minimum_guess = opts.minimum_guess,
                        .initial_guess = opts.initial_guess orelse target_radius,
                    },
                );
            }

            return ctf;
        }

        /// The same as `calculateRadius` but also appends it to the set of
        /// transfer functions owned by this table. Returns a pointer to the
        /// new element which may be safely discarded if not needed.
        ///
        /// This function checks that the resulting table is still sorted in
        /// radius, and, if it is not, sorts the table.
        pub fn calculateRadiusAndAppend(
            self: *Self,
            allocator: std.mem.Allocator,
            target_radius: T.T,
            opts: Context.Options,
        ) !*CTF {
            var ctf = try self.calculateRadius(allocator, target_radius, opts);
            errdefer ctf.deinit(allocator);

            var ptr = try self.transfer_functions.addOne(allocator);
            ptr.* = ctf;

            if (!self.isSorted()) {
                std.sort.heap(CTF, self.transfer_functions.items, {}, CTF.sortByRadius);
                // Find the transfer function again after sorting:
                // TODO: some assertion should be added to avoid having two
                // transfer functions with exactly the same radius.
                for (self.transfer_functions.items) |*tf| {
                    if (tf.target_radius == target_radius) {
                        ptr = tf;
                    }
                }
            }

            return ptr;
        }

        /// Return a utility instance for integrating the transfer functions
        /// over a particular energy grid.
        ///
        /// Caller must call `deinit` on the TableIntegrator. It will only
        /// clear integration cache resources, and not the table itself.
        pub fn integrator(
            self: Self,
            allocator: std.mem.Allocator,
            g_grid: []const T.T,
            g_fine_len: usize,
        ) !TableIntegrator(T) {
            std.debug.assert(self.isSorted());
            return try .init(allocator, self, g_grid, g_fine_len);
        }

        /// Serialise the transfer function as a FITS binary table HDU.
        ///
        /// Note, this does not include every field in the `Trace`, but rather
        /// a common subset. Use `toFITSSelected` to select desired fields.
        pub fn toFITS(self: *const Self, allocator: std.mem.Allocator) !zfits.Hdu {
            return self.toFITSSelected(allocator, &.{ .f, .delta_t, .local_theta });
        }

        /// Serialise the selected columns of the transfer function table to a
        /// FITS HDU.
        pub fn toFITSSelected(
            self: *const Self,
            allocator: std.mem.Allocator,
            selected_fields: []const SelectedField,
        ) !zfits.Hdu {
            const N = 40;
            var hdu = zfits.Hdu.init(allocator, .{ .binary_table = .empty });
            errdefer hdu.deinit();

            try hdu.setName(
                "CTF_TABLE",
                "This is a Cunningham transfer function table",
            );
            try utils.addKerrzFITSInfo(&hdu);

            // Add some useful bits of information to the header.
            try hdu.appendHeaderRecord("GSTAR_H", .{
                .comment = "The offset either side of g_star limits",
                .value = .{ .float = EDGE_TOLERANCE },
            });

            // Metric information:
            try self.tracer.metric.addToHdu(&hdu);
            // Observer information:
            try utils.addObserverInformation(&hdu, self.tracer.x_obs);

            // Setup the columns
            try hdu.data.binary_table.appendColumn(.{
                .label = "radius",
                .comment = "The radius on the accretion disc",
                .units = "rg",
                .units_comment = "In units of rg = GM/c^2",
            });
            try hdu.data.binary_table.appendColumn(.{
                .label = "gmin",
                .comment = "The smallest energyshift from this radius",
            });
            try hdu.data.binary_table.appendColumn(.{
                .label = "gmax",
                .comment = "The largest energyshift from this radius",
            });

            for (selected_fields) |field| {
                try hdu.data.binary_table.appendColumn(
                    .{
                        .label = field.nameLower(),
                        .col_type = .{ .repeat = N },
                        .units = field.getUnits(),
                        .comment = field.getShortDescriptor(),
                    },
                );
                try hdu.data.binary_table.appendColumn(
                    .{
                        .label = field.nameUpper(),
                        .col_type = .{ .repeat = N },
                        .units = field.getUnits(),
                        .comment = field.getShortDescriptor(),
                    },
                );
            }

            for (self.transfer_functions.items) |ctf| {
                const row = try hdu.data.binary_table.addRow();
                row.cols[0].one.float_32 = @floatCast(ctf.target_radius);
                row.cols[1].one.float_32 = @floatCast(ctf.g_min);
                row.cols[2].one.float_32 = @floatCast(ctf.g_max);

                var g_itt = gstarIterator(f32, N);
                for (0..N) |i| {
                    const g_star = g_itt.next().?;
                    const trace = ctf.interpolateTraceGstar(g_star).?;
                    for (selected_fields, 0..) |field, index| {
                        const lower_index = 3 + (index * 2);
                        const upper_index = lower_index + 1;
                        row.cols[lower_index].many[i].float_32 = @floatCast(
                            field.getField(T, &trace.lower),
                        );
                        row.cols[upper_index].many[i].float_32 = @floatCast(
                            field.getField(T, &trace.upper),
                        );
                    }
                }
            }

            return hdu;
        }

        /// Used to check if a particular field of the `Trace` is present in
        /// this transfer function table.
        pub fn hasField(self: *const Self, field: SelectedField) bool {
            return self.active_fields.isSet(field.toIndex(T));
        }

        /// Initialise a transfer function table from a FITS HDU.
        pub fn fromSingleFITS(
            allocator: std.mem.Allocator,
            hdu: *zfits.Hdu,
        ) !Self {
            // Read the metric and observer information from the HDU records:
            const mass = hdu.getRecord("MASS").?.value.float;
            const spin = hdu.getRecord("SPIN").?.value.float;
            const r_obs = hdu.getRecord("R_OBS").?.value.float;
            const r_incl = hdu.getRecord("INCL_OBS").?.value.float;
            const gstar_h = hdu.getRecord("GSTAR_H").?.value.float;

            const metric: KerrMetric(T) = .init(
                .promote(@floatCast(mass)),
                .promote(@floatCast(spin)),
            );
            const x_obs: FourVector(T) = .{
                .t = .zero,
                .r = .promote(@floatCast(r_obs)),
                .th = .promote(@floatCast(std.math.degreesToRadians(r_incl))),
                .ph = .zero,
            };

            // TODO: read the disc information from the FITS file.
            var self: Self = .init(metric, x_obs, .equatorial_plane);
            errdefer self.deinit(allocator);

            const data = &hdu.data.binary_table;

            // Set all fields inactive.
            self.active_fields.mask = 0;

            // Now determine which fields are in the table:
            for (data.column_headers[3..]) |header| {
                const field = SelectedField.fromName(header.label).?;
                self.active_fields.set(field.toIndex(T));
            }

            const all_transfer_functions = try allocator.alloc(CTF, data.num_rows);
            errdefer allocator.free(all_transfer_functions);

            for (all_transfer_functions) |*ctf| {
                // Initialise an empty ctf
                ctf.* = .{
                    .traces = &.{},
                    .g_max = 0,
                    .g_max_index = 0,
                    .g_min = 0,
                    .g_min_index = 0,
                    .target_radius = 0,
                };
            }

            errdefer for (all_transfer_functions) |*tf| {
                if (tf.traces.len > 0) {
                    allocator.free(tf.traces);
                }
            };

            for (0.., all_transfer_functions) |row_index, *ctf| {
                // TODO: reuse a container since all the columns will be the
                // same shape.
                const row = try data.getOrParseRow(allocator, row_index);

                // Allocate the traces. Twice as many as the number of items
                // are needed, as it will be the lower and upper branches
                // concatenated.
                const N = row.cols[3].many.len;
                const traces = try allocator.alloc(CTF.Trace, 2 * N);
                errdefer allocator.free(traces);

                for (traces) |*trace| trace.* = .zero;

                // Do the fixed known columns:
                ctf.target_radius = @floatCast(row.cols[0].one.float_32);
                ctf.g_min = @floatCast(row.cols[1].one.float_32);
                ctf.g_max = @floatCast(row.cols[2].one.float_32);
                ctf.g_max_index = N;
                ctf.g_min_index = 2 * N - 1;

                // Set the g_star fields:
                var g_itt = gstarIteratorTol(f32, N, gstar_h);
                for (0..N) |index| {
                    const t1 = &traces[index];
                    const t2 = &traces[index + N];
                    t1.g_star = g_itt.next().?;
                    t1.g = CTF.gFrom(t1.g_star, ctf.g_min, ctf.g_max);
                    t2.g_star = t1.g_star;
                    t2.g = t1.g;
                }

                // Read in the rest of the data:
                for (data.column_headers[3..], 3..) |header, col_index| {
                    const field = SelectedField.fromName(header.label).?;
                    const offset = if (std.mem.endsWith(u8, header.label, "_lower")) N else 0;
                    // Read in the corresponding value into each trace:
                    for (0..N) |trace_index| {
                        const i = trace_index + offset;
                        field.getFieldPtr(
                            T,
                            &traces[i],
                        ).* = @floatCast(row.cols[col_index].many[trace_index].float_32);
                    }
                }

                // Reverse the lower branch indices:
                std.mem.reverse(CTF.Trace, traces[N..]);

                ctf.traces = traces;
            }

            self.transfer_functions = .fromOwnedSlice(all_transfer_functions);

            return self;
        }

        /// Create a duplicate of all the memory this table uses.
        fn dupe(self: *const Self, allocator: std.mem.Allocator) !Self {
            // Duplicate the first table for the interpolation chace.
            const duped_tables = try allocator.dupe(CTF, self.transfer_functions.items);
            errdefer allocator.free(duped_tables);

            // Copy each CTF
            var dupe_count: usize = 0;
            errdefer for (0..dupe_count) |i| {
                duped_tables[i].deinit(allocator);
            };

            for (duped_tables, self.transfer_functions.items) |*d, orig| {
                d.* = orig;
                // And duplicate it's table
                d.traces = try allocator.dupe(CTF.Trace, orig.traces);
                dupe_count += 1;
            }

            return .{
                .active_fields = self.active_fields,
                .tracer = self.tracer,
                .transfer_functions = .fromOwnedSlice(duped_tables),
            };
        }
    };
}

test "serialising and de-serialising" {
    const Dual = DualNumber(f64, 0);
    const metric: KerrMetric(Dual) = .init(.promote(1.0), .promote(0.998));
    const x_obs: FourVector(Dual) = .{
        .t = .zero,
        .r = .promote(1e7),
        .th = .promote(std.math.degreesToRadians(30.0)),
        .ph = .zero,
    };

    const allocator = std.testing.allocator;

    var table = CunninghamTransferFunctionTable(Dual).init(
        metric,
        x_obs,
        .equatorial_plane,
    );
    defer table.deinit(allocator);
    // Add a single table:
    _ = try table.calculateRadiusAndAppend(allocator, 10.0, .{});

    var hdu = try table.toFITS(allocator);
    defer hdu.deinit();

    var parsed_table = try CunninghamTransferFunctionTable(Dual).fromSingleFITS(
        allocator,
        &hdu,
    );
    defer parsed_table.deinit(allocator);

    // Were the g_star, g, and f values correctly read back in?
    const parsed1 = parsed_table.transfer_functions.items[0];
    const ctf1 = table.transfer_functions.items[0];

    try std.testing.expect(parsed_table.hasField(.f));
    try std.testing.expect(!parsed_table.hasField(.stokes_x));

    try std.testing.expect(table.hasField(.f));
    try std.testing.expect(table.hasField(.stokes_x));

    try std.testing.expectApproxEqAbs(
        ctf1.g_max,
        parsed1.g_max,
        TEST_TOLERANCE,
    );

    try std.testing.expectApproxEqAbs(
        ctf1.g_min,
        parsed1.g_min,
        TEST_TOLERANCE,
    );

    // Check a sample of the traces
    for (0..10) |trace_index| {
        const parsed_trace = parsed1.traces[trace_index];
        const true_trace = ctf1.interpolateTraceGstar(parsed_trace.g_star).?;

        try std.testing.expectApproxEqAbs(
            true_trace.upper.g,
            parsed_trace.g,
            TEST_TOLERANCE,
        );

        try std.testing.expectApproxEqAbs(
            true_trace.upper.local_theta,
            parsed_trace.local_theta,
            TEST_TOLERANCE,
        );

        try std.testing.expectApproxEqAbs(
            true_trace.upper.f,
            parsed_trace.f,
            TEST_TOLERANCE,
        );
    }
}

test "transfer function tables" {
    const Dual = DualNumber(f64, 0);
    const metric: KerrMetric(Dual) = .init(.promote(1.0), .promote(0.998));
    const x_obs: FourVector(Dual) = .{
        .t = .zero,
        .r = .promote(1e7),
        .th = .promote(std.math.degreesToRadians(30.0)),
        .ph = .zero,
    };

    const allocator = std.testing.allocator;

    var table = CunninghamTransferFunctionTable(Dual).init(
        metric,
        x_obs,
        .equatorial_plane,
    );
    defer table.deinit(allocator);

    // Calculate two different transfer functions
    (try table.calculateRadiusAndAppend(allocator, 10.0, .{})).normaliseAngles();
    (try table.calculateRadiusAndAppend(allocator, 12.0, .{})).normaliseAngles();

    // Now try the interpolation
    {
        const trace = table.interpolateDiscCoordinates(11.0, 1.0);
        try std.testing.expectApproxEqAbs(
            1.0,
            trace.phi,
            TEST_TOLERANCE,
        );
    }
    {
        const trace = table.interpolateDiscCoordinates(11.0, 1e-3);
        try std.testing.expectApproxEqAbs(
            1e-3,
            trace.phi,
            TEST_TOLERANCE,
        );
    }
}

/// An integrator for integrating lineprofiles or impulse response functions
/// from transfer functions.
pub fn TableIntegrator(comptime T: type) type {
    return struct {
        const Table = CunninghamTransferFunctionTable(T);
        const Branch = Table.CTF.InterpolatedBranches;
        const Self = @This();

        pub const IntegrationOptions = struct {
            num_r_steps: usize = 3000,
            r_step_grid: Grid(T.T) = .log10,
            r_min: T.T,
            r_max: T.T,
            emissivity: EmissivityProfile(T) = .{ .powerlaw = .{} },
            zero_time: T.T = 0,
            /// Normalise the lineprofile so that the sum under the profile is
            /// unity. If this is `false`, the line profile has physical units
            /// of flux.
            normalise: bool = true,
        };

        table: Table,
        g_grid: []const T.T,

        /// These are caches used to store the g_star interpolated transfer
        /// function branches.
        branches_1: []Branch, // TODO: Struct-of-Arrays here?
        branches_2: []Branch, // TODO: Struct-of-Arrays here?

        // An allocated temporary flux grid that is calculated per radius. This
        // is needed for convolving the lineprofiles with some underlying
        // reflection spectrum.
        flux_cache: []T.T,

        pub fn init(
            allocator: std.mem.Allocator,
            table: Table,
            g_grid: []const T.T,
            g_fine_len: usize,
        ) !Self {
            // TODO: this is only really needed for the spectral convolution.
            // Does it always need to be allocated?
            const flux_cache = try allocator.alloc(T.T, g_grid.len);
            errdefer allocator.free(flux_cache);

            const branches_1 = try allocator.alloc(Branch, g_fine_len);
            errdefer allocator.free(branches_1);

            const branches_2 = try allocator.alloc(Branch, g_fine_len);
            errdefer allocator.free(branches_2);

            return .{
                .table = table,
                .g_grid = g_grid,
                .branches_1 = branches_1,
                .branches_2 = branches_2,
                .flux_cache = flux_cache,
            };
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.flux_cache);
            allocator.free(self.branches_1);
            allocator.free(self.branches_2);
            self.* = undefined;
        }

        inline fn adjustIndex(grid: []const T.T, index: usize, target: T.T) usize {
            var i: usize = index;
            while (i != 0 and grid[i] >= target) {
                i -= 1;
            }
            while (i != grid.len - 1 and grid[i] < target) {
                i += 1;
            }
            return i;
        }

        const LineprofileContext = struct {
            parent: *Self,
            flux: []T.T,
            g_index: usize = 0,

            inline fn addFluxToEnergy(
                self: *LineprofileContext,
                g: T.T,
                flux_weight: T.T,
                branch_1: Branch,
                branch_2: Branch,
                r_weight: T.T,
                t: T.T,
            ) void {
                _ = t;
                const lower = interpolations.lerpValue(
                    T.T,
                    r_weight,
                    branch_1.lower,
                    branch_2.lower,
                );
                const upper = interpolations.lerpValue(
                    T.T,
                    r_weight,
                    branch_1.upper,
                    branch_2.upper,
                );

                // Find the bin to put this value in. Assuming the
                // indexes never change very much and are approximately
                // only a bin or two off.
                self.g_index = adjustIndex(self.parent.g_grid, self.g_index, g);

                const flux_element = flux_weight * (lower + upper);
                self.flux[self.g_index] += flux_element;
            }
        };

        const ImpulseResponseContext = struct {
            g_grid: []const T.T,
            t_grid: []const T.T,
            flux: Matrix(T.T),
            g_index: usize = 0,
            t_index_lower: usize = 0,
            t_index_upper: usize = 0,
            zero_time: T.T,

            inline fn addFluxToResponse(
                self: *ImpulseResponseContext,
                g: T.T,
                flux_weight: T.T,
                branch_1: Branch,
                branch_2: Branch,
                r_weight: T.T,
                t: T.T,
            ) void {
                const lower = interpolations.lerpValue(
                    T.T,
                    r_weight,
                    branch_1.lower,
                    branch_2.lower,
                );
                const upper = interpolations.lerpValue(
                    T.T,
                    r_weight,
                    branch_1.upper,
                    branch_2.upper,
                );

                const delta_t_lower = interpolations.lerpValue(
                    T.T,
                    r_weight,
                    branch_1.delta_t_lower,
                    branch_2.delta_t_lower,
                );
                const delta_t_upper = interpolations.lerpValue(
                    T.T,
                    r_weight,
                    branch_1.delta_t_upper,
                    branch_2.delta_t_upper,
                );

                self.g_index = adjustIndex(self.g_grid, self.g_index, g);
                self.t_index_lower = adjustIndex(
                    self.t_grid,
                    self.t_index_lower,
                    t + delta_t_lower - self.zero_time,
                );
                self.t_index_upper = adjustIndex(
                    self.t_grid,
                    self.t_index_upper,
                    t + delta_t_upper - self.zero_time,
                );

                // Sum the contributions into the impulse response.
                self.flux.getPtr(
                    self.g_index,
                    self.t_index_lower,
                ).* += lower * flux_weight;

                self.flux.getPtr(
                    self.g_index,
                    self.t_index_upper,
                ).* += upper * flux_weight;
            }
        };

        /// Integrate the transfer function table with some control provided by
        /// `IntegrationOptions` into a lineprofile.
        ///
        /// The output is written into `flux`.
        pub fn lineprofile(
            self: *Self,
            flux: []T.T,
            opts: IntegrationOptions,
        ) void {
            std.debug.assert(flux.len == self.g_grid.len);

            // Zero the flux output.
            @memset(flux, 0);
            var ctx: LineprofileContext = .{
                .parent = self,
                .flux = flux,
            };
            self.integrateImpl2(&ctx, opts, LineprofileContext.addFluxToEnergy);

            const dg = self.g_grid[1] - self.g_grid[0];

            // Normalise the output
            if (opts.normalise) {
                var total_flux: T.T = 0;
                for (flux) |f| total_flux += f;
                if (total_flux > 0) {
                    for (flux) |*f| f.* /= (total_flux * dg);
                }
            } else {
                for (flux) |*f| f.* /= dg;
            }
        }

        /// Integrate the (2D) impulse response function.
        ///
        /// The output is written into the given matrix. The time-axis must be
        /// supplied by the caller.
        pub fn impulseResponse(
            self: *Self,
            t_grid: []const T.T,
            flux: Matrix(T.T),
            opts: IntegrationOptions,
        ) void {
            std.debug.assert(flux.n_rows == self.g_grid.len);
            std.debug.assert(flux.n_cols == t_grid.len);

            // Zero the flux output.
            @memset(flux.values, 0);

            var ctx: ImpulseResponseContext = .{
                .g_grid = self.g_grid,
                .t_grid = t_grid,
                .flux = flux,
                .zero_time = opts.zero_time,
            };
            self.integrateImpl2(&ctx, opts, ImpulseResponseContext.addFluxToResponse);

            // TODO: normalise the impulse response.
        }

        /// Convolve a disc spectrum. The flux grid must be setup with the
        /// correct dimensions, i.e. one less than the length of the disc
        /// spectrum energy grid.
        pub fn convolve(
            self: Self,
            flux: []T.T,
            disc_spectrum: spectra.DiscSpectrum,
            opts: IntegrationOptions,
        ) void {
            const n_flux = disc_spectrum.energy_grid.len - 1;
            std.debug.assert(flux.len == n_flux);
            self.integrateImpl(flux, opts, disc_spectrum, .with_spectrum);
        }

        /// The prototype for functions that are used to sum together
        /// infinitesimal flux elements.
        fn SumFluxPrototype(comptime Ctx: type) type {
            return fn (
                Ctx,
                g: T.T,
                flux_weight: T.T,
                branch_1: Branch,
                branch_2: Branch,
                r_weight: T.T,
                t: T.T, // The corona-to-disc time.
            ) callconv(.@"inline") void;
        }

        fn integrateImpl2(
            self: Self,
            ctx: anytype,
            opts: IntegrationOptions,
            comptime sumFlux: SumFluxPrototype(@TypeOf(ctx)),
        ) void {
            const tfs = self.table.transfer_functions.items;

            const n_radii = @divFloor(opts.num_r_steps, tfs.len);

            outer: for (1..tfs.len) |tf_index| {
                const tf1 = tfs[tf_index - 1];
                const tf2 = tfs[tf_index];

                // Skip those below the minimum radius
                if (opts.r_min > tf2.target_radius) continue;
                // Do not go past the outermost radius
                if (opts.r_max < tf1.target_radius) break;

                var g_star_grid = gstarIterator(T.T, self.branches_1.len);

                var branch_index: usize = 0;
                while (g_star_grid.next()) |g_star| {
                    // Cache all of the interpolated branches for these
                    // transfer functions on the fine g_star grid.
                    self.branches_1[branch_index] = tf1.interpolateBranchesGstar(g_star);
                    self.branches_2[branch_index] = tf2.interpolateBranchesGstar(g_star);
                    branch_index += 1;
                }

                // Construct a radial grid.
                var r_itt = opts.r_step_grid.iterator(
                    @max(opts.r_min, tf1.target_radius),
                    @min(opts.r_max, tf2.target_radius),
                    n_radii,
                );

                var r = r_itt.next() orelse continue :outer;
                while (r_itt.next()) |r_next| {
                    const delta_r = r_next - r;

                    // TODO: precompute this before this loop?
                    const em_vals = opts.emissivity.radialValues(.promote(r));

                    const r_weight = interpolations.lerpWeight(
                        T.T,
                        r,
                        tf1.target_radius,
                        tf2.target_radius,
                    );

                    // Interpolate the extremal redshift values.
                    const g_min = interpolations.lerpValue(
                        T.T,
                        r_weight,
                        tf1.g_min,
                        tf2.g_min,
                    );
                    const g_max = interpolations.lerpValue(
                        T.T,
                        r_weight,
                        tf1.g_max,
                        tf2.g_max,
                    );
                    const weight = r * em_vals.em.x * std.math.pi * delta_r / (g_max - g_min);

                    // Reset the g_star grid:
                    g_star_grid = gstarIterator(T.T, self.branches_1.len);
                    branch_index = 0;

                    var g_star_prev = g_star_grid.next().?;

                    while (g_star_grid.next()) |g_star| {
                        const g_star_delta = @abs(g_star - g_star_prev);

                        const branches_1 = self.branches_1[branch_index];
                        const branches_2 = self.branches_2[branch_index];
                        branch_index += 1;

                        const inv_denom = 1 / @sqrt(g_star * (1 - g_star));

                        const g = Table.CTF.gFrom(g_star, g_min, g_max);
                        const delta_g = Table.CTF.gFrom(g_star + g_star_delta, g_min, g_max) - g;

                        const flux_weight = g * g * weight * inv_denom * delta_g;
                        sumFlux(
                            ctx,
                            g,
                            flux_weight,
                            branches_1,
                            branches_2,
                            r_weight,
                            em_vals.t.x,
                        );

                        g_star_prev = g_star;
                    }

                    // Integrate the edges of the g_star grid, i.e. [0, eps) and (1-eps, 1].
                    // const g_low = Table.CTF.gFrom(eps, g_min, g_max);
                    // const g_high = Table.CTF.gFrom(1 - eps, g_min, g_max);

                    // self.integrateEdge(0, g_min, g_low, r_weight, weight);
                    // self.integrateEdge(self.branches_1.len - 1, g_high, g_max, r_weight, weight);

                    r = r_next;
                }
            }
        }

        fn integrateEdge(
            self: Self,
            branch_index: usize,
            g_low: T.T,
            g_high: T.T,
            r_weight: T.T,
            weight: T.T,
        ) void {
            std.debug.assert(g_low <= g_high);
            const branches_1 = self.branches_1[branch_index];
            const branches_2 = self.branches_2[branch_index];

            const lower = interpolations.lerpValue(
                T.T,
                r_weight,
                branches_1.lower,
                branches_2.lower,
            );
            const upper = interpolations.lerpValue(
                T.T,
                r_weight,
                branches_1.upper,
                branches_2.upper,
            );

            const edge_weight = weight * (@sqrt(g_high) - @sqrt(g_low));
            const edge_f = lower + upper;

            const g = (g_low + g_high) / 2;

            // Find the g bin index:
            var g_index: usize = 0;
            while (g_index != 0 and self.g_grid[g_index] > g) g_index -= 1;
            while (g_index != self.g_grid.len - 1 and self.g_grid[g_index] < g) g_index += 1;

            self.flux_cache[g_index] += edge_f * g * g * edge_weight;
        }

        fn integrateImpl(
            self: Self,
            flux: []T.T,
            opts: IntegrationOptions,
            disc_spectrum: ?spectra.DiscSpectrum,
            comptime spectrum: enum { with_spectrum, no_spectrum },
        ) void {
            const tfs = self.table.transfer_functions.items;

            // Zero the output flux and the flux buffer.
            @memset(flux, 0);
            @memset(self.flux_cache, 0);

            const delta_g = (self.g_grid[1] - self.g_grid[0]) / @as(
                T.T,
                @floatFromInt(opts.g_refine + 1),
            );

            // Construct a radial grid.
            var r_itt = opts.r_step_grid.iterator(
                T.T,
                tfs[0].target_radius,
                tfs[tfs.len - 1].target_radius,
                opts.num_r_steps,
            );

            var tf_index: usize = 1;
            var r = r_itt.next().?;

            // The current spectrum that is being considered. Only applies for
            // convolutions.
            var spectrum_index: usize = 0;

            radial: while (r_itt.next()) |r_next| {
                switch (spectrum) {
                    .no_spectrum => {},
                    .with_spectrum => {
                        // If we are convolving, skip those that do not contribute.
                        if (r < disc_spectrum.?.radii[spectrum_index]) continue;
                    },
                }

                const delta_r = r_next - r;

                // Find the corresponding transfer function pair. Loop in reverse over
                // the transfer functions to find the first that has a target radius
                // equal to or higher than the current radius.
                while (tfs[tf_index].target_radius <= r) {
                    tf_index += 1;
                    if (tf_index == tfs.len) break :radial;
                }

                const tf1 = tfs[tf_index - 1];
                const tf2 = tfs[tf_index];

                const r1 = tf1.target_radius;
                const r2 = tf2.target_radius;

                // Radial interpolation weight.
                const r_weight = interpolations.lerpWeight(T.T, r, r1, r2);

                const em = opts.emissivity.radialEmissivity(.promote(r)).x;

                const g_min = interpolations.lerpValue(
                    T.T,
                    r_weight,
                    tf1.g_min,
                    tf2.g_min,
                );
                const g_max = interpolations.lerpValue(
                    T.T,
                    r_weight,
                    tf1.g_max,
                    tf2.g_max,
                );

                const weight = r * em * std.math.pi * delta_r / (g_max - g_min);

                energy: for (self.flux_cache, self.g_grid) |*f, g0| {
                    f.* = 0;
                    for (0..opts.g_refine) |n| {
                        const g = g0 + delta_g * @as(T.T, @floatFromInt(n));
                        if (g < g_min) continue;
                        // There is no point continuing further.
                        if (g > g_max) break :energy;

                        const g_star = Table.CTF.gstarFrom(g, g_min, g_max);

                        const branches_1 = tf1.interpolateBranchesGstar(g_star);
                        const branches_2 = tf2.interpolateBranchesGstar(g_star);

                        const lower = interpolations.lerpValue(
                            T.T,
                            r_weight,
                            branches_1.lower,
                            branches_2.lower,
                        );
                        const upper = interpolations.lerpValue(
                            T.T,
                            r_weight,
                            branches_1.upper,
                            branches_2.upper,
                        );

                        const denom = @sqrt(g_star * (1 - g_star));
                        const total_f = lower + upper;
                        // TODO: there is technically a dg term here too but since it's all
                        // equal it can be normalised away later.
                        f.* += total_f * g * g * weight / denom;
                    }
                }

                switch (spectrum) {
                    .no_spectrum => {
                        // Write the output
                        for (flux, self.flux_cache) |*fout, fline| {
                            fout.* += fline;
                        }
                    },
                    .with_spectrum => {
                        const ds = disc_spectrum.?;
                        if (r + delta_r > ds.radii[spectrum_index]) {
                            std.debug.print(">> Convolving at {d}\n", .{r});
                            const spec = ds.fluxes[spectrum_index];

                            // Perform the convolution and sum into the output buffer.
                            spectra.convolve(
                                T.T,
                                flux,
                                self.g_grid,
                                self.flux_cache,
                                ds.energy_grid,
                                spec[0 .. spec.len - 1],
                            );

                            // Advance the disc spectrum.
                            while (r + delta_r > ds.radii[spectrum_index]) {
                                spectrum_index += 1;
                                if (spectrum_index >= ds.radii.len) break :radial;
                            }
                        }
                    },
                }

                r = r_next;
            }
        }
    };
}

/// Calculates the positive difference between two angles.
fn angularDifference(comptime T: type, a1: T, a2: T) T {
    const two_pi = std.math.pi * 2.0;
    const _a1 = @mod(a1, two_pi);
    const _a2 = @mod(a2, two_pi);
    const diff = @mod(_a2 - _a1 + std.math.pi, two_pi) - std.math.pi;
    return @abs(diff);
}

pub const Heuristic = enum {
    /// A heuristic based on the arclengths between two points, refining always
    /// those points that have the largest arclength difference.
    arclength,
    /// Refine those points who have the largest l2 distance between their
    /// impact parameters.
    impact,
    /// Refine based on the energyshift parameter.
    redshift,
    /// Refine based on the Jacobian parameter.
    jacobian,

    /// Parse a heuristic from a string.
    pub fn fromString(s: []const u8) !Heuristic {
        return std.meta.stringToEnum(Heuristic, s) orelse
            return error.NoSuchHeuristic;
    }

    /// Similar to `fromString` but will return `null` if `null` is given.
    pub fn fromOptionalString(string: ?[]const u8) !?Heuristic {
        if (string) |s| return try .fromString(s);
        return null;
    }
};

/// The refinement options for the edges of the transfer function calculations.
pub const Refine = struct {
    /// How many points around the edges to consider when determining the
    /// limits of the refinement.
    N: usize,
    /// How many additional points to trace between the limits.
    M: usize,

    pub fn parse(refine: ?[]const u8) !Refine {
        if (refine) |ref| {
            if (std.mem.indexOfScalar(u8, ref, ',')) |comma| {
                const N = try std.fmt.parseInt(usize, ref[0..comma], 10);
                const M = try std.fmt.parseInt(usize, ref[comma + 1 .. ref.len], 10);
                return .{ .N = N, .M = M };
            }
            return error.MalformedRefine;
        }
        return .{ .N = 0, .M = 0 };
    }
};

fn parseRefine(refine: ?[]const u8) !Refine {
    return Refine.parse(refine);
}

/// The context used for calculating Cunningham's transfer function.
///
/// This is an implementation detail, and is not part of the public interface.
/// It exists principally to seperate the logic that steers how the transfer
/// functions are computed in terms of the spacing and heuristics.
///
/// It is solely used by the `CunninghamTransferFunctionTable`.
fn TransferFunctionContext(comptime T: type) type {
    return struct {
        const TransferFunction = CunninghamTransferFunction(T);
        const Tracer = TransferFunction.Tracer;

        const eps = 1e-6;

        pub const Options = struct {
            minimum_guess: T.T = 5.0,
            max_points: usize = 200,
            initial_guess: ?T.T = null,
            heuristic: ?Heuristic = null,
            refine: Refine = .{ .N = 0, .M = 0 },
            optimise: usize = 17,

            fn toHeuristicOptions(self: Options, initial_guess: T.T) TransferFunction.HeuristicOptions {
                return .{
                    .max_points = self.max_points,
                    .initial_guess = initial_guess,
                    .minimum_guess = self.minimum_guess,
                };
            }
        };

        /// Calculate a transfer function for a particular radius.
        fn calculateRadius(
            allocator: std.mem.Allocator,
            tracer: Tracer,
            target_radius: T.T,
            args: Options,
        ) !TransferFunction {
            var tf = try calculateRadiusImpl(allocator, tracer, target_radius, args);
            try tf.optimiseExtremaAlloc(allocator, tracer, args.optimise);
            return tf;
        }

        fn calculateRadiusImpl(
            allocator: std.mem.Allocator,
            tracer: Tracer,
            target_radius: T.T,
            args: Options,
        ) !TransferFunction {
            const init_guess = args.initial_guess orelse target_radius;
            if (args.heuristic) |heuristic| {
                switch (heuristic) {
                    inline else => |h| {
                        const _heuristicTrace = switch (h) {
                            .redshift => redshiftHeuristic,
                            .jacobian => jacobianHeuristic,
                            .arclength => arclengthHeuristic,
                            .impact => impactHeuristic,
                        };
                        return try TransferFunction.traceAllHeuristicAlloc(
                            allocator,
                            tracer,
                            target_radius,
                            _heuristicTrace,
                            args.toHeuristicOptions(init_guess),
                        );
                    },
                }
            } else {
                var angle_range = RangeIterator(T.T).init(
                    eps,
                    std.math.pi * 2.0 - eps,
                    args.max_points,
                );
                const angles = try angle_range.drain(allocator);
                defer allocator.free(angles);

                const opts: TransferFunction.TracerOptions = .{
                    .initial_guess = init_guess,
                    .minimum_guess = args.minimum_guess,
                };
                return try TransferFunction.traceAllAlloc(
                    allocator,
                    tracer,
                    angles,
                    target_radius,
                    opts,
                );
            }
        }

        fn arclengthHeuristic(left: TransferFunction.Trace, right: TransferFunction.Trace) T.T {
            const delta_theta =
                angularDifference(T.T, left.image_angle, right.image_angle);
            const average_r = (left.image_radius + right.image_radius) / 2.0;
            const score = delta_theta * average_r;
            return score;
        }

        fn impactHeuristic(left: TransferFunction.Trace, right: TransferFunction.Trace) T.T {
            const delta_beta = left.beta - right.beta;
            const delta_alpha = left.alpha - right.alpha;
            return delta_beta * delta_beta + delta_alpha * delta_alpha;
        }

        fn redshiftHeuristic(left: TransferFunction.Trace, right: TransferFunction.Trace) T.T {
            return @abs(left.g - right.g);
        }

        fn jacobianHeuristic(left: TransferFunction.Trace, right: TransferFunction.Trace) T.T {
            return @abs(left.jacobian - right.jacobian);
        }
    };
}

pub fn TableOfTables(comptime T: type) type {
    return struct {
        const Self = @This();
        pub const Table = CunninghamTransferFunctionTable(T);

        /// This determines whether the memory is owned by this table of tables
        /// or not.
        needs_free: bool = true,

        /// In the storage scheme, this is the first or slower varying axis.
        spins: []T.T,
        /// In the storage scheme, this is the second or faster varying axis.
        observer_incls: []T.T,
        /// All of the transfer function tables.
        tables: []Table,

        /// The cache for the interpolated table.
        table: Table,

        pub fn init(allocator: std.mem.Allocator, spins: []T.T, incls: []T.T, tables: []Table) !Self {
            // TODO: check that the dimensions of all of the tables are the same.

            // Create a dupe of the first table to use as the interpolation
            // cache.
            var duped_table = try tables[0].dupe(allocator);
            errdefer duped_table.deinit(allocator);

            return .{
                .spins = spins,
                .observer_incls = incls,
                .tables = tables,
                .table = duped_table,
            };
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            for (self.tables) |*table| {
                table.deinit(allocator);
            }
            self.table.deinit(allocator);
            allocator.free(self.spins);
            allocator.free(self.observer_incls);
            allocator.free(self.tables);
        }

        const Interpolator = dinterp.HypercubeInterpolator(2, T.T);

        /// Interpolate a transfer function table for a particular spin and
        /// inclination. The table memory is owned by the table of tables.
        ///
        /// This function is not thread safe.
        ///
        /// Use `interpolateAlloc` to allocate a new table which the caller
        /// owns. Use `interpolateInto` to interpolate into a pre-allocated
        /// table.
        pub fn interpolate(self: *Self, spin: T, incl: T) Table {
            self.interpolateInto(&self.table, spin, incl);
            return self.table;
        }

        /// Same as `interpolate`, but the table returned is owned by the
        /// caller.
        pub fn interpolateAlloc(
            self: *const Self,
            allocator: std.mem.Allocator,
            spin: T,
            incl: T,
        ) !Table {
            var out = try self.table.dupe(allocator);
            errdefer out.deinit(allocator);
            self.interpolateInto(&out, spin, incl);
            return out;
        }

        /// Same as `interpolate` but the output memory is provided by the
        /// caller. This is threadsafe, provided each thread has it's own `out`
        /// table.
        pub fn interpolateInto(self: *const Self, out: *Table, spin: T, incl: T) void {
            const interp = Interpolator.interpolate(
                .fromSlices(.{ self.spins, self.observer_incls }),
                .{ spin.x, incl.x },
            );
            // Zero the interpolated cache table
            self.zeroInitInterpolatedTable(out, spin.x, incl.x);
            // Do the interpolation
            interp.applyContext(Table, out, self.tables, applyInterpolation);
        }

        fn zeroInitInterpolatedTable(self: *const Self, out: *Table, spin: T.T, incl: T.T) void {
            // Clamp the interpolation variables:
            const _spin = std.math.clamp(
                spin,
                self.spins[0],
                self.spins[self.spins.len - 1],
            );
            const _incl = std.math.clamp(
                incl,
                self.observer_incls[0],
                self.observer_incls[self.observer_incls.len - 1],
            );

            // Zero the cached table before interpolating.
            for (out.transfer_functions.items) |*tf| {
                tf.g_min = 0;
                tf.g_max = 0;
                tf.target_radius = 0;
                for (tf.traces) |*trace| {
                    trace.* = std.mem.zeroes(Table.CTF.Trace);
                }
            }

            // TODO: use the correct disc tracer
            out.tracer = .init(
                .init(.one, .promote(_spin)),
                .{
                    .t = .zero,
                    .r = self.tables[0].tracer.x_obs.r,
                    .th = .promote(std.math.degreesToRadians(_incl)),
                    .ph = .zero,
                },
                .equatorial_plane,
            );
        }

        fn applyInterpolation(out: *Table, op: Interpolator.Op, tables: []const Table) void {
            const table_tfs = tables[op.index].transfer_functions.items;
            for (out.transfer_functions.items, table_tfs) |*out_tf, tf| {
                out_tf.g_min += tf.g_min * op.weight;
                out_tf.g_max += tf.g_max * op.weight;
                out_tf.target_radius += tf.target_radius * op.weight;

                // TODO: maybe skip the inactive fields?
                for (out_tf.traces, tf.traces) |*out_trace, trace| {
                    inline for (@typeInfo(Table.CTF.Trace).@"struct".fields) |field| {
                        @field(out_trace, field.name) +=
                            op.weight * @field(trace, field.name);
                    }
                }
            }
        }
    };
}

/// Read a single transfer function from a file.
pub fn readSingleFromFile(
    comptime T: type,
    allocator: std.mem.Allocator,
    file_path: []const u8,
) !TableOfTables(T).Table {
    var fits = try zfits.FitsFile.open(allocator, .{ .path = file_path });
    defer fits.deinit();
    return parseSingleFromFITS(T, allocator, fits);
}

/// Parse a single transfer function from an opened FITS file.
pub fn parseSingleFromFITS(
    comptime T: type,
    allocator: std.mem.Allocator,
    fits: *zfits.FitsFile,
) !TableOfTables(T).Table {
    const Tables = TableOfTables(T);
    const table_type = fits.hdus[0].getRecord("KZTYPE").?.value.string;
    std.debug.assert(std.mem.eql(u8, table_type, "ctf"));
    // Single CTF table
    return try Tables.Table.fromSingleFITS(allocator, &fits.hdus[1]);
}

/// Read a grid of transfer function tables from a file.
///
pub fn readFromFile(
    comptime T: type,
    allocator: std.mem.Allocator,
    file_path: []const u8,
) !TableOfTables(T) {
    var fits = try zfits.FitsFile.open(allocator, .{ .path = file_path });
    defer fits.deinit();
    return parseFromFITS(T, allocator, fits);
}

/// Read a grid of transfer function tables from an opened FITS file.
///
/// TODO: maybe make this lazy loading?
pub fn parseFromFITS(
    comptime T: type,
    allocator: std.mem.Allocator,
    fits: *zfits.FitsFile,
) !TableOfTables(T) {
    const Tables = TableOfTables(T);

    // First need to determine what kind of table we have. Is it a single
    // transfer function table or a grid of them?
    const table_type = fits.hdus[0].getRecord("KZTYPE").?.value.string;
    std.debug.assert(std.mem.eql(u8, table_type, "ctfgrid"));

    // Read the spins and inclination grids.
    const spin_table = &fits.hdus[1].data.binary_table;
    const incl_table = &fits.hdus[2].data.binary_table;

    const spins = try allocator.alloc(T.T, spin_table.num_rows);
    errdefer allocator.free(spins);

    for (spins, 0..) |*spin, row_index| {
        const row = try spin_table.getOrParseRow(allocator, row_index);
        spin.* = @floatCast(row.cols[0].one.float_32);
    }

    const incls = try allocator.alloc(T.T, incl_table.num_rows);
    errdefer allocator.free(incls);

    for (incls, 0..) |*incl, row_index| {
        const row = try incl_table.getOrParseRow(allocator, row_index);
        incl.* = @floatCast(row.cols[0].one.float_32);
    }

    // Now read all of the tables:
    const tables = try allocator.alloc(Tables.Table, spins.len * incls.len);
    errdefer allocator.free(tables);

    var tables_read: usize = 0;
    errdefer for (tables[0..tables_read]) |*table| table.deinit(allocator);

    for (tables, 3..) |*table, hdu_index| {
        table.* = try Tables.Table.fromSingleFITS(allocator, &fits.hdus[hdu_index]);
        tables_read += 1;
    }

    return try .init(allocator, spins, incls, tables);
}
