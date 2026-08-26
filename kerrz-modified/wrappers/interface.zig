/// This file provides functional bindings for using kerrz in other languages.
const std = @import("std");
const kerrz = @import("kerrz");
const c = @cImport({
    @cInclude("kerrz.h");
});

const Dual0 = kerrz.DualNumber(f64, 0);
const EmissivityCache = kerrz.EmissivityCache(Dual0);
const allocator = std.heap.c_allocator;
const CTF = kerrz.transfer_tables.CunninghamTransferFunction(Dual0);

export fn krz_ThreadPool_init(pool: *c.krz_ThreadPool, n_threads: usize) callconv(.c) c_int {
    const num_threads = if (n_threads == 0) std.Thread.getCpuCount() catch 1 else n_threads;
    const threads = kerrz.ThreadMap.init(allocator, .{ .num_threads = num_threads }) catch {
        return c.RETCODE_ALLOCATION_FAILED;
    };
    pool.pool = @ptrCast(threads);
    pool.num_threads = num_threads;
    return c.RETCODE_SUCCESS;
}

export fn krz_ThreadPool_deinit(pool: *c.krz_ThreadPool) callconv(.c) void {
    const threads: *kerrz.ThreadMap = @ptrCast(@alignCast(pool.pool));
    threads.deinit();
}

export fn krz_kerrMetric(M: f64, a: f64) callconv(.c) c.krz_KerrMetric {
    const metric = kerrz.KerrMetric(Dual0).init(
        .promote(M),
        // TODO: remove me when a=0 is implemented.
        .promote(if (std.math.approxEqAbs(f64, a, 0, 1e-3)) 1e-3 else a),
    );
    return .{
        .M = metric.M.x,
        .a = metric.a.x,
        .horizon_radius = metric.horizon_radius.x,
        .horizon_radius_negative = metric.horizon_radius_negative.x,
        .isco = metric.isco.x,
    };
}

fn unwrap(comptime T: type, comptime DualType: type, a: anytype) T {
    var out: T = undefined;
    const info = @typeInfo(T).@"struct";
    inline for (info.fields) |f| {
        @field(out, f.name) = DualType.promote(@field(a, f.name));
    }
    return out;
}

fn packInitialConditions(comptime T: type, geod: kerrz.NullGeodesic(T)) c.krz_InitialConditions {
    // Assign values back.
    // TODO: use memcpy or related
    return .{
        .E = geod.E.x,
        .L = geod.L.x,
        .Q = geod.Q.x,
        .lambda = geod.lambda.x,
        .eta = geod.eta.x,
        .x_init = .{
            .t = geod.x_init.t.x,
            .th = geod.x_init.th.x,
            .ph = geod.x_init.ph.x,
            .r = geod.x_init.r.x,
        },
        .theta_sign = @floatFromInt(geod.theta_sign),
        .radial_sign = @floatFromInt(geod.radial_sign),
        .windings = @floatFromInt(geod.windings),
    };
}

export fn krz_fromImpactParameters(
    metric: c.krz_KerrMetric,
    x_init: c.krz_FourVector,
    alpha: f64,
    beta: f64,
) callconv(.c) c.krz_InitialConditions {
    const _metric = unwrap(kerrz.KerrMetric(Dual0), Dual0, metric);
    const _x_init = unwrap(kerrz.FourVector(Dual0), Dual0, x_init);
    const geod = kerrz.NullGeodesic(Dual0).fromImpactParameters(
        _metric,
        _x_init,
        .promote(alpha),
        .promote(beta),
    );
    return packInitialConditions(Dual0, geod);
}

fn translateStatusCode(status: kerrz.Status) c.krz_STATUS {
    return switch (status) {
        .event_horizon => c.STATUS_EVENT_HORIZON,
        .infinity => c.STATUS_AT_INFINITY,
        .intersected_disc => c.STATUS_INTERSECTED_DISC,
        .no_status => c.STATUS_NONE,
    };
}

fn unwrapInitialConditions(
    comptime T: type,
    init_conds: c.krz_InitialConditions,
) kerrz.NullGeodesic(T) {
    return .{
        .E = .promote(init_conds.E),
        .L = .promote(init_conds.L),
        .Q = .promote(init_conds.Q),
        .lambda = .promote(init_conds.lambda),
        .eta = .promote(init_conds.eta),
        .x_init = .{
            .t = .promote(init_conds.x_init.t),
            .th = .promote(init_conds.x_init.th),
            .ph = .promote(init_conds.x_init.ph),
            .r = .promote(init_conds.x_init.r),
        },
        .theta_sign = @intFromFloat(init_conds.theta_sign),
        .radial_sign = @intFromFloat(init_conds.radial_sign),
        .windings = @intFromFloat(init_conds.windings),
    };
}

fn wrapResult(
    comptime T: type,
    metric: kerrz.KerrMetric(T),
    geod: kerrz.NullGeodesic(T),
    result: kerrz.TraceResult(T),
) c.krz_TraceResult {
    const x_final: kerrz.FourVector(T) = switch (result.status) {
        .event_horizon => .zeros,
        else => b: {
            const total = result.totalAntiderivatives(metric, geod);
            break :b .{
                .t = total.coordinateTime(metric, geod),
                .r = result.r,
                .th = result.theta,
                .ph = total.coordinateAzimuth(metric, geod),
            };
        },
    };
    return .{
        .status = translateStatusCode(result.status),
        .mino_time = result.mino_time.x,
        .x_final = .{
            .t = x_final.t.x,
            .r = x_final.r.x,
            .th = x_final.th.x,
            .ph = x_final.ph.x,
        },
        .winding = @floatFromInt(result.winding),
    };
}

export fn krz_traceToAngle(
    metric: c.krz_KerrMetric,
    init_conds: c.krz_InitialConditions,
    angle: f64,
) callconv(.c) c.krz_TraceResult {
    const _metric = unwrap(kerrz.KerrMetric(Dual0), Dual0, metric);
    const _geod = unwrapInitialConditions(Dual0, init_conds);
    const result = _geod.traceToAngle(_metric, .promote(angle), .{});
    return wrapResult(Dual0, _metric, _geod, result);
}

export fn krz_traceToRadius(
    metric: c.krz_KerrMetric,
    init_conds: c.krz_InitialConditions,
    radius: f64,
) callconv(.c) c.krz_TraceResult {
    const _metric = unwrap(kerrz.KerrMetric(Dual0), Dual0, metric);
    const _geod = unwrapInitialConditions(Dual0, init_conds);
    const result = _geod.traceToRadius(_metric, .promote(radius), .{});
    return wrapResult(Dual0, _metric, _geod, result);
}

const EmissivityContext = struct {
    cache: EmissivityCache,
    bin_cache: EmissivityCache.BinningCache,
    table: ?kerrz.emissivity.TableEmissivity(Dual0.T) = null,
};

export fn krz_EmissivityCache_init(
    cache: *c.krz_EmissivityCache,
    num_traces: usize,
) callconv(.c) c.krz_RETCODE {
    const ctx_pointer = allocator.create(EmissivityContext) catch {
        return c.RETCODE_ALLOCATION_FAILED;
    };

    var em_cache = EmissivityCache.init(
        allocator,
        .{ .max_traces = num_traces },
    ) catch {
        allocator.destroy(ctx_pointer);
        return c.RETCODE_ALLOCATION_FAILED;
    };

    const bin_cache = EmissivityCache.BinningCache.init(
        allocator,
        .{},
    ) catch {
        em_cache.deinit(allocator);
        allocator.destroy(ctx_pointer);
        return c.RETCODE_ALLOCATION_FAILED;
    };

    ctx_pointer.* = .{
        .cache = em_cache,
        .bin_cache = bin_cache,
    };

    cache.context = @ptrCast(ctx_pointer);

    return c.RETCODE_SUCCESS;
}

export fn krz_EmissivityCache_deinit(
    ctx: *c.krz_EmissivityCache,
) void {
    const context: *EmissivityContext = @ptrCast(@alignCast(ctx.context));
    context.cache.deinit(allocator);
    context.bin_cache.deinit(allocator);
    allocator.destroy(context);
}

const CoronalModel = union(enum) {
    ring: kerrz.emissivity.Ring(Dual0),
    lp: kerrz.emissivity.Lamppost(Dual0),

    /// Convert the C interface coronal abstraction to a Zig native one.
    fn fromC(
        metric: c.krz_KerrMetric,
        model: c.krz_CoronaModel,
    ) CoronalModel {
        switch (model.tag) {
            c.RING_CORONA => {
                return .{
                    .ring = .init(
                        unwrap(kerrz.KerrMetric(Dual0), Dual0, metric),
                        .{
                            .height = .promote(model.as.ring.height),
                            .radius = .promote(model.as.ring.radius),
                        },
                    ),
                };
            },
            else => unreachable,
        }
    }
};

export fn krz_emissivity_ring(
    pool: *c.krz_ThreadPool,
    cache: *c.krz_EmissivityCache,
    metric: c.krz_KerrMetric,
    ring: c.krz_RingCorona,
) callconv(.c) c.krz_RETCODE {
    return krz_emissivity(
        pool,
        cache,
        metric,
        .{
            .tag = c.RING_CORONA,
            .as = .{ .ring = ring },
        },
    );
}

export fn krz_emissivity(
    pool: *c.krz_ThreadPool,
    cache: *c.krz_EmissivityCache,
    metric: c.krz_KerrMetric,
    model: c.krz_CoronaModel,
) callconv(.c) c.krz_RETCODE {
    const ctx: *EmissivityContext = @ptrCast(@alignCast(cache.context));
    const threads: *kerrz.ThreadMap = @ptrCast(@alignCast(pool.pool));

    const corona = CoronalModel.fromC(metric, model);

    const ring_ctx: RingThreadContext = .{
        .ctx = ctx,
        .ring = corona.ring,
    };

    // TODO: make the emissivity tool compatible with an in-place calculation
    // to avoid repeating code here.

    threads.map(
        kerrz.EmissivityTrace(Dual0),
        ctx.cache.traces,
        ring_ctx,
        RingThreadContext.tracer,
        .{
            .chunk_size = 4096,
        },
    ) catch {
        return c.RETCODE_THREAD_ERROR;
    };
    threads.blockUntilDone();

    // Sort them all by radius, then can perform the 2d rebinning without
    // needing to a 2d grid.
    ctx.cache.sortByRadius();
    ctx.table = ctx.cache.rebinTimeDependentInplace(&ctx.bin_cache, .keep);
    ctx.table.?.interpolateMissing();

    return c.RETCODE_SUCCESS;
}

const RingThreadContext = struct {
    ring: kerrz.emissivity.Ring(Dual0),
    ctx: *EmissivityContext,

    fn tracer(
        ctx: RingThreadContext,
        out: *kerrz.EmissivityTrace(Dual0),
        index: usize,
        _: kerrz.ThreadMap.ThreadId,
    ) void {
        out.* = ctx.ctx.cache.traceRing(ctx.ring, index);
        ctx.ctx.cache.traces[index] = out.*;
    }
};

export fn krz_interpolate_emissivity(
    cache: *c.krz_EmissivityCache,
    r: f64,
    phi: f64,
) callconv(.c) c.krz_EmissivityTrace {
    const context: *EmissivityContext = @ptrCast(@alignCast(cache.context));
    const trace = context.table.?.values(r, phi);
    return .{
        .t = trace.t,
        .g = trace.g,
        .em = trace.em,
        .local_theta = trace.local_theta,
    };
}

fn projection_trace(
    metric: kerrz.KerrMetric(Dual0),
    obs: kerrz.FourVector(Dual0),
    r_disc: Dual0,
    angle: Dual0,
    initial_guess: Dual0.T,
) !struct {
    alpha: Dual0.T,
    beta: Dual0.T,
    phi: Dual0.T,
    r_offset: Dual0.T,
    time: Dual0.T,
} {
    const sol = try kerrz.transfer_tables.impactOffsetForRadius(
        Dual0,
        metric,
        obs,
        .promote(angle),
        .promote(r_disc),
        .{ .initial_guess = initial_guess },
    );

    const params = kerrz.transfer_tables.toImpactParameters(
        Dual0,
        .promote(sol.r_offset),
        .promote(sol.theta_image_plane),
    );

    const geod = kerrz.NullGeodesic(Dual0).fromImpactParameters(
        metric,
        obs,
        params.alpha,
        params.beta,
    );
    const res = geod.traceToAngle(metric, .promote(std.math.pi / 2.0), .{});
    const total = res.totalAntiderivatives(metric, geod);
    const phi = total.coordinateAzimuth(metric, geod);
    const time = total.coordinateTime(metric, geod);

    return .{
        .alpha = params.alpha.x,
        .beta = params.beta.x,
        .phi = phi.x,
        .r_offset = sol.r_offset,
        .time = time,
    };
}

const TransferFunctionTable = kerrz.transfer_tables.CunninghamTransferFunctionTable(Dual0);

const CunninghamCache = struct {
    arena: std.heap.ArenaAllocator,
    table: TransferFunctionTable,
    radii: []const f64,
};

export fn krz_TransferFunctionCache_free(cache: *c.krz_TransferFunctionCache) callconv(.c) void {
    const ctf_cache: *CunninghamCache = @ptrCast(@alignCast(cache.cache));
    ctf_cache.arena.deinit();
}

export fn krz_transfer_functions(
    cache: *c.krz_TransferFunctionCache,
    metric: c.krz_KerrMetric,
    observer: c.krz_FourVector,
    r_min: f64,
    r_max: f64,
    num_radii: usize,
) c.krz_STATUS {
    const ctf_cache = allocator.create(CunninghamCache) catch {
        return c.RETCODE_ALLOCATION_FAILED;
    };

    const dual_metric = unwrap(kerrz.KerrMetric(Dual0), Dual0, metric);
    const dual_observer: kerrz.FourVector(Dual0) = .{
        .t = .promote(observer.t),
        .r = .promote(observer.r),
        .th = .promote(observer.th),
        .ph = .promote(observer.ph),
    };

    ctf_cache.* = calculate_transfer_functions(
        dual_metric,
        dual_observer,
        r_min,
        r_max,
        num_radii,
    ) catch {
        krz_TransferFunctionCache_free(cache);
        return c.RETCODE_CONVERGENCE_FAILED;
    };

    cache.* = .{ .cache = ctf_cache };

    return c.RETCODE_SUCCESS;
}

fn sortTraceAngle(_: void, lhs: CTF.Trace, rhs: CTF.Trace) bool {
    return lhs.phi < rhs.phi;
}

const CTFTable = kerrz.CunninghamTransferFunctionTable(Dual0);

fn calculate_transfer_functions(
    metric: kerrz.KerrMetric(Dual0),
    x_obs: kerrz.FourVector(Dual0),
    r_min: f64,
    r_max: f64,
    n_radii: usize,
) !CunninghamCache {
    var arena = std.heap.ArenaAllocator.init(allocator);
    errdefer arena.deinit();
    const alloc = arena.allocator();

    var table = TransferFunctionTable.init(metric, x_obs, .equatorial_plane);
    errdefer table.deinit(alloc);

    const radii = try alloc.alloc(f64, n_radii);
    errdefer alloc.free(radii);

    const options: CTFTable.Options = .{
        .max_points = 200,
        .optimise = 17,
        .minimum_guess = 5.0,
    };

    const log_r_max = std.math.log10(r_max);
    const log_r_min = std.math.log10(r_min);

    for (0..n_radii) |i| {
        const f = @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(n_radii - 1));
        const log_r = (log_r_max - log_r_min) * f + log_r_min;
        const r = std.math.pow(f64, 10, log_r);

        radii[i] = r;
        const ptr = try table.calculateRadiusAndAppend(allocator, r, options);
        ptr.normaliseAngles();
    }

    return .{
        .arena = arena,
        .table = table,
        .radii = radii,
    };
}

export fn krz_interpolate_disc_coordinates(
    cache: *c.krz_TransferFunctionCache,
    r: f64,
    phi: f64,
) c.krz_InterpolatedPoint {
    const ctf_cache: *CunninghamCache = @ptrCast(@alignCast(cache.cache));

    const phi_mod = @mod(phi, std.math.pi * 2.0);
    const trace = ctf_cache.table.interpolateDiscCoordinates(r, phi_mod);

    return .{
        .alpha = trace.alpha,
        .beta = trace.beta,
        .delta_t = trace.delta_t,
    };
}

export fn krz_shadow(
    metric: c.krz_KerrMetric,
    obs_theta: f64,
    alpha: [*]f64,
    beta: [*]f64,
    num_points: usize,
) c.krz_STATUS {
    const dual_metric = unwrap(kerrz.KerrMetric(Dual0), Dual0, metric);
    calc_shadow(
        dual_metric,
        .promote(obs_theta),
        alpha[0..num_points],
        beta[0..num_points],
    ) catch {
        return c.RETCODE_CONVERGENCE_FAILED;
    };
    return c.RETCODE_SUCCESS;
}

fn calc_shadow(
    metric: kerrz.KerrMetric(Dual0),
    obs_theta: Dual0,
    alpha: []f64,
    beta: []f64,
) !void {
    const dual_observer: kerrz.FourVector(Dual0) = .{
        .t = .zero,
        .r = .promote(1e6),
        .th = obs_theta,
        .ph = .zero,
    };

    const Context = struct {
        m: kerrz.KerrMetric(Dual0),
        x: kerrz.FourVector(Dual0),
        theta: Dual0.T,

        const Self = @This();

        pub fn impactParameters(self: Self, r: Dual0.T) struct {
            alpha: f64,
            beta: f64,
        } {
            return .{
                .alpha = r * @cos(self.theta),
                .beta = r * @sin(self.theta),
            };
        }

        pub fn isHorizon(self: Self, r: Dual0.T) bool {
            const params = self.impactParameters(r);

            const geod = kerrz.NullGeodesic(Dual0).fromImpactParameters(
                self.m,
                self.x,
                .promote(params.alpha),
                .promote(params.beta),
            );

            const radial_roots = kerrz.potentials.rootsOfRadialPotential(
                Dual0,
                self.m,
                geod.eta,
                geod.lambda,
            );

            const radial_case = radial_roots.determineCase(
                self.m.horizon_radius,
                geod.x_init.r,
            );

            return radial_case == .case_III;
        }
    };

    const num_r: usize = 2000;
    const max_r: f64 = 9.0;

    outer: for (0..alpha.len) |i| {
        const f = @as(f64, @floatFromInt(i)) /
            @as(f64, @floatFromInt(alpha.len));
        const theta = f * std.math.pi * 2.0 + 1e-4;
        const ctx = Context{
            .m = metric,
            .x = dual_observer,
            .theta = theta,
        };

        for (0..num_r) |j| {
            const fr = @as(f64, @floatFromInt(j)) /
                @as(f64, @floatFromInt(num_r));
            const r = (1 - fr) * max_r;
            if (ctx.isHorizon(r)) {
                const params = ctx.impactParameters(r);
                alpha[i] = params.alpha;
                beta[i] = params.beta;
                continue :outer;
            }
        }

        unreachable;
    }
}

const PathBuilderState = struct {
    const PathBuilder = kerrz.NullGeodesic(Dual0).PathBuilder;

    state: PathBuilder,
    geodesic: kerrz.NullGeodesic(Dual0),
    metric: kerrz.KerrMetric(Dual0),

    pub fn init(
        metric: kerrz.KerrMetric(Dual0),
        geod: kerrz.NullGeodesic(Dual0),
    ) PathBuilderState {
        return .{
            .state = geod.traceBuilder(metric, .{}),
            .geodesic = geod,
            .metric = metric,
        };
    }
};

export fn krz_PathBuilder_init(
    pb: *c.krz_PathBuilder,
    metric: c.krz_KerrMetric,
    init_conds: c.krz_InitialConditions,
) callconv(.c) c.krz_STATUS {
    const ptr = allocator.create(PathBuilderState) catch
        return c.RETCODE_ALLOCATION_FAILED;

    const _metric = unwrap(kerrz.KerrMetric(Dual0), Dual0, metric);
    const _geod = unwrapInitialConditions(Dual0, init_conds);
    ptr.* = PathBuilderState.init(_metric, _geod);
    pb.state = @ptrCast(ptr);
    return c.RETCODE_SUCCESS;
}

export fn krz_PathBuilder_deinit(
    pb: *c.krz_PathBuilder,
) callconv(.c) void {
    const ptr: *PathBuilderState = @ptrCast(@alignCast(pb.state));
    allocator.destroy(ptr);
    pb.state = null;
}

export fn krz_at_mino_time(
    pb: *c.krz_PathBuilder,
    mino_time: f64,
) callconv(.c) c.krz_TraceResult {
    const ptr: *PathBuilderState = @ptrCast(@alignCast(pb.state));
    const result = ptr.state.atMinoTime(.promote(mino_time));
    return wrapResult(Dual0, ptr.metric, ptr.geodesic, result);
}

export fn krz_mino_time_to_turning_points(
    pb: *c.krz_PathBuilder,
) callconv(.c) c.krz_TurningPoints {
    const ptr: *PathBuilderState = @ptrCast(@alignCast(pb.state));
    return .{
        .theta_0 = ptr.state.theta_0_time.x,
        .theta_1 = ptr.state.theta_1_time.x,
        .r_0 = ptr.state.r_0_time.x,
        .r_1 = ptr.state.r_1_time.x,
    };
}

export fn krz_tool_TransferFunction_defaults() callconv(.c) c.krz_tool_TransferFunction {
    return .{
        .tag = c.TOOL_TRANSFER_FUNCTION,
        .metric = krz_kerrMetric(1.0, 0.998),
        .x_obs = .{
            .t = 0,
            .r = 1e6,
            .th = std.math.degreesToRadians(60),
            .ph = 0,
        },
        .r_target = 10.0,
        .options = .{
            .max_points = 200,
            .refine_N = 0,
            .refine_M = 0,
            .optimise = 17,
            .minimum_guess = 5.0,
            .heuristic = c.TF_HEURISTIC_NONE,
        },
    };
}
export fn krz_CunninghamTransferFunction_deinit(
    ctf: *c.krz_CunninghamTransferFunction,
) callconv(.c) void {
    const traces = ctf.traces[0..ctf.num_traces];
    allocator.free(traces);
}

export fn krz_tool_TransferFunction_run(
    tool: c.krz_tool_TransferFunction,
    result: *c.krz_CunninghamTransferFunction,
) callconv(.c) c.krz_RETCODE {
    const Tool = kerrz.tools.TransferFunction(Dual0);
    const _metric = unwrap(kerrz.KerrMetric(Dual0), Dual0, tool.metric);
    const _x_obs = unwrap(kerrz.FourVector(Dual0), Dual0, tool.x_obs);

    const _tool: Tool = .{
        .metric = _metric,
        .x_obs = _x_obs,
        .r_target = tool.r_target,
        .opts = .{
            .max_points = tool.options.max_points,
            .optimise = tool.options.optimise,
            .refine = .{
                .M = tool.options.refine_M,
                .N = tool.options.refine_N,
            },
            .heuristic = switch (tool.options.heuristic) {
                c.TF_HEURISTIC_NONE => null,
                c.TF_HEURISTIC_ARCLEN => .arclength,
                c.TF_HEURISTIC_IMPACT => .impact,
                else => unreachable,
            },
            .minimum_guess = tool.options.minimum_guess,
        },
    };

    var _result = _tool.run(allocator) catch {
        // TODO: handle these properly
        return c.RETCODE_ALLOCATION_FAILED;
    };
    defer _result.deinit(allocator);

    const traces = allocator.alloc(c.krz_CunninghamTrace, _result.ctf.traces.len) catch {
        return c.RETCODE_ALLOCATION_FAILED;
    };

    // By doing this copy explicitly instead of casting the pointer the
    // implementation detail of the trace can be decoupled from what is in the
    // wrapper.
    for (traces, _result.ctf.traces) |*out, t| {
        out.* = .{
            .image_angle = t.image_angle,
            .image_radius = t.image_radius,
            .alpha = t.alpha,
            .beta = t.beta,
            .g = t.g,
            .g_star = t.g_star,
            .phi = t.phi,
            .delta_t = t.delta_t,
            .jacobian = t.jacobian,
            .f = t.f,
            .r_err = t.r_err,
        };
    }

    result.* = .{
        .g_min = _result.ctf.g_min,
        .g_max = _result.ctf.g_max,
        .traces = traces.ptr,
        .num_traces = traces.len,
    };

    return c.RETCODE_SUCCESS;
}

export fn krz_stationaryFrame(
    metric: c.krz_KerrMetric,
    x_init: c.krz_FourVector,
) callconv(.c) c.krz_OrthonormalFrame {
    const m = unwrap(kerrz.KerrMetric(Dual0), Dual0, metric);
    const x = unwrap(kerrz.FourVector(Dual0), Dual0, x_init);
    const ts = m.tangentSpace(x);
    const v_src = kerrz.orbits.stationary(Dual0, ts);
    const frame = ts.localFrame(v_src);
    return .{
        .x = x_init,
        .metric_components = .{
            ts.metric_components.tt.x,
            ts.metric_components.rr.x,
            ts.metric_components.thth.x,
            ts.metric_components.phph.x,
            ts.metric_components.tph.x,
        },
        .matrix = @bitCast(frame.m.v),
    };
}

export fn krz_frame(
    metric: c.krz_KerrMetric,
    x_init: c.krz_FourVector,
    v_frame: c.krz_FourVector,
) callconv(.c) c.krz_OrthonormalFrame {
    const m = unwrap(kerrz.KerrMetric(Dual0), Dual0, metric);
    const x = unwrap(kerrz.FourVector(Dual0), Dual0, x_init);
    const v = unwrap(kerrz.FourVector(Dual0), Dual0, v_frame);
    const ts = m.tangentSpace(x);
    const frame = ts.localFrame(v);
    return .{
        .x = x_init,
        .metric_components = .{
            ts.metric_components.tt.x,
            ts.metric_components.rr.x,
            ts.metric_components.thth.x,
            ts.metric_components.phph.x,
            ts.metric_components.tph.x,
        },
        .matrix = @bitCast(frame.m.v),
    };
}

export fn krz_lnrFrame(
    metric: c.krz_KerrMetric,
    x_init: c.krz_FourVector,
) callconv(.c) c.krz_OrthonormalFrame {
    const m = unwrap(kerrz.KerrMetric(Dual0), Dual0, metric);
    const x = unwrap(kerrz.FourVector(Dual0), Dual0, x_init);
    const ts = m.tangentSpace(x);
    const v_src = kerrz.orbits.lnr(Dual0, ts);
    const frame = ts.localFrame(v_src);
    return .{
        .x = x_init,
        .metric_components = .{
            ts.metric_components.tt.x,
            ts.metric_components.rr.x,
            ts.metric_components.thth.x,
            ts.metric_components.phph.x,
            ts.metric_components.tph.x,
        },
        .matrix = @bitCast(frame.m.v),
    };
}

export fn krz_circularOrbitVelocity(
    metric: c.krz_KerrMetric,
    r: f64,
) callconv(.c) c.krz_FourVector {
    const m = unwrap(kerrz.KerrMetric(Dual0), Dual0, metric);
    const v = kerrz.orbits.keplerianPlungingAlt(Dual0, m, .promote(r));
    return .{
        .t = v.t.x,
        .r = v.r.x,
        .th = v.th.x,
        .ph = v.ph.x,
    };
}

export fn krz_fromSkyAngles(
    metric: c.krz_KerrMetric,
    frame: c.krz_OrthonormalFrame,
    theta: f64,
    phi: f64,
) callconv(.c) c.krz_InitialConditions {
    const m = unwrap(kerrz.KerrMetric(Dual0), Dual0, metric);
    const x = unwrap(kerrz.FourVector(Dual0), Dual0, frame.x);

    const f: kerrz.geometry.TetradFrame(Dual0) = .{
        .m = .{ .v = @bitCast(frame.matrix) },
    };
    const ts: kerrz.KerrMetric(Dual0).TangentSpace = .{
        .x = x,
        .metric_components = .{
            .tt = .promote(frame.metric_components[0]),
            .rr = .promote(frame.metric_components[1]),
            .thth = .promote(frame.metric_components[2]),
            .phph = .promote(frame.metric_components[3]),
            .tph = .promote(frame.metric_components[4]),
        },
    };

    const velocity = kerrz.geodesic.skyAnglesToVelocityFrame(
        Dual0,
        x,
        f,
        .promote(std.math.clamp(
            @rem(@abs(theta), std.math.pi),
            1e-4,
            std.math.pi - 1e-4,
        )),
        .promote(phi),
    );
    const geod = kerrz.NullGeodesic(Dual0).fromNullVelocityTangentSpace(
        m,
        ts,
        velocity,
    );
    return packInitialConditions(Dual0, geod);
}

export fn krz_traceContinuumLamppost(
    metric: c.krz_KerrMetric,
    x_obs: c.krz_FourVector,
    height: f64,
) callconv(.c) c.krz_ContinuumLamppost {
    const m = unwrap(kerrz.KerrMetric(Dual0), Dual0, metric);
    const x = unwrap(kerrz.FourVector(Dual0), Dual0, x_obs);

    const sol = kerrz.continuum.continuumGeodesic(
        Dual0,
        m,
        x,
        .{ .lamppost = .init(m, .{ .height = .promote(height) }) },
    ) catch |err| {
        const str = std.fmt.allocPrint(allocator, ">> {any}\n", .{err}) catch unreachable;
        defer allocator.free(str);
        @panic(str);
    };

    return .{
        .angle_delta = sol.local_theta.x,
        .dcosd_dcosth = sol.jacobianAsCosine(),
        .res = wrapResult(Dual0, m, sol.geod, sol.result),
        .alpha = sol.alpha.x,
        .beta = sol.beta.x,
    };
}

export fn krz_traceContinuumRing(
    out: *c.krz_ContinuumRing,
    metric: c.krz_KerrMetric,
    x_obs: c.krz_FourVector,
    ring: c.krz_RingCorona,
) callconv(.c) c.krz_RETCODE {
    const dual_metric = unwrap(kerrz.KerrMetric(Dual0), Dual0, metric);
    const dual_observer: kerrz.FourVector(Dual0) = .{
        .t = .promote(x_obs.t),
        .r = .promote(x_obs.r),
        .th = .promote(x_obs.th),
        .ph = .promote(x_obs.ph),
    };
    const corona = CoronalModel.fromC(metric, .{ .tag = c.RING_CORONA, .as = .{ .ring = ring } });
    const ptr = calculate_continuum(dual_metric, dual_observer, corona.ring) catch {
        return c.RETCODE_ALLOCATION_FAILED;
    };

    out.* = .{ .cache = @ptrCast(@alignCast(ptr)) };

    return c.RETCODE_SUCCESS;
}

export fn krz_ContinuumRing_deinit(
    continuum: *c.krz_ContinuumRing,
) callconv(.c) void {
    const data: *kerrz.continuum.ContinuumTransfer(Dual0) = @ptrCast(
        @alignCast(continuum.cache),
    );
    data.deinit(allocator);
}

export fn krz_interpolate_continuum(
    continuum: *c.krz_ContinuumRing,
    phi: f64,
) c.krz_ContinuumRingPoint {
    const data: *kerrz.continuum.ContinuumTransfer(Dual0) = @ptrCast(
        @alignCast(continuum.cache),
    );
    const point = data.interpolateAzimuth(phi);
    const dcosth_dcosd = point.jac *
        (@sin(data.theta_obs.x) / @sin(point.local_theta));
    return .{
        .dcosd_dcosth = 1.0 / dcosth_dcosd,
        .energyshift = point.g,
        .delta_t = point.delta_t,
    };
}

fn calculate_continuum(
    metric: kerrz.KerrMetric(Dual0),
    x_obs: kerrz.FourVector(Dual0),
    ring: kerrz.emissivity.Ring(Dual0),
) !*kerrz.continuum.ContinuumTransfer(Dual0) {
    const transfer = try kerrz.continuum.transferFunction(
        Dual0,
        allocator,
        metric,
        x_obs,
        .{ .ring = ring },
        .{
            .max_points = 200,
        },
    );
    errdefer transfer.deinit(allocator);

    const ptr = try allocator.create(kerrz.continuum.ContinuumTransfer(Dual0));
    ptr.* = transfer;
    return ptr;
}

fn loadAxisymmetricEmissivity(path: []const u8) !kerrz.emissivity.AxisymmetricEmissivity(Dual0) {
    const fits = try kerrz.FitsFile.open(allocator, .{ .path = path });
    defer fits.deinit();
    return try kerrz.emissivity.AxisymmetricEmissivity(Dual0).fromFITS(
        allocator,
        &fits.hdus[1],
    );
}

export fn krz_tool_Lineprofile_defaults() callconv(.c) c.krz_tool_Lineprofile {
    return .{
        .metric = krz_kerrMetric(1.0, 0.998),
        .x_obs = .{
            .t = 0,
            .r = 1e7,
            .th = std.math.degreesToRadians(30),
            .ph = 0,
        },
        .emissivity_path = null,
        .r_in = 0,
        .r_out = 400,
        .nradii = 100,
        .nangles = 200,
        .ng = 1000,
        .ngstar = 2800,
        .nrsteps = 3000,
        .normalise = 1,
    };
}

export fn krz_Lineprofile_deinit(lp: *c.krz_Lineprofile) callconv(.c) void {
    if (lp.g_grid != null and lp.num_g > 0) {
        allocator.free(lp.g_grid[0..lp.num_g]);
    }
    if (lp.flux != null and lp.num_g > 0) {
        allocator.free(lp.flux[0..lp.num_g]);
    }
    lp.* = .{
        .g_grid = null,
        .flux = null,
        .num_g = 0,
    };
}

export fn krz_tool_Lineprofile_run(
    pool: *c.krz_ThreadPool,
    tool: c.krz_tool_Lineprofile,
    result: *c.krz_Lineprofile,
) callconv(.c) c.krz_RETCODE {
    if (tool.emissivity_path == null) {
        return c.RETCODE_INVALID_ARGUMENT;
    }
    const path = std.mem.span(tool.emissivity_path);

    const threads: *kerrz.ThreadMap = @ptrCast(@alignCast(pool.pool));
    const _metric = unwrap(kerrz.KerrMetric(Dual0), Dual0, tool.metric);
    const _x_obs = unwrap(kerrz.FourVector(Dual0), Dual0, tool.x_obs);

    var em_table = loadAxisymmetricEmissivity(path) catch {
        return c.RETCODE_INVALID_ARGUMENT;
    };
    defer em_table.deinit(allocator);

    const TFTableTool = kerrz.tools.TransferFunctionTable(Dual0);
    const tf_tool: TFTableTool = .{
        .metric = _metric,
        .x_obs = _x_obs,
        .table_opts = .{
            .num_radii = tool.nradii,
            .r_in = if (tool.r_in == 0) null else tool.r_in,
            .r_out = tool.r_out,
        },
        .tf_opts = .{
            .max_points = tool.nangles,
        },
    };

    var tf_result = tf_tool.run(allocator, threads, .{}) catch {
        return c.RETCODE_CONVERGENCE_FAILED;
    };
    defer tf_result.deinit(allocator);

    const r_min = if (tool.r_in == 0) _metric.isco.x else tool.r_in;

    const LP = kerrz.tools.Lineprofile(Dual0);
    const lp_tool: LP = .{
        .table = tf_result.table,
        .num_g = tool.ng,
        .num_fine_g = tool.ngstar,
        .lp_opts = .{
            .emissivity = .{ .axisymmetric = em_table },
            .num_r_steps = tool.nrsteps,
            .r_min = r_min,
            .r_max = tool.r_out,
            .normalise = tool.normalise != 0,
        },
    };

    var lp_result = lp_tool.run(allocator, .{}) catch {
        return c.RETCODE_ALLOCATION_FAILED;
    };
    defer lp_result.deinit(allocator);

    const g_grid = allocator.dupe(f64, lp_result.profile.g_grid) catch {
        return c.RETCODE_ALLOCATION_FAILED;
    };
    errdefer allocator.free(g_grid);

    const flux = allocator.dupe(f64, lp_result.profile.flux) catch {
        return c.RETCODE_ALLOCATION_FAILED;
    };

    result.* = .{
        .g_grid = g_grid.ptr,
        .flux = flux.ptr,
        .num_g = g_grid.len,
    };

    return c.RETCODE_SUCCESS;
}
