const std = @import("std");
const clippy = @import("clippy");
const kerrz = @import("kerrz");

const utils = @import("./utils.zig");

const MAX_MAPPERS = 64;

const Dual = kerrz.DualNumber(f64, 2);

pub const short_description = "Trace a single geodesic and print information.";
pub const description =
    \\Trace a single geodesic specified by impact parameters on an image plane.
    \\This command can then either serialise information about the geodesic, or
    \\print diagnostic information about the solver, such as the
    \\anti-derivatives or potentials.
    \\
    \\If no `--to-angle` or `--disc` argument is given, the default behaviour
    \\is the equivalent of `--to-angle 90`.
    \\
;

pub const Args = clippy.Arguments([_]clippy.ArgumentDescriptor{
    .{
        .arg = "--spin spin",
        .argtype = f64,
        .default = "0.998",
        .help = "The black hole spin.",
    },
    .{
        .arg = "-t/--to-angle angle",
        .argtype = f64,
        .help = "The opening angle of the cone to trace the geodesic to, in degrees. A value of 90 is the equatorial plane.",
    },
    utils.Position.ArgDescriptorDefaults(1e6, 45, 0),
    .{
        .arg = "--winding winding",
        .argtype = usize,
        .default = "0",
        .help = "How many windings around the equatorial plane to trace. Setting this paramter to 0 is the same as `--no-false-image` in the `image` command.",
    },
    .{
        // TODO: make this support the usual `--velocity` argument
        .arg = "--co-rotate",
        .help = "Set the origin to be co-rotating with a Keplerian orbit at the equivalent radius.",
    },
    .{
        .arg = "--verbose",
        .help = "Print more information about the geodesic than usual.",
    },
    .{
        .arg = "-o/--output filename",
        .help = "Save the full trajectory to a FITS file.",
    },
    .{
        .arg = "--num-points num-points",
        .argtype = usize,
        .default = "1024",
        .help = "For saving: the number of points to save. This is equivalently the number of points to linearly parameterise the Mino time interval into.",
    },
    utils.MapperArg(""),
    .{
        .arg = "--min-time min",
        .argtype = f64,
        .default = "0",
        .help = "The minimum (Mino) time to calculate.",
    },
    .{
        .arg = "--interpolate n",
        .argtype = usize,
        .help = "How many points to interpolate between each pair of calculated points. This may be used to upscale the trajectory. It uses a cubic-spline interpolation over the trajectory. This flag is currently only used when saving the traced trajectory.",
    },
    .{
        .arg = "--max-time start",
        .argtype = f64,
        .help = "When saving the full path to file, this flag controls the maximum (Mino) time that will be calcualted. If not set, this is determined by tracing to `--angle`.",
    },
} ++ utils.DiscParameters.ArgDescriptorList() ++ utils.InitialParameters.ArgDescriptorList());

const DefaultParams: utils.InitialParameters = .{ .defaults = .{
    .impact_parameters = .{ .alpha = 1.0, .beta = 1.0 },
} };

pub fn run(out: *std.Io.Writer, allocator: std.mem.Allocator, itt: *clippy.ArgumentIterator) !void {
    const args = try Args.initParseAll(itt, .{});

    const mappers = try utils.parseMapperValues(allocator, args.map);
    defer allocator.free(mappers);

    if (mappers.len > MAX_MAPPERS) {
        return error.TooManyMappers;
    }

    const metric: kerrz.KerrMetric(Dual) = .init(.one, .promote(args.spin));

    const pos = try utils.Position.fromString(args.position);

    const x: kerrz.FourVector(Dual) = .{
        .t = .zero,
        .r = .promote(pos.r),
        .th = .promote(pos.theta),
        .ph = .promote(pos.phi),
    };

    const params = try DefaultParams.parse(args);

    if (args.@"co-rotate" and params == .impact_parameters) {
        return error.IncompatibleArguments;
    }

    const ts = metric.tangentSpace(x);

    const v_source = if (args.@"co-rotate")
        kerrz.orbits.coRotating(Dual, metric, ts)
    else
        kerrz.orbits.stationary(Dual, ts);

    const geod = switch (params) {
        .impact_parameters => |p| kerrz.NullGeodesic(Dual).fromImpactParameters(
            metric,
            x,
            Dual.promote(p.alpha).diff(0),
            Dual.promote(p.beta).diff(1),
        ),
        .sky_angles => |p| kerrz.NullGeodesic(Dual).fromSkyAnglesTangentSpace(
            metric,
            ts,
            v_source,
            Dual.promote(p.theta).diff(0),
            Dual.promote(p.phi).diff(1),
        ),
    };

    const result = try traceSingleGeodesic(args, metric, geod);

    const total = result.totalAntiderivatives(metric, geod);

    const eshift = kerrz.redshift.keplerianRedshiftResult(
        Dual,
        metric,
        geod,
        result,
        v_source,
    );

    const ctx: TraceContext = .{
        .geometry = metric,
        .v_source = v_source,
    };

    try printInfo(
        out,
        ctx,
        geod,
        result,
        total,
        eshift,
        args.verbose,
        params,
        v_source,
        mappers,
    );

    const max_time = args.@"max-time" orelse result.mino_time.x;

    if (args.output) |output_path| {
        const min_time = args.@"min-time";
        const num_points = args.@"num-points";

        const total_mino_interval = max_time - min_time;
        const dt = total_mino_interval / @as(Dual.T, @floatFromInt(num_points));

        var path = std.ArrayList(PathPoint).empty;
        defer path.deinit(allocator);

        const path_builder = geod.traceBuilder(metric, .{});

        for (0..num_points) |i| {
            const t = @min(
                max_time,
                @as(Dual.T, @floatFromInt(i)) * dt + min_time,
            );

            const p = path_builder.atMinoTime(.promote(t));

            if (p.r.x < metric.horizon_radius.x) {
                break;
            }

            const p_total = p.totalAntiderivatives(metric, geod);

            const point = try path.addOne(allocator);
            point.* = .{
                .mino_time = t,
                .t = p_total.coordinateTime(metric, geod).x,
                .r = p.r.x,
                .theta = p.theta.x,
                .phi = p_total.coordinateAzimuth(metric, geod).x,
            };

            for (mappers) |m| {
                point.extra_name[point.extra_count] = @tagName(m.mapper.value);
                point.extra[point.extra_count] = m.mapper.calculate(
                    Dual,
                    .{
                        .em = ctx.em,
                        .geometry = ctx.geometry,
                        .v_source = ctx.v_source,
                    },
                    geod,
                    p,
                );
                point.extra_count += 1;
            }
        }

        if (args.interpolate) |n_interpolate| {
            try interpolateExtraPoints(allocator, &path, n_interpolate);
        }

        try out.print(
            "Saving {d} points to file: '{s}'.\n",
            .{ path.items.len, output_path },
        );

        const fits = try kerrz.FitsFile.init(allocator);
        defer fits.deinit();

        const hdu = try fits.newHdu("GEODPATH", .binary_table);
        try kerrz.addKerrzFITSInfo(hdu);
        try metric.addToHdu(hdu);

        const table = &hdu.data.binary_table;
        try table.appendColumn(.{
            .label = "mino",
            .comment = "The Mino time at this point along the geodesic",
        });
        try table.appendColumn(.{
            .label = "t",
            .comment = "The coordinate time",
            .units = "tg",
            .units_comment = "In GM/c^3 units",
        });
        try table.appendColumn(.{
            .label = "r",
            .comment = "The radial Boyer-Lindquist coordinate",
            .units = "rg",
            .units_comment = "In GM/c^2 units",
        });
        try table.appendColumn(.{
            .label = "theta",
            .comment = "The polar coordinate angle",
            .units = "rad",
        });
        try table.appendColumn(.{
            .label = "phi",
            .comment = "The azimuthal coordinate angle",
            .units = "rad",
        });

        for (mappers) |mapper| {
            try table.appendColumn(.{
                .label = mapper.name,
                .comment = "Extra field. See kerrz --help mapper",
            });
        }

        for (path.items) |p| {
            const row = try table.addRow();
            row.cols[0].one.float_32 = @floatCast(p.mino_time);
            row.cols[1].one.float_32 = @floatCast(p.t);
            row.cols[2].one.float_32 = @floatCast(p.r);
            row.cols[3].one.float_32 = @floatCast(p.theta);
            row.cols[4].one.float_32 = @floatCast(p.phi);

            for (row.cols[5..], p.extra[0..p.extra_count]) |*col, extra| {
                col.one.float_32 = @floatCast(extra);
            }
        }

        try fits.save(.{ .path = output_path });
    }
}

fn traceSingleGeodesic(
    args: Args.Parsed,
    metric: kerrz.KerrMetric(Dual),
    geod: kerrz.NullGeodesic(Dual),
) !kerrz.TraceResult(Dual) {
    // Make sure no incompatible arguments are given
    var active: usize = 0;

    // Parse the disc parameters
    const disc = try utils.DiscParameters.parse(args);

    if (disc != .no_disc) active += 1;
    if (args.@"to-angle" != null) active += 1;
    if (args.@"max-time" != null) active += 1;

    if (active > 1) {
        return error.IncompatibleArguments;
    }

    if (disc != .no_disc) {
        return geod.traceDisc(
            metric,
            disc.toAccretionDisc(Dual, metric),
            .{ .winding = args.winding },
        );
    }

    if (args.@"max-time") |max_time| {
        const builder = geod.traceBuilder(metric, .{});
        return builder.atMinoTime(.promote(max_time));
    }

    // fallback
    const angle = args.@"to-angle" orelse 90;
    const target_angle: Dual = .promote(
        std.math.degreesToRadians(angle),
    );
    return geod.traceToAngle(
        metric,
        target_angle,
        .{ .winding = args.winding },
    );
}

const TraceContext = struct {
    geometry: kerrz.KerrMetric(Dual),
    em: kerrz.emissivity.EmissivityProfile(Dual) = .{ .powerlaw = .{} },
    v_source: kerrz.FourVector(Dual),
};

const PathPoint = struct {
    mino_time: Dual.T,
    t: Dual.T,
    r: Dual.T,
    theta: Dual.T,
    phi: Dual.T,

    /// Optional extra fields
    extra_name: [MAX_MAPPERS][]const u8 = undefined,
    extra: [MAX_MAPPERS]Dual.T = undefined,
    extra_count: usize = 0,
};

fn printInfo(
    out: *std.Io.Writer,
    ctx: TraceContext,
    geod: kerrz.NullGeodesic(Dual),
    result: kerrz.TraceResult(Dual),
    total: kerrz.antiderivatives.TotalAntiderivatives(Dual),
    eshift: Dual,
    verbose: bool,
    params: utils.InitialParameters.Parameters,
    v_source: kerrz.FourVector(Dual),
    mappers: []const utils.NamedMapper,
) !void {
    const x_sym = if (params == .impact_parameters) "α" else "θ";
    const y_sym = if (params == .impact_parameters) "β" else "φ";

    const x_val = if (params == .impact_parameters)
        params.impact_parameters.alpha
    else
        std.math.radiansToDegrees(params.sky_angles.theta);
    const y_val = if (params == .impact_parameters)
        params.impact_parameters.beta
    else
        std.math.radiansToDegrees(params.sky_angles.phi);

    const x_unit = if (params == .impact_parameters) "" else "°";

    try out.writeAll("Initial conditions:\n");
    try out.print("r_obs             = {d}\n", .{geod.x_init.r.x});
    try out.print("Θ_obs             = {d}°\n", .{std.math.radiansToDegrees(geod.x_init.th.x)});
    try out.print("φ_obs             = {d}°\n", .{std.math.radiansToDegrees(geod.x_init.ph.x)});
    try out.print("{s}_sky             = {d}{s}\n", .{ x_sym, x_val, x_unit });
    try out.print("{s}_sky             = {d}{s}\n", .{ y_sym, y_val, x_unit });

    if (params != .impact_parameters) {
        try out.print("v_t               = {d}\n", .{v_source.t.x});
        try out.print("v_phi             = {d}\n", .{v_source.ph.x});
    }

    try out.writeAll("\nGeneral values:\n");
    try out.print("Mino time (angle) = {d}\n", .{result.mino_time.x});
    try out.print("Status            = {s}\n", .{@tagName(result.status)});
    try out.print("Case              = {s}\n", .{@tagName(result.state.radial_case)});
    try out.print("Angular Case      = {s}\n", .{@tagName(result.state.angular_case)});
    try out.print("Winding           = {d}\n", .{result.winding});
    try out.print("η                 = {d}\n", .{geod.eta.x});
    try out.print("λ                 = {d}\n", .{geod.lambda.x});
    try out.print("t                 = {d}\n", .{total.coordinateTime(ctx.geometry, geod).x});
    try out.print("r                 = {d}\n", .{result.r.x});
    try out.print("Θ                 = {d}\n", .{result.theta.x});
    try out.print("φ                 = {d}\n", .{total.coordinateAzimuth(ctx.geometry, geod).x});
    try out.print("g (E/E0)          = {d}\n", .{eshift.x});

    if (mappers.len > 0) {
        try out.writeAll("\nMapper values:\n");
        for (mappers) |m| {
            try out.print(
                "{s: <18} = {d}\n",
                .{ m.name, m.mapper.calculate(
                    Dual,
                    .{
                        .em = ctx.em,
                        .geometry = ctx.geometry,
                        .v_source = ctx.v_source,
                    },
                    geod,
                    result,
                ) },
            );
        }
    }

    try out.writeAll("\nDerivatives:\n");
    try out.print("dr/d{s}             = {d}\n", .{ x_sym, result.r.dx[0] });
    try out.print("dr/d{s}             = {d}\n", .{ x_sym, result.r.dx[1] });
    try out.print("dg/d{s}             = {d}\n", .{ y_sym, eshift.dx[0] });
    try out.print("dg/d{s}             = {d}\n", .{ y_sym, eshift.dx[1] });

    if (verbose) {
        try out.writeAll("\nRadial roots:\n");
        try out.print(
            "r1                = {d} + {d}i\n",
            .{ result.state.radial_roots.r1.real().x, result.state.radial_roots.r1.imag().x },
        );
        try out.print(
            "r2                = {d} + {d}i\n",
            .{ result.state.radial_roots.r2.real().x, result.state.radial_roots.r2.imag().x },
        );
        try out.print(
            "r3                = {d} + {d}i\n",
            .{ result.state.radial_roots.r3.real().x, result.state.radial_roots.r3.imag().x },
        );
        try out.print(
            "r4                = {d} + {d}i\n",
            .{ result.state.radial_roots.r4.real().x, result.state.radial_roots.r4.imag().x },
        );
        try out.writeAll("\nAngular roots:\n");
        try out.print("u plus            = {d}\n", .{result.state.angular_roots.u_plus.x});
        try out.print("u minus           = {d}\n", .{result.state.angular_roots.u_minus.x});
        // G&L, Equation (20) to (23)
        try out.print(
            "θ1                = {d}\n",
            .{std.math.acos(@sqrt(result.state.angular_roots.u_plus.x))},
        );
        try out.print(
            "θ2                = {d}\n",
            .{std.math.acos(@sqrt(result.state.angular_roots.u_minus.x))},
        );
        try out.print(
            "θ3                = {d}\n",
            .{std.math.acos(-@sqrt(result.state.angular_roots.u_minus.x))},
        );
        try out.print(
            "θ4                = {d}\n",
            .{std.math.acos(-@sqrt(result.state.angular_roots.u_plus.x))},
        );

        const turn_times = result.state.minoTimeToAngularTurns(
            geod.windings,
            @floatFromInt(geod.theta_sign),
        );
        try out.print(
            "θ1 Mino time      = {d}\n",
            .{turn_times.tau_0.x},
        );
        try out.print(
            "θ4 Mino time      = {d}\n",
            .{turn_times.tau_1.x},
        );
        try out.print(
            "Half period       = {d}\n",
            .{result.state.angular_cache.G_theta_half_libration.x},
        );

        try out.writeAll("\nAntiderivatives:\n");
        try out.print("I_0               = {d}\n", .{total.radial.I_0.x});
        try out.print("I_1               = {d}\n", .{total.radial.I_1.x});
        try out.print("I_2               = {d}\n", .{total.radial.I_2.x});
        try out.print("I_plus            = {d}\n", .{total.radial.I_plus.x});
        try out.print("I_minus           = {d}\n", .{total.radial.I_minus.x});
        try out.print("G_φ (init)        = {d}\n", .{total.angular.G_phi_init.x});
        try out.print("G_φ (final)       = {d}\n", .{total.angular.G_phi_final.x});
        try out.print("G_t (init)        = {d}\n", .{total.angular.G_t_init.x});
        try out.print("G_t (final)       = {d}\n", .{total.angular.G_t_final.x});
        try out.print("G_φ (halflib)     = {d}\n", .{total.angular.G_phi_half.x});
        try out.print("G_t (halflib)     = {d}\n", .{total.angular.G_t_half.x});
        try out.print("G_θ (init)        = {d}\n", .{result.state.angular_cache.G_theta_init.x});
        try out.print("G_θ (halflib)     = {d}\n", .{result.state.angular_cache.G_theta_half_libration.x});
        const angular_values = kerrz.antiderivatives.angularFromCache(
            Dual,
            result.state.angular_case,
            result.state.angular_cache,
            result.theta,
        );
        try out.print("G_θ (final)       = {d}\n", .{angular_values.G_theta_final.x});
    }
}

fn interpolateExtraPoints(
    allocator: std.mem.Allocator,
    path: *std.ArrayList(PathPoint),
    n_interpolate: usize,
) !void {
    const RangeIterator = kerrz.iterators.RangeIterator(Dual.T);
    // Unpack into the required format for the interpolation
    var pathpoints: kerrz.geodesic.PathPoints(Dual) = .empty;
    defer pathpoints.deinit(allocator);

    for (path.items) |p| {
        try pathpoints.addPoint(
            allocator,
            p.mino_time,
            p.r,
            p.theta,
        );
    }

    const interp = try pathpoints.interpolator(allocator);
    defer interp.deinit(allocator);

    const N = path.items.len;
    for (0..N - 1) |i| {
        const p = path.items[i];
        var itt = RangeIterator.init(
            p.mino_time,
            path.items[i + 1].mino_time,
            n_interpolate + 2,
        );

        // Ignore the first and last point.
        _ = itt.next();
        itt.remaining -|= 1;

        while (itt.next()) |mino| {
            const new = interp.interpolate(mino);
            try path.append(allocator, .{
                .extra = p.extra,
                .extra_count = p.extra_count,
                .extra_name = p.extra_name,
                .mino_time = new.mino,
                .phi = 0,
                .r = new.r,
                .t = 0,
                .theta = new.theta,
            });
        }
    }
}
