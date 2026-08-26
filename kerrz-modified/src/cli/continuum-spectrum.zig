const std = @import("std");
const clippy = @import("clippy");
const kerrz = @import("kerrz");

const utils = @import("./utils.zig");

pub const short_description = "Calculate continuum profiles.";
pub const description = short_description ++ "\n";

pub const Args = clippy.Arguments([_]clippy.ArgumentDescriptor{
    .{
        .arg = "--spin spin",
        .argtype = f64,
        .default = "0.998",
        .help = "The black hole spin.",
    },
    .{
        .arg = "--incl inclination",
        .argtype = f64,
        .default = "80",
        .help = "The observer inclination in degrees.",
    },
    .{
        .arg = "--radius r",
        .argtype = f64,
        .default = "1e7",
        .help = "The observer radial distance in rg.",
    },
    .{
        .arg = "--verbose",
        .help = "Print additional diagnostic information.",
    },
    .{
        .arg = "-o/--output path",
        .default = "continuum.fits",
        .help = "The path to write the transfer function to.",
    },
    .{
        .arg = "--nangles n",
        .argtype = usize,
        .default = "200",
        .help = "The number of points along the ring to trace. This does not apply to the lamppost geometry.",
    },
    utils.VelocityProfile.ArgDescriptor,
} ++ utils.CoronalParameters.ArgDescriptorList());

const Dual = kerrz.DualNumber(f64, 0);

pub fn run(
    out: *std.Io.Writer,
    allocator: std.mem.Allocator,
    itt: *clippy.ArgumentIterator,
) !void {
    const args = try Args.initParseAll(itt, .{});

    const metric: kerrz.KerrMetric(Dual) = .init(.one, .promote(args.spin));
    const coronal_parameters = try utils.CoronalParameters.parse(args);
    const v_prof = try utils.VelocityProfile.fromArgs(args);
    const x_obs: kerrz.FourVector(Dual) = .{
        .t = .zero,
        .r = .promote(args.radius),
        .th = .promote(std.math.degreesToRadians(args.incl)),
        .ph = .zero,
    };

    switch (coronal_parameters) {
        .lamppost => |p| {
            try out.print("Lamppost model h:{d:.5}, vr:{d:.5}\n", .{
                p.height,
                p.radial_velocity,
            });
        },
        .ring => |p| {
            try out.print("Ring model h:{d:.5}, x:{d:.5} '{s}' using {d} angular points\n", .{
                p.height,
                p.radius,
                @tagName(v_prof),
                args.nangles,
            });
        },
        // TODO: implement me
        .disc, .umbrella => unreachable,
    }

    // TODO: combine this with the `emissivity` command
    const corona_opts: kerrz.CoronalModelOptions(Dual) = switch (coronal_parameters) {
        .lamppost => |p| .{
            .lamppost = .{
                .height = .promote(p.height),
                .radial_velocity = .promote(p.radial_velocity),
            },
        },
        .ring => |p| .{
            .ring = .{
                .height = .promote(p.height),
                .radius = .promote(p.radius),
                .velocity = v_prof,
            },
        },
        // TODO: implement me
        .disc, .umbrella => unreachable,
    };

    const corona = try corona_opts.toCoronalModel(allocator, metric);
    defer corona.deinit(allocator);

    switch (corona) {
        .lamppost => try calculateAndPrintCentroidTime(
            out,
            metric,
            x_obs,
            corona,
            args.verbose,
        ),
        .ring => try calculateAndPrintTransferFunction(
            out,
            allocator,
            metric,
            x_obs,
            corona_opts,
            corona,
            args,
        ),
        // TODO: implement me
        .umbrella, .disc => unreachable,
    }
}

fn calculateAndPrintCentroidTime(
    out: *std.Io.Writer,
    metric: kerrz.KerrMetric(Dual),
    x_obs: kerrz.FourVector(Dual),
    corona: kerrz.CoronalModel(Dual),
    verbose: bool,
) !void {
    const sol = try kerrz.continuum.continuumGeodesic(
        Dual,
        metric,
        x_obs,
        corona,
    );

    const g = sol.calculateEnergyshift(metric, corona);

    try out.writeAll("\n");
    try out.print("          x_obs.r : {d:.6}\n", .{x_obs.r.x});
    try out.print(
        "         x_obs.th : {d:.6}°\n",
        .{std.math.radiansToDegrees(x_obs.th.x)},
    );
    try out.print("     Time - r_obs : {d:.6}\n", .{sol.delta_t.x});
    try out.print("        Mino time : {d:.6}\n", .{sol.result.mino_time.x});
    try out.print("      Energyshift : {d:.6}\n", .{g.x});
    try out.print("        |∂θ / ∂Y| : {d:.6}\n", .{sol.jac});
    try out.print("  |∂cosθ / ∂cosY| : {d:.6}\n", .{sol.jacobianAsCosine()});
    if (verbose) {
        try out.print(
            "                Y : {d:.6}°\n",
            .{std.math.radiansToDegrees(sol.local_theta.x)},
        );
        try out.print("                α : {d:.6}\n", .{sol.alpha.x});
        try out.print("                β : {d:.6}\n", .{sol.beta.x});
        try out.print("      Error (abs) : {e:.3}\n", .{sol.err});
        try out.print("          f calls : {d}\n", .{sol.f_calls});
        try out.print("     Angular case : {s}\n", .{@tagName(sol.result.state.angular_case)});
        try out.print("      Radial case : {s}\n", .{@tagName(sol.result.state.radial_case)});
    }
    try out.writeAll("\n");
}

fn calculateAndPrintTransferFunction(
    out: *std.Io.Writer,
    allocator: std.mem.Allocator,
    metric: kerrz.KerrMetric(Dual),
    x_obs: kerrz.FourVector(Dual),
    corona_options: kerrz.CoronalModelOptions(Dual),
    corona: kerrz.CoronalModel(Dual),
    args: Args.Parsed,
) !void {
    const transfer = try kerrz.continuum.transferFunction(
        Dual,
        allocator,
        metric,
        x_obs,
        corona,
        .{
            .max_points = args.nangles,
        },
    );
    defer transfer.deinit(allocator);

    const delta_angle = 1.0 / @as(Dual.T, @floatFromInt(transfer.traces.len));

    // Calculate some summary statistics:
    var g_min: Dual.T = std.math.floatMax(Dual.T);
    var g_max: Dual.T = std.math.floatMin(Dual.T);
    var g_avg: Dual.T = std.math.floatMin(Dual.T);
    var delta_t_min: Dual.T = std.math.floatMax(Dual.T);
    var delta_t_max: Dual.T = std.math.floatMin(Dual.T);
    var delta_t_avg: Dual.T = std.math.floatMin(Dual.T);
    var jac_min: Dual.T = std.math.floatMax(Dual.T);
    var jac_max: Dual.T = std.math.floatMin(Dual.T);
    var jac_avg: Dual.T = std.math.floatMin(Dual.T);
    var Y_min: Dual.T = std.math.floatMax(Dual.T);
    var Y_max: Dual.T = std.math.floatMin(Dual.T);
    var Y_avg: Dual.T = std.math.floatMin(Dual.T);
    var psi_min: Dual.T = std.math.floatMax(Dual.T);
    var psi_max: Dual.T = std.math.floatMin(Dual.T);
    var psi_avg: Dual.T = std.math.floatMin(Dual.T);
    var alpha_min: Dual.T = std.math.floatMax(Dual.T);
    var alpha_max: Dual.T = std.math.floatMin(Dual.T);
    var alpha_avg: Dual.T = std.math.floatMin(Dual.T);
    var beta_min: Dual.T = std.math.floatMax(Dual.T);
    var beta_max: Dual.T = std.math.floatMin(Dual.T);
    var beta_avg: Dual.T = std.math.floatMin(Dual.T);

    var f_calls_total: usize = 0;
    for (transfer.traces) |trace| {
        g_min = @min(trace.g.x, g_min);
        g_max = @max(trace.g.x, g_max);
        delta_t_min = @min(trace.delta_t.x, delta_t_min);
        delta_t_max = @max(trace.delta_t.x, delta_t_max);
        jac_min = @min(trace.jac, jac_min);
        jac_max = @max(trace.jac, jac_max);
        Y_min = @min(trace.local_theta.x, Y_min);
        Y_max = @max(trace.local_theta.x, Y_max);
        psi_min = @min(trace.local_phi.x, psi_min);
        psi_max = @max(trace.local_phi.x, psi_max);
        alpha_min = @min(trace.alpha.x, alpha_min);
        alpha_max = @max(trace.alpha.x, alpha_max);
        beta_min = @min(trace.beta.x, beta_min);
        beta_max = @max(trace.beta.x, beta_max);

        g_avg += trace.g.x * delta_angle;
        delta_t_avg += trace.delta_t.x * delta_angle;
        jac_avg += trace.jac * delta_angle;
        Y_avg += trace.local_theta.x * delta_angle;
        psi_avg += trace.local_phi.x * delta_angle;
        alpha_avg += trace.alpha.x * delta_angle;
        beta_avg += trace.beta.x * delta_angle;

        f_calls_total += trace.f_calls;
    }

    try out.writeAll("\n");
    try out.print("          x_obs.r : {d:.6}\n", .{x_obs.r.x});
    try out.print(
        "         x_obs.th : {d:.6}°\n",
        .{std.math.radiansToDegrees(x_obs.th.x)},
    );
    try out.print("     Time - r_obs : {d:.6} [{d:.6} - {d:.6}]\n", .{ delta_t_avg, delta_t_min, delta_t_max });
    try out.print("      Energyshift : {d:.6} [{d:.6} - {d:.6}]\n", .{ g_avg, g_min, g_max });
    try out.print("        |∂θ / ∂Y| : {d:.6} [{d:.6} - {d:.6}]\n", .{ jac_avg, jac_min, jac_max });
    if (args.verbose) {
        try out.print("          f calls : {d}\n", .{f_calls_total});
        try out.print("                Y : {d:.6}° [{d:.6}° - {d:.6}°]\n", .{
            std.math.radiansToDegrees(Y_avg),
            std.math.radiansToDegrees(Y_min),
            std.math.radiansToDegrees(Y_max),
        });
        try out.print("                Ψ : {d:.6}° [{d:.6}° - {d:.6}°]\n", .{
            std.math.radiansToDegrees(psi_avg),
            std.math.radiansToDegrees(psi_min),
            std.math.radiansToDegrees(psi_max),
        });
        try out.print("                α : {d:.6} [{d:.6} - {d:.6}]\n", .{ alpha_avg, alpha_min, alpha_max });
        try out.print("                β : {d:.6} [{d:.6} - {d:.6}]\n", .{ beta_avg, beta_min, beta_max });
    }
    try out.writeAll("\n");

    const fits = try kerrz.FitsFile.init(allocator);
    defer fits.deinit();

    const hdu_ptr = try fits.addHdu();
    // Serialise to a FITS file.
    hdu_ptr.* = try transfer.toFITS(allocator);
    try corona_options.addToHdu(hdu_ptr);
    try metric.addToHdu(hdu_ptr);

    try out.print("Writing continuum transfer data to '{s}'\n", .{args.output});
    try fits.save(.{ .path = args.output });
}
