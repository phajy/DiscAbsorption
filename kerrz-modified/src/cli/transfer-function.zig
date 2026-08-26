const std = @import("std");
const clippy = @import("clippy");
const kerrz = @import("kerrz");

const utils = @import("utils.zig");

pub const short_description = "Calculate transfer functions.";
pub const description =
    \\Calculate various relativistic transfer functions that can be used as tables or
    \\integration kernels for relativistic effects. This command only calculates
    \\individual transfer functions, and is for exploring different methods and
    \\heuristics for their calculation. For tables and batch operations, use other
    \\commands, such as `table`.
    \\
    \\Cunningham's transfer function (Cunningham 1975) is a particularly useful
    \\formulation of a numerically stable transfer function. They can be integrated
    \\to transfer any quantity that can be expressed in coordinates of the disc back
    \\to the observer. The parameterisation is in terms of the radius on the disc and
    \\the redshift along that annulus.
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
        .arg = "--incl inclination",
        .argtype = f64,
        .default = "80",
        .help = "The observer inclination in degrees.",
    },
    .{
        .arg = "-r/--target-radius r",
        .argtype = f64,
        .default = "5.0",
        .help = "The target radius on the disc to calculate Cunningham's transfer function for.",
    },
    .{
        .arg = "-o/--output filename",
        .default = "ctf.fits",
        .help = "The path to write the output transfer function to.",
    },
} ++ &CTFArguments);

/// General Cunningham transfer function arguments.
pub const CTFArguments = [_]clippy.ArgumentDescriptor{
    .{
        .arg = "--nangles n",
        .argtype = usize,
        .default = "200",
        .help = "The number of points (angles on the image plane) to calculate the transfer function for.",
    },
    .{
        .arg = "--min-guess value",
        .default = "5",
        .argtype = f64,
        .help = "The minimum initial offset value that the solver should use. A larger value can make the solver more stable at the cost of performance. Note that different heuristics may fail to converge for different values!",
    },
    .{
        .arg = "--initial-guess value",
        .argtype = f64,
        .help = "The initial guess that the solver should try when root finding a radius on the disc. If not supplied, uses the target radius as the initial guess.",
    },
    .{
        .arg = "--heuristic f",
        .help = "If set, use a heuristic function to refine the angles traces for the transfer function. The possible values are `arclength` or `impact`.",
    },
    .{
        .arg = "--refine ref",
        .display_name = "--refine N,M",
        .help = "Refine `M` times around the last `N` points about the extremal energyshift. This can help when the transfer function looks ill-determined because of a poor estimate of the extremal values. If not set, no refine step will occur.",
    },
    .{
        .arg = "--optimise n",
        .argtype = usize,
        .default = "17",
        .help = "How many evaluations should be used for optimsing the extremal energyshift values. This is seperate from `refine`, as it uses a Golden-section search to try to extremise the energyshift values, whereas refine is simple a refinement.",
    },
};

const Dual = kerrz.DualNumber(f64, 0);
const TF = kerrz.tools.TransferFunction(Dual);
const Table = kerrz.transfer_tables.CunninghamTransferFunctionTable(Dual);

pub fn run(out: *std.Io.Writer, allocator: std.mem.Allocator, itt: *clippy.ArgumentIterator) !void {
    const args = try Args.initParseAll(itt, .{});

    const metric: kerrz.KerrMetric(Dual) = .init(.one, .promote(args.spin));
    const x_obs: kerrz.FourVector(Dual) = .{
        .t = .zero,
        .r = .promote(1e7),
        .th = .promote(std.math.degreesToRadians(args.incl)),
        .ph = .zero,
    };

    const tool: TF = .{
        .metric = metric,
        .x_obs = x_obs,
        .r_target = args.@"target-radius",
        .opts = .{
            .minimum_guess = args.@"min-guess",
            .initial_guess = args.@"initial-guess",
            .heuristic = try .fromOptionalString(args.heuristic),
            .max_points = args.nangles,
            .refine = try .parse(args.refine),
            .optimise = args.optimise,
        },
    };

    const time_now = std.time.milliTimestamp();

    var result = try tool.run(allocator);
    defer result.deinit(allocator);

    const duration = std.time.milliTimestamp() - time_now;
    try out.print("Time elapsed: {D}\n", .{duration * std.time.ns_per_ms});

    const ctf = result.ctf;

    var largest_error: f64 = 0;
    for (ctf.traces) |trace| {
        if (trace.r_err > largest_error) {
            largest_error = trace.r_err;
        }
    }

    try out.print("Target radius     = {d}\n", .{ctf.target_radius});
    try out.print(
        "Traces[0] angle   = {d:.3}°\n",
        .{std.math.radiansToDegrees(ctf.traces[0].image_angle)},
    );
    try out.print(
        "Traces[end] angle = {d:.3}°\n",
        .{std.math.radiansToDegrees(ctf.traces[ctf.traces.len - 1].image_angle)},
    );
    try out.print("g_min             = {d:.6}\n", .{ctf.g_min});
    try out.print("g_max             = {d:.6}\n", .{ctf.g_max});
    try out.print("g_min index       = {d}\n", .{ctf.g_min_index});
    try out.print("g_max index       = {d}\n", .{ctf.g_max_index});
    try out.print(
        "g_min angle       = {d:.3}°\n",
        .{std.math.radiansToDegrees(ctf.traces[ctf.g_min_index].image_angle)},
    );
    try out.print(
        "g_max angle       = {d:.3}°\n",
        .{std.math.radiansToDegrees(ctf.traces[ctf.g_max_index].image_angle)},
    );
    try out.print("Largest error     = {e:.4}\n", .{largest_error});

    const save_path = args.output;
    try out.print("Saving {d} points to file: '{s}'.\n", .{ ctf.traces.len, save_path });

    var tmp_ctf: [1]Table.CTF = .{ctf};

    var table: Table = .{
        .tracer = .init(metric, x_obs, .equatorial_plane),
        .transfer_functions = .fromOwnedSlice(&tmp_ctf),
    };

    const fits = try kerrz.FitsFile.init(allocator);
    defer fits.deinit();
    (try fits.addHdu()).* = try table.toFITS(allocator);

    try fits.save(.{ .path = save_path });
}
