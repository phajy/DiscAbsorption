const std = @import("std");
const clippy = @import("clippy");
const kerrz = @import("kerrz");

const TableArgs = @import("tables.zig").TableArgs;
const transfer_function = @import("transfer-function.zig");

const utils = @import("utils.zig");

pub const short_description = "Calculate line profiles.";
pub const description =
    \\Calculates line profiles using Cunningham's transfer function.
    \\
    \\This command outputs a single table with two columns: the first is the
    \\energyshift grid, and the second is the observed flux.
    \\
;

pub const Args = clippy.Arguments(&[_]clippy.ArgumentDescriptor{
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
        .arg = "--nthreads nthreads",
        .argtype = usize,
        .help = "The number of CPU threads to use.",
    },
    .{
        .arg = "-o/--output filename",
        .default = "lineprofile.dat",
        .help = "The name of the file to write the lineprofile to.",
    },
    .{
        .arg = "--convolve filename",
        .help = "A disc-spectrum file to convolve as a function of radius on the disc.",
    },
    .{
        .arg = "--no-norm",
        .help = "Whether to normalise the lineprofile so that it integrates to unity. If this is flag is passed, no normalisation is performed and the units of the lineprofile are directly in flux.",
    },
    .{
        .arg = "--using-table path",
        .help = "Instead of caluclating the transfer function table, read and parse it from a file given by `path`. Note this has a few limitations, and will silently ignore other arguments, such as the spin or table arguments.",
    },
    utils.EmissivityFunction.ArgDescriptor,
} ++ utils.LineprofileIntegrationArguments ++ TableArgs);

const Dual = kerrz.DualNumber(f64, 0);
const Tool = kerrz.tools.Lineprofile(Dual);
const Table = kerrz.tools.TransferFunctionTable(Dual);

pub fn run(out: *std.Io.Writer, allocator: std.mem.Allocator, itt: *clippy.ArgumentIterator) !void {
    const args = try Args.initParseAll(itt, .{});

    var emissivity = try utils.EmissivityFunction.fromArgs(Dual, allocator, args);
    defer emissivity.deinit(allocator);

    const num_threads = utils.getNumThreads(args.nthreads);

    var threads = try kerrz.ThreadMap.init(
        allocator,
        .{ .num_threads = num_threads },
    );
    defer threads.deinit();

    var table = b: {
        if (args.@"using-table") |table_path| {
            break :b try utils.readOrReadAndInterpolateTable(
                Dual,
                out,
                allocator,
                table_path,
                args.spin,
                args.incl,
            );
        } else {
            const metric: kerrz.KerrMetric(Dual) = .init(.one, .promote(args.spin));
            const x_obs: kerrz.FourVector(Dual) = .{
                .t = .zero,
                .r = .promote(1e7),
                .th = .promote(std.math.degreesToRadians(args.incl)),
                .ph = .zero,
            };
            var table_tool: Table = .{
                .metric = metric,
                .x_obs = x_obs,
                .table_opts = .{
                    .num_radii = args.nradii,
                    .r_in = args.rin,
                    .r_out = args.rout,
                    .r_grid = try .fromString(args.rgrid),
                },
                .tf_opts = .{
                    .minimum_guess = args.@"min-guess",
                    .initial_guess = args.@"initial-guess",
                    .heuristic = try .fromOptionalString(args.heuristic),
                    .max_points = args.nangles,
                    .refine = try .parse(args.refine),
                    .optimise = args.optimise,
                },
            };
            break :b try table_tool.run(
                allocator,
                threads,
                .{ .out = out, .thread_chunk_size = 1 },
            );
        }
    };

    defer table.deinit(allocator);

    if (!table.table.hasField(.f)) {
        return itt.throwError(
            error.MissingField,
            "Transfer function is missing `f` field\n",
            .{},
        );
    }

    const tool: Tool = .{
        .lp_opts = .{
            .emissivity = emissivity.profile,
            .num_r_steps = args.nrsteps,
            .r_step_grid = try .fromString(args.rstepgrid),
            .r_min = args.rin orelse 0,
            .r_max = args.rout,
            .normalise = !args.@"no-norm",
        },
        .num_g = args.ng,
        .num_fine_g = args.ngstar,
        .table = table.table,
    };

    var result = try tool.run(allocator, .{ .out = out });
    defer result.deinit(allocator);

    try out.print(
        "Integration time on single thread: {D}\n",
        .{result.integration_time * std.time.ns_per_ms},
    );

    const dg = result.profile.g_grid[1] - result.profile.g_grid[0];
    var sum_under_curve: Dual.T = 0;
    for (result.profile.flux) |f| {
        sum_under_curve += f * dg;
    }
    try out.print("Sum under curve: {d:.5}\n", .{sum_under_curve});

    const filename = args.output;
    try out.print("Writing lineprofile to '{s}'\n", .{filename});
    try utils.writeFileColumns(f64, filename, &.{ result.profile.g_grid, result.profile.flux });
}
