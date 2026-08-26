const std = @import("std");
const clippy = @import("clippy");
const kerrz = @import("kerrz");

const TableArgs = @import("tables.zig").TableArgs;
const transfer_function = @import("transfer-function.zig");

const utils = @import("utils.zig");

pub const short_description = "Calculate impulse responses.";
pub const description =
    \\Calculates relativistic impulse response of an illumated accretion disc.
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
        .arg = "--dist distance",
        .argtype = f64,
        .default = "1e7",
        .help = "The observer distance in rg",
    },
    .{
        .arg = "--nthreads nthreads",
        .argtype = usize,
        .help = "The number of CPU threads to use.",
    },
    .{
        .arg = "-o/--output filename",
        .default = "impulse.dat",
        .help = "The name of the file to write the impulse response to.",
    },
    .{
        .arg = "--squash",
        .help = "Sum the impulse response in energy, such that it becomes a one-dimensional function of flux over time.",
    },
    .{
        .arg = "--zero-time time",
        .argtype = f64,
        .help = "The zero-time for the time axis of the impulse repsonse. If this is not provided, the default depends on the context: if an emissivity file is provided, the zero-time is calculated as the continuum-to-observer time for the coronal parameters in the table metadata. If no emissivity is given, then the zero time is simply the radial coordinate value.",
    },
    .{
        .arg = "--using-table path",
        .help = "Instead of caluclating the transfer function table, read and parse it from a file given by `path`. Note this has a few limitations, and will silently ignore other arguments, such as the spin or table arguments.",
    },
    .{
        .arg = "--nt n",
        .argtype = usize,
        .default = "1000",
        .help = "Number of time bins.",
    },
    .{
        .arg = "--tmin t",
        .argtype = f64,
        .default = "1.0",
        .help = "The minimum time bin.",
    },
    .{
        .arg = "--tmax t",
        .argtype = f64,
        .default = "1000.0",
        .help = "The maximum time bin.",
    },
    .{
        .arg = "--tgrid grid",
        .default = "log10",
        .help = "The grid spacing to use for the time grid. Possible options are " ++ utils.makeList(kerrz.iterators.GridSpacing),
    },
    utils.EmissivityFunction.ArgDescriptor,
} ++ utils.LineprofileIntegrationArguments ++ TableArgs);

const Dual = kerrz.DualNumber(f64, 0);
const Tool = kerrz.tools.ImpulseResponse(Dual);
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
                .r = .promote(args.dist),
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

    if (!table.table.hasField(.delta_t)) {
        return itt.throwError(
            error.MissingField,
            "Transfer function is missing `delta_t` field\n",
            .{},
        );
    }

    const zero_time = args.@"zero-time" orelse b: {
        if (emissivity.params) |params| {
            const corona = try params.toCoronalModel(allocator, table.table.tracer.metric);
            defer corona.deinit(allocator);
            const sol = try kerrz.continuum.continuumGeodesic(
                Dual,
                table.table.tracer.metric,
                table.table.tracer.x_obs,
                corona,
            );
            break :b sol.delta_t.x;
        }
        break :b 0;
    };

    try out.print("Zero time set to {d}\n", .{zero_time});

    const tool: Tool = .{
        .table = table.table,
        .lp_opts = .{
            .emissivity = emissivity.profile,
            .num_r_steps = args.nrsteps,
            .r_step_grid = try .fromString(args.rstepgrid),
            .zero_time = zero_time,
            .r_min = args.rin orelse 0,
            .r_max = args.rout,
        },
        .num_g = args.ng,
        .num_fine_g = args.ngstar,
        .num_t = args.nt,
        .t_grid = try .fromString(args.tgrid),
        .t_max = args.tmax,
        .t_min = args.tmin,
    };

    var result = try tool.run(
        allocator,
        .{ .out = out, .thread_chunk_size = 1 },
    );
    defer result.deinit(allocator);

    try out.print(
        "Integration time on single thread: {D}\n",
        .{result.integration_time * std.time.ns_per_ms},
    );

    const filename = args.output;

    if (args.squash) {
        const squashed = try result.impulse_reponse.sumColumns(allocator);
        defer allocator.free(squashed);

        // Sanity check that we summed over the correct axis:
        std.debug.assert(squashed.len == result.t_grid.len);

        // Normalise so that the repsonse sums to one:
        var total: Dual.T = 0;
        for (squashed) |s| total += s;
        if (total > 0) {
            for (squashed) |*s| s.* /= total;
        }

        try out.print("Writing squashed impulse response to '{s}'\n", .{filename});
        try utils.writeFileColumns(f64, filename, &.{ result.t_grid, squashed });
    } else {
        try out.print("Writing impulse response matrix to '{s}'\n", .{filename});
        try utils.writeFile(filename, MatrixWriter{ .matrix = result.impulse_reponse });
    }
}

const MatrixWriter = struct {
    matrix: kerrz.Matrix(f64),

    pub fn writeAll(self: MatrixWriter, out: *std.Io.Writer) !void {
        for (0..self.matrix.n_rows) |row| {
            for (0..self.matrix.n_cols) |col| {
                if (col == 0) {
                    try out.print("{d}", .{self.matrix.get(row, col)});
                } else {
                    try out.print(" {d}", .{self.matrix.get(row, col)});
                }
            }
            try out.writeAll("\n");
        }
    }
};
