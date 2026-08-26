const std = @import("std");
const clippy = @import("clippy");
const kerrz = @import("kerrz");

const transfer_function = @import("./transfer-function.zig");
const CTFT = kerrz.TransferFunctionContext(Dual);

const utils = @import("utils.zig");

pub const short_description = "Calculate non-trivial tables.";
pub const description =
    \\Calculate non-trivial tables, such as the Cunningham's transfer function
    \\tables. These are then saved to file for use in other places.
    \\
    \\If only photon trajectories and results are needed, many of the commands in
    \\kerrz have various `--output` options for exporting trivial tables. See for
    \\example `--help emissivity`.
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
        .arg = "-o/--output path",
        .help = "The filepath to write the table to.",
        .default = "ctftable.fits",
    },
    .{
        .arg = "--grid path",
        .help = "A path to a grid of spins and observer inclinations (in degrees) for which to calculate the transfer functions for. The grid should contains two rows (spin and inclination), with values seperated by commas. If this argument is passed, `--spin` and `--incl` are ignored.",
    },
} ++ TableArgs);

pub const TableArgs = [_]clippy.ArgumentDescriptor{
    .{
        .arg = "--nradii n",
        .argtype = usize,
        .default = "100",
        .help = "The number of radii to calculate transfer functions for.",
    },
    .{
        .arg = "--rin r_in",
        .argtype = f64,
        .help = "The inner radius. Defaults to the event horizon.",
    },
    .{
        .arg = "--rout r_out",
        .argtype = f64,
        .default = "1000",
        .help = "The outer radius.",
    },
    .{
        .arg = "--rgrid grid",
        .default = "log_sigmoid",
        .help = "The grid spacing to use for the radial coordinates. Possible options are " ++ utils.makeList(kerrz.iterators.GridSpacing),
    },
    .{
        .arg = "--fields names",
        .default = "f,delta_t,local_theta",
        .help = "A comma-seperated list of fields that should be serialised in the fits table. Possible values are " ++ utils.makeList(kerrz.transfer_tables.SelectedField),
    },
} ++ transfer_function.CTFArguments;

const Dual = kerrz.DualNumber(f64, 0);
const Tool = kerrz.tools.TransferFunctionTable(Dual);

pub fn run(out: *std.Io.Writer, allocator: std.mem.Allocator, itt: *clippy.ArgumentIterator) !void {
    const args = try Args.initParseAll(itt, .{});
    const num_threads = utils.getNumThreads(args.nthreads);

    const selected_fields = try utils.parseCommaSeperatedEnum(
        allocator,
        kerrz.transfer_tables.SelectedField,
        args.fields,
    );
    defer allocator.free(selected_fields);

    const fits = try kerrz.FitsFile.init(allocator);
    defer fits.deinit();

    var threads = try kerrz.ThreadMap.init(
        allocator,
        .{ .num_threads = num_threads },
    );
    defer threads.deinit();

    var arena = std.heap.ArenaAllocator.init(allocator);
    defer arena.deinit();
    const name_allocator = arena.allocator();

    if (args.grid) |grid_path| {
        const grid = try TableGrid.parse(allocator, grid_path);
        defer grid.deinit(allocator);

        try out.print(
            "Read {d} x {d} grid from '{s}'\n",
            .{ grid.spins.len, grid.inclinations.len, grid_path },
        );

        // Write the HDUs containing the parameters.
        const a_hdu = try fits.newHdu("spins", .binary_table);
        try a_hdu.data.binary_table.appendColumn(.{
            .label = "a",
            .comment = "Spin of the black hole",
            .units = "J/M",
            .units_comment = "Dimensionless standard units",
        });
        for (grid.spins) |spin| {
            const row = try a_hdu.data.binary_table.addRow();
            row.cols[0].one.float_32 = @floatCast(spin);
        }
        const mu = try fits.newHdu("incl", .binary_table);
        try mu.data.binary_table.appendColumn(.{
            .label = "incl",
            .comment = "The observer inclination",
            .units = "degree",
        });
        for (grid.inclinations) |incl| {
            const row = try mu.data.binary_table.addRow();
            row.cols[0].one.float_32 = @floatCast(incl);
        }

        var key_buffer: [16]u8 = undefined;

        const total_number = grid.spins.len * grid.inclinations.len;

        const time_now = std.time.milliTimestamp();

        for (grid.spins, 0..) |spin, j| {
            for (grid.inclinations, 0..) |incl, i| {
                const incl_deg = incl;
                const index = i + j * grid.inclinations.len + 1;
                if (index != 1) {
                    try out.writeAll("\r\x1b[K");
                }
                try out.print(
                    "{d: >4} of {d}: a = {d:.5}, incl = {d:.5}",
                    .{ index, total_number, spin, incl_deg },
                );
                try evaluateTableHdu(
                    null,
                    allocator,
                    fits,
                    threads,
                    spin,
                    incl_deg,
                    args,
                    selected_fields,
                );
                const name = try std.fmt.bufPrint(
                    &key_buffer,
                    "CTF_{d}_{d}",
                    .{ i, j },
                );
                try fits.hdus[fits.hdus.len - 1].setName(
                    try name_allocator.dupe(u8, name),
                    null,
                );
            }
        }

        const duration = std.time.milliTimestamp() - time_now;

        // Clear the indicator line
        _ = try out.write("\r\x1b[K");
        try out.print(
            "Calculated {d} tables in {D} ({D} per table)\n",
            .{
                total_number,
                duration * std.time.ns_per_ms,
                @divFloor(duration * std.time.ns_per_ms, @as(i64, @intCast(total_number))),
            },
        );

        // Set the table type in the primary HDU:
        try fits.hdus[0].appendHeaderRecord(
            "KZTYPE",
            .{
                .value = .{ .string = "ctfgrid" },
                .comment = "This file contains a grid of CTF tables",
            },
        );
    } else {
        try evaluateTableHdu(
            out,
            allocator,
            fits,
            threads,
            args.spin,
            args.incl,
            args,
            selected_fields,
        );

        // Set the table type in the primary HDU:
        try fits.hdus[0].appendHeaderRecord(
            "KZTYPE",
            .{
                .value = .{ .string = "ctf" },
                .comment = "This file contains a single CTF table",
            },
        );
    }

    try out.print(
        "Writing transfer function table to '{s}'\n",
        .{args.output},
    );
    try fits.save(.{ .path = args.output });
}

const TableGrid = struct {
    spins: []const f64,
    inclinations: []const f64,

    pub fn parse(allocator: std.mem.Allocator, path: []const u8) !TableGrid {
        const dir = std.fs.cwd();
        const stat = try dir.statFile(path);

        const contents = try dir.readFileAlloc(allocator, path, stat.size);
        defer allocator.free(contents);

        var line_itt = std.mem.tokenizeScalar(u8, contents, '\n');

        const spin_line = line_itt.next() orelse {
            return error.MissingSpinRow;
        };
        const incl_line = line_itt.next() orelse {
            return error.MissingInclinationRow;
        };

        const spins = try parseLine(allocator, spin_line);
        errdefer allocator.free(spins);

        const incl = try parseLine(allocator, incl_line);
        errdefer allocator.free(incl);

        return .{
            .spins = spins,
            .inclinations = incl,
        };
    }

    fn parseLine(allocator: std.mem.Allocator, line: []const u8) ![]const f64 {
        var itt = std.mem.tokenizeAny(u8, line, ", ");

        var list = std.ArrayList(f64).empty;
        defer list.deinit(allocator);

        while (itt.next()) |token| {
            const value = try std.fmt.parseFloat(f64, token);
            try list.append(allocator, value);
        }

        return try list.toOwnedSlice(allocator);
    }

    pub fn deinit(self: *const TableGrid, allocator: std.mem.Allocator) void {
        allocator.free(self.spins);
        allocator.free(self.inclinations);
    }
};

fn evaluateTableHdu(
    out: ?*std.Io.Writer,
    allocator: std.mem.Allocator,
    fits: *kerrz.FitsFile,
    threads: *kerrz.ThreadMap,
    spin: f64,
    incl: f64,
    args: Args.Parsed,
    selected_fields: []const kerrz.transfer_tables.SelectedField,
) !void {
    const metric: kerrz.KerrMetric(Dual) = .init(.one, .promote(spin));
    const x_obs: kerrz.FourVector(Dual) = .{
        .t = .zero,
        .r = .promote(args.dist),
        .th = .promote(std.math.degreesToRadians(incl)),
        .ph = .zero,
    };

    const tool: Tool = .{
        .metric = metric,
        .x_obs = x_obs,
        .tf_opts = .{
            .minimum_guess = args.@"min-guess",
            .initial_guess = args.@"initial-guess",
            .heuristic = try .fromOptionalString(args.heuristic),
            .max_points = args.nangles,
            .refine = try .parse(args.refine),
            .optimise = args.optimise,
        },
        .table_opts = .{
            .num_radii = args.nradii,
            .r_in = args.rin,
            .r_out = args.rout,
            .r_grid = try .fromString(args.rgrid),
        },
    };

    var result = try tool.run(
        allocator,
        threads,
        .{
            .out = out,
            .thread_chunk_size = 1,
            .show_progress = out != null,
        },
    );
    defer result.deinit(allocator);

    const hdu_ptr = try fits.addHdu();
    // Serialise to a FITS file.
    hdu_ptr.* = try result.table.toFITSSelected(allocator, selected_fields);
    // Add some information about how the radial grid was computed
    try hdu_ptr.appendHeaderRecord("RGRID", .{
        .value = .{ .string = args.rgrid },
        .comment = "The kerrz grid used to compute the radial axis",
    });
}
