const std = @import("std");
const clippy = @import("clippy");
const kerrz = @import("kerrz");

const utils = @import("./utils.zig");

pub const short_description = "Generate grids.";
pub const description =
    \\Generate various different types of grids.
    \\
    \\This tool is designed for generating the parameter grids used as the argument
    \\to e.g. `--grid` in the `table` command. It does not perform any complex
    \\calculations and is primarily a utility command.
    \\
;

pub const Args = clippy.Arguments(&[_]clippy.ArgumentDescriptor{
    .{
        .arg = "--spacing spacing",
        .default = "linear",
        .help = "The grid spacing to use. Possible options are " ++ utils.makeList(kerrz.iterators.GridSpacing),
    },
    .{
        .arg = "-n/--len len",
        .argtype = usize,
        .default = "10",
        .help = "The length of the grid to generate.",
    },
    .{
        .arg = "--min min",
        .argtype = f64,
        .default = "0",
        .help = "The minimum value of the grid.",
    },
    .{
        .arg = "--max max",
        .argtype = f64,
        .default = "1",
        .help = "The maximum value of the grid.",
    },
    .{
        .arg = "--invert",
        .help = "Invert or reverse the grid spacing. Does not change the order of the elements. To change the order, switch `min` and `max`.",
    },
});

pub fn run(out: *std.Io.Writer, allocator: std.mem.Allocator, itt: *clippy.ArgumentIterator) !void {
    _ = allocator;
    const args = try Args.initParseAll(itt, .{});

    const grid = try kerrz.iterators.Grid(f64).fromString(args.spacing);

    var grid_itt: kerrz.iterators.GridIterator(f64) = undefined;
    if (args.invert) {
        grid_itt = grid.iterator(args.max, args.min, args.len);
        grid_itt.invert = args.invert;
    } else {
        grid_itt = grid.iterator(args.min, args.max, args.len);
    }

    while (grid_itt.next()) |item| {
        try out.print("{d} ", .{item});
    }
    try out.writeAll("\n");
}
