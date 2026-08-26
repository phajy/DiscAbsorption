const std = @import("std");
const kerrz = @import("kerrz");
const clippy = @import("clippy");
const utils = @import("utils.zig");

pub const short_description = "Print information about the mapper.";
pub const description =
    \\Below are listed all of the values that may be given to a `--map`
    \\argument. Various operations may be applied to the values, which are
    \\specified via a `name.op`, for example, `time.log` will give the
    \\logarithm (base 10) of the time coordinate.
    \\
    \\Derivatives are also supported with respect to the coordinate axes, i.e.
    \\the impact parameters of `image` and the sky angles for `sky`. These are
    \\specified before the operator, such as `redshift.dx.log` will give the
    \\logarithm of the derivative of the redshift with respect to the `x`
    \\coordinate axes.
    \\
    \\The derivative-based values will calculate the derivative
    \\components of each geodesic, which can almost double the runtime under
    \\certain circumstances.
;

pub const Args = clippy.Arguments(&.{});

pub fn printHelp(writer: *std.Io.Writer) !void {
    try writer.writeAll("\nAvailable operations:\n");
    inline for (@typeInfo(kerrz.Mapper.Operator).@"enum".fields) |field| {
        try writer.print("- {s}\n", .{field.name});
    }
    try writer.writeAll("\nMapper values:\n");
    for (kerrz.Mapper.Fields) |item| {
        try writer.print("- {s: <15} {s}\n", .{ item.name, item.description });
    }
}

pub fn run(out: *std.Io.Writer, allocator: std.mem.Allocator, itt: *clippy.ArgumentIterator) !void {
    _ = allocator;
    const args = try Args.initParseAll(itt, .{});
    _ = args;
    try out.print("{s}\n", .{description});
    try printHelp(out);
}
