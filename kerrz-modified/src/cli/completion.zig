const std = @import("std");
const clippy = @import("clippy");

pub const short_description = "Generate shell completion helpers.";
pub const description = short_description ++ "\n";
pub const Args = clippy.Arguments(&.{});

const Commands = @import("../main.zig").Command;

fn wrapCommands() type {
    const Field = std.builtin.Type.UnionField;
    comptime var fields: []const Field = &.{};
    inline for (@typeInfo(Commands).@"union".fields) |f| {
        const field: Field = .{
            .alignment = f.alignment,
            .name = f.name,
            .type = f.type.Args,
        };
        fields = fields ++ .{field};
    }
    return clippy.Commands(@Type(.{
        .@"union" = .{
            .decls = &.{},
            .fields = fields,
            .layout = .auto,
            .tag_type = @typeInfo(Commands).@"union".tag_type,
        },
    }));
}

const WrappedCommands = wrapCommands();

/// Execute the completion helper
pub fn run(out: *std.Io.Writer, allocator: std.mem.Allocator, _: *clippy.ArgumentIterator) !void {
    const completion = try WrappedCommands.generateCompletion(
        allocator,
        .{ .function_name = "kerrz" },
    );
    defer allocator.free(completion);

    try out.writeAll("#compdef _kerrz kerrz\n");
    try out.writeAll(completion);
    try out.writeAll("\n_arguments_kerrz\n");
}
