const std = @import("std");
const clippy = @import("clippy");
const kerrz = @import("kerrz");
const el = kerrz.elliptic_integrals;

const utils = @import("./cli/utils.zig");

test "all" {
    _ = Command;
    _ = utils;
}

const HelpArguments = clippy.Arguments(
    &.{
        .{
            .arg = "--help",
            .help = "Print this help message or help for a specific command.",
        },
        .{
            .arg = "--version",
            .help = "Print version information.",
        },
    },
);

pub const Command = union(enum) {
    calc: @import("cli/calculator.zig"),
    continuum: @import("cli/continuum-spectrum.zig"),
    emissivity: @import("cli/emissivity-curves.zig"),
    grid: @import("cli/grid.zig"),
    image: @import("cli/observer-image.zig"),
    impulse: @import("cli/impulse-reponse.zig"),
    lineprof: @import("cli/lineprofiles.zig"),
    sky: @import("cli/sky-image.zig"),
    table: @import("cli/tables.zig"),
    tf: @import("cli/transfer-function.zig"),
    trace: @import("cli/trace.zig"),
    mapper: @import("cli/mapper.zig"),
    completion: @import("cli/completion.zig"),
};

const CommandName = std.meta.FieldEnum(Command);

fn printVersion(writer: *std.Io.Writer) !void {
    try writer.writeAll("kerrz GPL 3.0 version: ");
    try kerrz.version.format(writer);
    try writer.writeAll("\nhttps://git.sr.ht/~fjebaker/kerrz\n");
}

fn printHelp(writer: *std.Io.Writer, command: ?CommandName) !void {
    try printVersion(writer);
    try writer.print(
        "\nUsage: kerrz {{{s}}} [--flags and positional arguments]",
        .{if (command) |cmd| @tagName(cmd) else "command"},
    );
    const command_fields = @typeInfo(Command).@"union".fields;
    if (command) |cmd| {
        const cmd_string = @tagName(cmd);
        inline for (command_fields) |field| {
            if (std.mem.eql(u8, field.name, cmd_string)) {
                try writer.print("\n\n{s}\n", .{field.type.description});
            }
        }

        if (cmd == .mapper) {
            const Mapper = command_fields[std.meta.fieldIndex(Command, "mapper").?];
            try Mapper.type.printHelp(writer);
        }

        inline for (command_fields) |field| {
            if (std.mem.eql(u8, field.name, cmd_string)) {
                if (@typeInfo(field.type.Args.Parsed).@"struct".fields.len > 0) {
                    try writer.writeAll("Arguments:\n");
                    try field.type.Args.writeHelp(writer, .{});
                }
            }
        }
    } else {
        try writer.writeAll("\n\nGeneral flags:\n");
        try HelpArguments.writeHelp(writer, .{});
        try writer.writeAll("\nCommands:\n");
        inline for (command_fields) |field| {
            try writer.print(
                "- {s: <11} {s}\n",
                .{ field.name, field.type.short_description },
            );
        }
    }
    try writer.flush();
}

const Dual = kerrz.DualNumber(f64, 0);

pub fn main() !void {
    var gpa: std.heap.DebugAllocator(.{}) = .init;
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const raw_args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, raw_args);

    var itt = clippy.ArgumentIterator.init(raw_args[1..]);
    var help_itt = itt;

    // first we parse to see if help was given
    const help_parsed = try HelpArguments.initParseAll(&help_itt, .{ .forgiving = true });

    var stdout = std.fs.File.stdout();
    var stdout_writer = stdout.writer(&.{});
    const writer = &stdout_writer.interface;

    if (help_parsed.help) {
        // parse all the remaining arguments to see if there is a positional
        // argument
        var cmd: ?CommandName = null;
        while (itt.next() catch null) |arg| {
            if (!arg.flag) {
                cmd = std.meta.stringToEnum(CommandName, arg.string);
                if (cmd != null) break;
            }
        }

        try printHelp(writer, cmd);
        return std.process.cleanExit();
    }

    // Only print version information if it's the only argument.
    if (itt.argCount() == 1 and help_parsed.version) {
        try printVersion(writer);
        return std.process.cleanExit();
    }

    const first_arg = (try itt.next()) orelse {
        utils.writeError(
            error.MissingCommand,
            "No command specified.\nUse `--help` to display all commands.",
            .{},
        ) catch {};
        return std.process.cleanExit();
    };

    if (first_arg.flag) {
        utils.writeError(
            error.UnknownCommand,
            "First argument must be a command, not a flag.\nUse `--help` to display all commands.",
            .{},
        ) catch {};
        return std.process.cleanExit();
    }

    // This acts as a cheap and chearful way to see if the argument is valid or
    // not.
    _ = std.meta.stringToEnum(
        CommandName,
        first_arg.string,
    ) orelse {
        utils.writeError(
            error.UnknownCommand,
            "Unknown command '{s}'.\nUse `--help` to display all commands.",
            .{first_arg.string},
        ) catch {};
        return std.process.cleanExit();
    };

    // Run the command.
    inline for (@typeInfo(Command).@"union".fields) |field| {
        if (std.mem.eql(u8, first_arg.string, field.name)) {
            try field.type.run(writer, allocator, &itt);
        }
    }

    // Don't forget to flush!
    try writer.flush();
    return std.process.cleanExit();
}
