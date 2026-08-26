const std = @import("std");
const test_helper = @import("test/build_helpers.zig");

fn getGitHash(allocator: std.mem.Allocator) ?[]const u8 {
    var proc = std.process.Child.init(&.{ "git", "rev-parse", "HEAD" }, allocator);
    proc.stdout_behavior = .Pipe;
    _ = proc.spawn() catch return null;
    const hash = proc.stdout.?.readToEndAlloc(allocator, 1024) catch return null;
    errdefer allocator.free(hash);
    _ = proc.wait() catch return null;
    return std.mem.trim(u8, hash, "\n \t");
}

fn getVersion(hash: []const u8) std.SemanticVersion {
    return std.SemanticVersion{
        .major = 0,
        .minor = 3,
        .patch = 4,
        .build = hash,
    };
}

fn addTracy(b: *std.Build, step: *std.Build.Step.Compile) !void {
    step.want_lto = false;
    step.addCSourceFile(
        .{
            .file = b.path("../tracy/public/TracyClient.cpp"),
            .flags = &.{
                "-DTRACY_ENABLE=1",
                "-fno-sanitize=undefined",
            },
        },
    );
    step.linkLibCpp();
}

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const tracy = b.option(
        []const u8,
        "tracy",
        "Enable Tracy integration. Supply path to Tracy source",
    );
    const llvm = b.option(
        bool,
        "llvm",
        "Compile with LLVM.",
    );

    const zad = b.dependency("zad", .{ .target = target, .optimize = optimize });
    const clippy = b.dependency("clippy", .{ .target = target, .optimize = optimize });
    const rootsolve = b.dependency("rootsolve", .{ .target = target, .optimize = optimize });
    const zfits = b.dependency("zfits", .{ .target = target, .optimize = optimize });
    const dinterp = b.dependency("dinterp", .{ .target = target, .optimize = optimize });

    const hash = getGitHash(b.allocator);
    defer if (hash) |h| b.allocator.free(h);

    const version = getVersion(hash orelse "no-git");

    var opts = b.addOptions();
    opts.addOption(f64, "test_numerical_tolerance", 1e-6);
    opts.addOption(std.SemanticVersion, "version", version);
    opts.addOption(bool, "enable_tracy", tracy != null);

    const mod = b.addModule("kerrz", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .imports = &.{
            .{ .name = "zad", .module = zad.module("zad") },
            .{ .name = "rootsolve", .module = rootsolve.module("rootsolve") },
            .{ .name = "zfits", .module = zfits.module("zfits") },
            .{ .name = "dinterp", .module = dinterp.module("dinterp") },
            .{ .name = "options", .module = opts.createModule() },
        },
    });

    const lib = b.addLibrary(.{
        .name = "kerrz",
        .linkage = .dynamic,
        .root_module = b.createModule(.{
            .root_source_file = b.path("./wrappers/interface.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "kerrz", .module = mod },
            },
            .pic = true,
        }),
    });
    lib.linkLibC();
    lib.addIncludePath(b.path("wrappers"));
    lib.installHeader(
        b.path("wrappers/kerrz.h"),
        "kerrz.h",
    );
    lib.installHeader(
        b.path("wrappers/kerrz.f90"),
        "kerrz.f90",
    );

    const lib_install_step = b.addInstallArtifact(lib, .{});
    const lib_step = b.step("lib", "Compile the shared library and wrappers.");
    lib_step.dependOn(&lib_install_step.step);

    const exe = b.addExecutable(.{
        .name = "kerrz",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "kerrz", .module = mod },
                .{ .name = "clippy", .module = clippy.module("clippy") },
                .{ .name = "options", .module = opts.createModule() },
            },
        }),
        .use_llvm = llvm,
    });

    if (tracy) |tracy_path| {
        const client_cpp = b.pathJoin(
            &[_][]const u8{ tracy_path, "public", "TracyClient.cpp" },
        );

        const tracy_c_flags: []const []const u8 = &.{
            "-DTRACY_ENABLE=1",
            "-fno-sanitize=undefined",
        };

        exe.root_module.addIncludePath(.{ .cwd_relative = tracy_path });
        exe.root_module.addCSourceFile(.{
            .file = .{ .cwd_relative = client_cpp },
            .flags = tracy_c_flags,
        });
        exe.root_module.link_libc = true;
        exe.root_module.link_libcpp = true;
    }

    b.installArtifact(exe);

    const run_step = b.step("run", "Run the app");

    const run_cmd = b.addRunArtifact(exe);
    run_step.dependOn(&run_cmd.step);

    run_cmd.step.dependOn(b.getInstallStep());

    if (b.args) |args| {
        run_cmd.addArgs(args);
    }

    const mod_tests = b.addTest(.{
        .root_module = mod,
    });

    const run_mod_tests = b.addRunArtifact(mod_tests);

    // Exe tests:
    const exe_tests = b.addTest(.{
        .root_module = exe.root_module,
    });

    const run_exe_tests = b.addRunArtifact(exe_tests);

    const test_step = b.step("test", "Run tests");
    test_step.dependOn(&run_mod_tests.step);
    test_step.dependOn(&run_exe_tests.step);

    // Regression testing:
    const regression_tests = b.addExecutable(.{
        .name = "kerrz-regressions",
        .root_module = b.createModule(.{
            .root_source_file = b.path("test/regression.zig"),
            .optimize = .ReleaseSafe,
            .target = target,
            .imports = &.{
                .{ .name = "kerrz", .module = mod },
            },
        }),
    });

    const run_regression_tests = b.addRunArtifact(regression_tests);

    const test_regression_step = b.step("test-regression", "Run regression tests");
    test_regression_step.dependOn(&run_regression_tests.step);

    // Documentation:
    const docs_step = b.step("docs", "Generate HTML documentation");
    const install_std_docs = b.addInstallDirectory(.{
        .source_dir = mod_tests.getEmittedDocs(),
        .install_dir = .prefix,
        .install_subdir = "docs",
    });

    docs_step.dependOn(&install_std_docs.step);

    // Wrapper tests
    test_helper.addWrapperTests(
        .{
            .b = b,
            .kerrz_lib = lib,
            .optimize = optimize,
            .target = target,
        },
        "lib-test-build",
        "lib-test",
    );
}
