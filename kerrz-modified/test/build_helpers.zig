const std = @import("std");

const ALL_WRAPPER_TESTS = [_]struct {
    subpath: []const u8,
    name: []const u8,
}{
    .{ .subpath = "test/test_shadow.c", .name = "test_shadow" },
    .{ .subpath = "test/test_emissivity.c", .name = "test_emissivity" },
    .{ .subpath = "test/test_transfer_function.c", .name = "test_transfer_function" },
    .{ .subpath = "test/test_lineprofile.c", .name = "test_lineprofile" },
    .{ .subpath = "test/test_continuum.c", .name = "test_continuum" },
};

pub const Helper = struct {
    target: std.Build.ResolvedTarget,
    optimize: std.builtin.OptimizeMode,
    kerrz_lib: *std.Build.Step.Compile,
    b: *std.Build,
};

fn buildC(h: Helper, subpath: []const u8, name: []const u8) *std.Build.Step.Compile {
    const exe = h.b.addExecutable(.{
        .name = name,
        .root_module = h.b.createModule(.{
            .target = h.target,
            .optimize = h.optimize,
            .sanitize_c = .full,
        }),
    });
    exe.linkLibC();
    exe.addCSourceFile(.{ .file = h.b.path(subpath) });
    exe.linkLibrary(h.kerrz_lib);
    return exe;
}

pub fn addWrapperTests(h: Helper, add_artifact_name: []const u8, run_name: []const u8) void {
    const build_all_wrapper_tests = h.b.step(add_artifact_name, "Build and install the C library wrapper tests.");
    const run_all_wrapper_tests = h.b.step(run_name, "Run the C library wrapper tests.");
    for (&ALL_WRAPPER_TESTS) |wt| {
        const compile = buildC(h, wt.subpath, wt.name);

        const install_step = h.b.addInstallArtifact(
            compile,
            .{ .dest_dir = .{ .override = .{ .custom = "bin/wrappers_tests" } } },
        );
        build_all_wrapper_tests.dependOn(&install_step.step);

        const run_step = h.b.addRunArtifact(compile);
        run_all_wrapper_tests.dependOn(&run_step.step);
    }
}
