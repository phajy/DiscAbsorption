/// These are regression tests that check that `kerrz` as the Zig library is
/// functioning correctly on things that are not strictly unit tests.
///
/// These are run in `ReleaseSafe` mode and can require slightly heavier
/// computations.
const kerrz = @import("kerrz");
const std = @import("std");

const logger = std.log.scoped(.regression);

const Dual0 = kerrz.DualNumber(f64, 0);

const kerr_metric: kerrz.KerrMetric(Dual0) = .init(.one, .promote(0.998));
const schwz_metric: kerrz.KerrMetric(Dual0) = .init(.one, .promote(0.01));

var allocator: std.mem.Allocator = undefined;
var threads: *kerrz.ThreadMap = undefined;

pub fn main() !void {
    var gpa = std.heap.DebugAllocator(.{}){};
    defer _ = gpa.deinit();
    allocator = gpa.allocator();

    threads = try kerrz.ThreadMap.init(
        allocator,
        .{ .num_threads = 4 },
    );
    defer threads.deinit();

    // Run the tests scripts.
    try runTransferFunctionTests();
    try runObserverImages();
}

fn distantObserver(comptime T: type, r: T.T, th: T.T) kerrz.FourVector(T) {
    return .{
        .t = .zero,
        .r = .promote(r),
        .th = .promote(std.math.degreesToRadians(th)),
        .ph = .zero,
    };
}

// Transfer functions ------------------------------------------------------- //

fn testTransferFunctions(tool: kerrz.tools.TransferFunction(Dual0)) !void {
    var result = try tool.run(allocator);
    defer result.deinit(allocator);
}

fn runTransferFunctionTests() !void {
    inline for (.{ kerr_metric, schwz_metric }) |metric| {
        // Testing the limits of the angles.
        const angles = [_]f64{ 2, 5, 45, 88 };
        const radii = [_]f64{ metric.isco.x, 10, 500.0, 1000.0 };
        for (angles) |angle| {
            for (radii) |r| {
                testTransferFunctions(
                    .{
                        .metric = metric,
                        .r_target = r,
                        .x_obs = distantObserver(Dual0, 1e6, angle),
                        .opts = .{},
                    },
                ) catch |err| {
                    logger.err(
                        "Failed for a={d}, r={d}, th={d}",
                        .{ metric.a.x, r, angle },
                    );
                    return err;
                };
            }
        }
    }
}

// Images ------------------------------------------------------------------- //

fn testImage(tool: kerrz.tools.Image(Dual0)) !void {
    const alpha = try kerrz.iterators.Grid(Dual0.T).linear.fill(
        allocator,
        -10,
        10,
        256,
    );
    defer allocator.free(alpha);
    const beta = try kerrz.iterators.Grid(Dual0.T).linear.fill(
        allocator,
        -10,
        10,
        256,
    );
    defer allocator.free(beta);

    var result = try tool.runImpactParameters(
        allocator,
        threads,
        alpha,
        beta,
        .{
            .thread_chunk_size = 512,
            .show_progress = false,
        },
    );
    defer result.deinit(allocator);
}

fn runObserverImages() !void {
    inline for (.{ kerr_metric, schwz_metric }) |metric| {
        // Testing the limits of the angles.
        const angles = [_]f64{ 2, 5, 45, 88 };
        for (angles) |angle| {
            testImage(.{
                .metric = metric,
                .x_obs = distantObserver(Dual0, 1e5, angle),
                .image_opts = .{},
            }) catch |err| {
                logger.err(
                    "Failed for a={d}, th={d}",
                    .{ metric.a.x, angle },
                );
                return err;
            };
        }
    }
}
