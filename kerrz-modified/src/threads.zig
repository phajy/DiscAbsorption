const std = @import("std");
const tracy = @import("tracy.zig");

/// Used to parallelize work over threads using a batching algorithm.
pub const ThreadMap = struct {
    pub const Options = struct {
        num_threads: usize = 1,
    };

    pool: std.Thread.Pool,
    allocator: std.mem.Allocator,
    shared_index: usize = 0,
    data_len: usize = 0,
    opts: Options,
    wg: std.Thread.WaitGroup = .{},

    /// Join all threads, cleanup allocated resources, and destroy self.
    pub fn deinit(self: *ThreadMap) void {
        self.pool.deinit();
        self.allocator.destroy(self);
    }

    /// Initialize a new thread map
    pub fn init(
        allocator: std.mem.Allocator,
        opts: Options,
    ) !*ThreadMap {
        const ptr = try allocator.create(ThreadMap);
        errdefer allocator.destroy(ptr);

        ptr.allocator = allocator;
        ptr.opts = opts;
        try ptr.pool.init(
            .{
                .allocator = allocator,
                .n_jobs = opts.num_threads,
            },
        );

        return ptr;
    }

    pub const MapOptions = struct {
        chunk_size: usize = 128,
    };

    pub const ThreadId = struct { usize };

    /// In-place map: maps the function `f` onto each element of `slice`
    pub fn map(
        self: *ThreadMap,
        comptime T: type,
        slice: []T,
        ctx: anytype,
        comptime f: fn (@TypeOf(ctx), *T, usize, ThreadId) void,
        opts: MapOptions,
    ) !void {
        const Context = @TypeOf(ctx);

        const Wrapper = struct {
            parent: *ThreadMap,
            user_ctx: Context,
            data: []T,
            index: usize = 0,
            chunk_size: usize,
            id: usize,

            fn getNextIndex(w: @This()) usize {
                const tracy_ctx = tracy.trace(@src());
                defer tracy_ctx.end();
                const mut = &w.parent.pool.mutex;

                mut.lock();
                defer mut.unlock();

                const i = w.parent.shared_index;
                w.parent.shared_index += w.chunk_size;
                return i;
            }

            fn doWork(w: @This()) void {
                const tracy_ctx = tracy.trace(@src());
                defer tracy_ctx.end();
                var ind = w.index;
                while (ind <= w.data.len) {
                    const s = w.data[ind..@min(w.data.len, ind + w.chunk_size)];
                    for (s, 0..) |*v, offset| f(w.user_ctx, v, ind + offset, ThreadId{w.id});
                    ind = w.getNextIndex();
                }

                w.parent.wg.finish();
            }
        };

        self.shared_index = 0;
        self.data_len = slice.len;
        self.wg = .{};
        for (0..self.opts.num_threads) |i| {
            self.wg.start();

            const w: Wrapper = .{
                .parent = self,
                .user_ctx = ctx,
                .data = slice,
                .index = self.shared_index,
                .chunk_size = opts.chunk_size,
                .id = i,
            };

            {
                // this could be a race condition so need to lock to increment
                self.pool.mutex.lock();
                defer self.pool.mutex.unlock();
                self.shared_index += opts.chunk_size;
            }

            try self.pool.spawn(Wrapper.doWork, .{w});
        }
    }

    /// Blocks the calling thread until the work group has finished
    pub fn blockUntilDone(self: *ThreadMap) void {
        self.pool.waitAndWork(&self.wg);
    }

    /// Display a progress indicator.
    pub fn blockWithProgress(self: *ThreadMap, options: ProgressOptions) !void {
        // TODO: in 0.16 replace this with IO concurrency
        var write_buffer: [128]u8 = undefined;
        var stdout = std.fs.File.stdout().writer(&write_buffer);
        const writer = &stdout.interface;

        var previous_time = std.time.milliTimestamp();
        while (!self.wg.isDone()) {
            const percent = @as(f64, @floatFromInt(self.shared_index)) /
                @as(f64, @floatFromInt(@max(1, self.data_len)));

            if (percent > 1) continue;

            const now = std.time.milliTimestamp();

            // update every 5th of a second
            if (now - previous_time > options.refresh) {
                previous_time = now;
                try writeProgress(writer, percent, options);
                try writer.flush();
            }
        }

        try writeProgress(writer, 1.0, options);
        _ = try writer.write("\r\x1b[K");
        try writer.flush();
    }
};

pub const ProgressOptions = struct {
    refresh: i64 = std.time.ms_per_s / 5,
    length: usize = 30,
};

fn writeProgress(out: *std.Io.Writer, fraction: f64, options: ProgressOptions) !void {
    _ = try out.write("\rProgress: [");
    const len: usize = @intFromFloat(fraction * @as(f64, @floatFromInt(options.length)));
    _ = try out.splatBytes("=", len);
    _ = try out.splatBytes(" ", options.length - 1 -| len);
    try out.print("] ({d}%)", .{@as(usize, @intFromFloat(fraction * 100))});
}
