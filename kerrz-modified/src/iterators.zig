const std = @import("std");

/// Linear iterator over a range of values.
pub fn RangeIterator(comptime T: type) type {
    return struct {
        const Self = @This();
        delta: T,
        current: T,
        remaining: usize,

        pub fn next(self: *Self) ?T {
            if (self.remaining > 0) {
                const v = self.current;
                self.current += self.delta;
                self.remaining -= 1;
                return v;
            } else return null;
        }

        pub fn init(min: T, max: T, N: usize) Self {
            const delta = (max - min) / @as(T, @floatFromInt(N - 1));
            return .{
                .delta = delta,
                .remaining = N,
                .current = min,
            };
        }

        pub fn drain(self: *Self, allocator: std.mem.Allocator) ![]T {
            var out = try allocator.alloc(T, self.remaining);
            var i: usize = 0;
            while (self.next()) |v| {
                out[i] = v;
                i += 1;
            }
            return out;
        }
    };
}

/// The different grid spacings that a `Grid` can be.
pub const GridSpacing = enum {
    log10,
    inverse,
    linear,
    geometric,
    sigmoid,
    log_sigmoid,
};

/// Grid types with their state.
pub fn Grid(comptime T: type) type {
    return struct {
        const Self = @This();

        pub const Arguments = union(enum) {
            no_arguments: void,
            sigmoid: struct {
                a: T = 3.0,
                b: T = 3.0,
                c: T = 1.0,
            },

            pub const none: Arguments = .{ .no_arguments = {} };
        };

        spacing: GridSpacing,
        args: Arguments,

        fn init(spacing: GridSpacing, args: Arguments) Self {
            return .{ .spacing = spacing, .args = args };
        }

        pub const log10: Self = .init(.log10, .none);
        pub const inverse: Self = .init(.inverse, .none);
        pub const linear: Self = .init(.linear, .none);
        pub const geometric: Self = .init(.geometric, .none);
        pub fn sigmoid(args: std.meta.fieldInfo(Arguments, .sigmoid).type) Self {
            return .init(.sigmoid, .{ .sigmoid = args });
        }
        pub fn log_sigmoid(args: std.meta.fieldInfo(Arguments, .sigmoid).type) Self {
            return .init(.log_sigmoid, .{ .sigmoid = args });
        }

        pub fn fromString(string: []const u8) !Self {
            var spacing_string = string;
            var arg_buffer: [3]?T = .{null} ** 3;
            if (std.mem.indexOfScalar(u8, string, '(')) |split| {
                if (string[string.len - 1] != ')') return error.InvalidString;
                spacing_string = string[0..split];
                var itt = std.mem.tokenizeAny(u8, string[split..], "(,)");

                var i_arg: usize = 0;
                while (itt.next()) |token| : (i_arg += 1) {
                    arg_buffer[i_arg] = try std.fmt.parseFloat(T, token);
                }
            }

            const spacing = std.meta.stringToEnum(GridSpacing, spacing_string) orelse
                return error.NoSuchGridSpacing;

            switch (spacing) {
                .log_sigmoid, .sigmoid => {
                    var self = Self{ .spacing = spacing, .args = .{ .sigmoid = .{} } };
                    if (arg_buffer[0]) |arg| self.args.sigmoid.a = arg;
                    if (arg_buffer[1]) |arg| self.args.sigmoid.b = arg;
                    if (arg_buffer[2]) |arg| self.args.sigmoid.c = arg;
                    return self;
                },
                else => return .init(spacing, .none),
            }
        }

        const State = union {
            default: Arguments,
            geometric: struct {
                base: T,
            },
        };

        fn initState(self: Self, low: T, high: T, len: usize) State {
            switch (self.spacing) {
                .geometric => {
                    // The constant multiplicative factor
                    const k = std.math.pow(
                        T,
                        (high / low),
                        1 / (@as(T, @floatFromInt(len))),
                    );
                    return .{ .geometric = .{ .base = k } };
                },
                else => return .{ .default = self.args },
            }
        }

        /// Get an iterator for this grid.
        pub fn iterator(
            self: Self,
            low: T,
            high: T,
            len: usize,
        ) GridIterator(T) {
            const state = self.initState(low, high, len);

            const range_iterator: RangeIterator(T) =
                switch (self.spacing) {
                    .linear => .init(
                        low,
                        high,
                        len,
                    ),
                    .geometric => .init(
                        0.0,
                        @floatFromInt(len),
                        len,
                    ),
                    .log_sigmoid, .sigmoid => .init(
                        0.0,
                        1.0,
                        len,
                    ),
                    .log10 => .init(
                        std.math.log10(low),
                        std.math.log10(high),
                        len,
                    ),
                    .inverse => .init(
                        1.0 / (low),
                        1.0 / (high),
                        len,
                    ),
                };

            return .{
                .spacing = self.spacing,
                .state = state,
                .itt = range_iterator,
                .min = low,
                .max = high,
            };
        }

        /// Allocate and fill the grid with values between `low` and `high`, with
        /// `size` elements.
        ///
        /// Caller owns memory.
        pub fn fill(
            self: Self,
            allocator: std.mem.Allocator,
            low: T,
            high: T,
            len: usize,
        ) ![]T {
            var range = self.iterator(low, high, len);
            return try range.drain(allocator);
        }
    };
}

pub fn GridIterator(comptime T: type) type {
    return struct {
        const Self = @This();

        itt: RangeIterator(T),
        state: Grid(T).State,
        spacing: GridSpacing,
        min: T,
        max: T,
        invert: bool = false,

        /// Get the next element of the grid.
        pub fn next(self: *Self) ?T {
            const val = self.itt.next() orelse
                return null;

            const grid_val = self.applyGridSpecifics(val);

            if (self.invert) {
                return self.max - grid_val + self.min;
            }

            return grid_val;
        }

        fn applyGridSpecifics(self: *Self, val: T) T {
            switch (self.spacing) {
                .linear => return val,
                .log10 => return std.math.pow(T, 10, val),
                .inverse => return 1 / val,
                .geometric => {
                    const state = self.state.geometric;
                    return self.min * std.math.pow(T, state.base, val);
                },
                .log_sigmoid, .sigmoid => {
                    const state = self.state.default.sigmoid;
                    const y = std.math.pow(
                        T,
                        1 - std.math.pow(
                            T,
                            1 - std.math.pow(T, std.math.clamp(val, 0, 1), state.a),
                            state.b,
                        ),
                        state.c,
                    );

                    if (self.spacing == .log_sigmoid) {
                        const grid_val = std.math.log10(
                            self.max / self.min,
                        ) * y + std.math.log10(self.min);
                        return std.math.pow(T, 10, grid_val);
                    }

                    return (self.max - self.min) * y + self.min;
                },
            }
        }

        /// Allocate an array and populate with the remaining grid elements.
        /// Caller owns the memory.
        pub fn drain(self: *Self, allocator: std.mem.Allocator) ![]T {
            var out = try allocator.alloc(T, self.itt.remaining);
            var i: usize = 0;
            while (self.next()) |v| : (i += 1) out[i] = v;
            return out;
        }
    };
}

test "grid" {
    var grid: Grid(f64) = .inverse;
    const values = try grid.fill(
        std.testing.allocator,
        0.1,
        10.0,
        10,
    );
    defer std.testing.allocator.free(values);

    try std.testing.expectEqualSlices(f64, &.{
        0.1,
        0.11235955056179775,
        0.1282051282051282,
        0.14925373134328357,
        0.17857142857142852,
        0.22222222222222213,
        0.2941176470588234,
        0.4347826086956519,
        0.8333333333333323,
        9.999999999999858,
    }, values);

    var grid_geom = Grid(f64).geometric;

    const values_geom = try grid_geom.fill(
        std.testing.allocator,
        0.1,
        10.0,
        10,
    );
    defer std.testing.allocator.free(values_geom);

    try std.testing.expectEqualSlices(f64, &.{
        0.1,
        0.1668100537200059,
        0.2782559402207125,
        0.46415888336127814,
        0.7742636826811276,
        1.2915496650148846,
        2.1544346900318847,
        3.593813663804628,
        5.994842503189411,
        10.000000000000004,
    }, values_geom);
}
