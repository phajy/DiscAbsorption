const std = @import("std");
const kerrz = @import("kerrz");
const clippy = @import("clippy");

/// The helper used to construct the CLI.
pub fn MapperArg(comptime default: []const u8) clippy.ArgumentDescriptor {
    return .{
        .arg = "--map mapper",
        .default = default,
        .help = "How to colour the image. Use `--help mapper` to see valid values and a brief description of each. This can also be a comma-seperated list of mappers when serialising data.",
    };
}

pub const NamedMapper = struct {
    mapper: kerrz.Mapper,
    name: []const u8,
};

/// Used to parse the `MapperArg`
pub fn parseMapperValues(allocator: std.mem.Allocator, string: []const u8) ![]NamedMapper {
    var list: std.ArrayList(NamedMapper) = .empty;
    defer list.deinit(allocator);

    var itt = std.mem.tokenizeScalar(u8, string, ',');
    while (itt.next()) |token| {
        const name = std.mem.trim(u8, token, "\n \t");

        const mapper = try kerrz.Mapper.fromString(name);
        try list.append(allocator, .{
            .mapper = mapper,
            .name = name,
        });
    }

    return try list.toOwnedSlice(allocator);
}

const DEFAULT_MAPPERS = [_]kerrz.Mapper{
    .{ .value = .mino },
    .{ .value = .time },
    .{ .value = .radius },
    .{ .value = .theta },
    .{ .value = .azimuth },
    .{ .value = .status },
};

pub fn configureMappers(
    allocator: std.mem.Allocator,
    selected_mappers: []NamedMapper,
    include_defaults: bool,
) ![]kerrz.Mapper {
    var mappers: std.ArrayList(kerrz.Mapper) = .empty;
    defer mappers.deinit(allocator);
    for (selected_mappers) |mapper| {
        try mappers.append(allocator, mapper.mapper);
    }
    if (include_defaults) {
        try mappers.appendSlice(allocator, &DEFAULT_MAPPERS);
    }
    return mappers.toOwnedSlice(allocator);
}

/// Construct a compile-time string that lists all of the fields in an
/// enumeration.
pub fn makeList(comptime T: type) []const u8 {
    var s: []const u8 = "";
    switch (@typeInfo(T)) {
        .@"enum" => |info| {
            inline for (info.fields, 0..) |field, i| {
                if (i == info.fields.len - 1) {
                    s = s ++ std.fmt.comptimePrint("and `{s}`.", .{field.name});
                } else {
                    s = s ++ std.fmt.comptimePrint("`{s}`, ", .{field.name});
                }
            }
        },
        else => @compileError("Cannot generate field list for given type"),
    }
    return s;
}

/// Parse a comma-seperated list to a list of enumerated values.
pub fn parseCommaSeperatedEnum(
    allocator: std.mem.Allocator,
    comptime T: type,
    string: []const u8,
) ![]const T {
    var list = std.ArrayList(T).empty;
    defer list.deinit(allocator);

    var itt = std.mem.tokenizeScalar(u8, string, ',');
    while (itt.next()) |token| {
        const value = std.meta.stringToEnum(T, token) orelse
            return error.InvalidProperty;
        try list.append(allocator, value);
    }

    return try list.toOwnedSlice(allocator);
}

fn createOrOverwrite(dir: std.fs.Dir, rel_path: []const u8) !std.fs.File {
    return try dir.createFile(rel_path, .{});
}

/// Write to file given that the data type has a and
///
///     data.writeAll(*std.Io.Writer) !void
///
/// method function.
pub fn writeFile(rel_path: []const u8, data: anytype) !void {
    const dir = std.fs.cwd();
    const f = try createOrOverwrite(dir, rel_path);
    defer f.close();

    var buffer: [512]u8 = undefined;
    var f_writer = f.writer(&buffer);
    const writer = &f_writer.interface;

    try data.writeAll(writer);

    try writer.flush();
}

/// Create or overwrite a file at `rel_path` with the supplied columns of data.
pub fn writeFileColumns(comptime T: type, rel_path: []const u8, columns: []const []const T) !void {
    const len = columns[0].len;
    for (columns) |col| {
        std.debug.assert(col.len == len);
    }

    const dir = std.fs.cwd();
    const f = try createOrOverwrite(dir, rel_path);
    defer f.close();

    var buffer: [512]u8 = undefined;
    var f_writer = f.writer(&buffer);
    const writer = &f_writer.interface;

    for (0..len) |i| {
        for (columns, 0..) |col, j| {
            if (j == 0) {
                // Don't print a comma on the first column.
                try writer.print("{d:.9}", .{col[i]});
            } else {
                try writer.print(", {d:.9}", .{col[i]});
            }
        }
        try writer.writeAll("\n");
    }

    try writer.flush();
}

pub fn writeError(err: anyerror, comptime fmt: []const u8, args: anytype) !void {
    var stderr = std.fs.File.stderr();
    var stderr_writer = stderr.writer(&.{});
    const writer = &stderr_writer.interface;

    try writer.print("{any}", .{err});
    try writer.writeAll(": ");
    try writer.print(fmt, args);
    try writer.writeAll("\n");

    return err;
}

fn WriteImageOptions(comptime T: type) type {
    return struct {
        clip_low: ?T = null,
        clip_high: ?T = null,
        width: usize,
        height: usize,
        stride: usize,
    };
}

fn valueClip(comptime T: type, v: T, min_value: T, max_value: T) u32 {
    return @intFromFloat(255 * @max(0, v - min_value) / (max_value - min_value));
}

pub fn writeImageFile(
    comptime T: type,
    filename: []const u8,
    output: []const T,
    opts: WriteImageOptions(T),
) !void {
    const file_path = "./output.pgm";
    std.fs.cwd().deleteFile(file_path) catch {};

    var f = try std.fs.cwd().createFile(filename, .{});
    defer f.close();

    var buffer: [1024]u8 = undefined;
    var f_writer = f.writer(&buffer);
    const writer = &f_writer.interface;

    // write the image header
    try writer.print("P2\n{d} {d}\n255\n", .{ opts.width, opts.height });
    for (0..opts.height) |reverse_i| {
        const i = opts.height - reverse_i - 1;
        for (0..opts.width) |j| {
            const index = ((i * opts.width) + j) * opts.stride;
            const v = if (output[index] == 0) 0.0 else output[index];
            if (opts.clip_low != null and opts.clip_high != null) {
                const grayscale = valueClip(T, v, opts.clip_low.?, opts.clip_high.?);
                try writer.print("{d} ", .{grayscale});
            } else {
                try writer.print("{d} ", .{@as(u32, @intFromFloat(v))});
            }
        }
        try writer.writeAll("\n");
    }

    try writer.flush();
}

/// A helper for parsing position as input.
pub const Position = struct {
    r: f64,
    theta: f64,
    phi: f64,

    /// The helper used to construct the CLI.
    pub fn ArgDescriptorDefaults(
        comptime r: f64,
        comptime theta: f64,
        comptime phi: f64,
    ) clippy.ArgumentDescriptor {
        const default = std.fmt.comptimePrint("{d},{d},{d}", .{ r, theta, phi });
        return .{
            .arg = "--position pos",
            .default = default,
            .help = "The three-position at which to image the sky. Conventionally should be `r,theta,phi`, comma-seperated values, where the two angular (Boyer-Lindquist) coordinates are in degrees, but may also be specified as `x:X,h:z` to match the ring-like corona.",
        };
    }

    /// Convert a string to a Position.
    pub fn fromString(string: []const u8) !Position {
        var out: [3]f64 = .{ 0, 0, 0 };
        var count: usize = 0;
        var all_cartesian: ?bool = null;
        var itt = std.mem.tokenizeScalar(u8, string, ',');

        while (itt.next()) |token| {
            var key: []const u8 = "";
            var value: []const u8 = token;
            var index = count;
            var is_cartesian: bool = false;
            if (std.mem.indexOfScalar(u8, token, ':')) |split| {
                key = token[0..split];
                value = token[split + 1 ..];
            }

            if (std.mem.eql(u8, key, "x")) {
                is_cartesian = true;
                index = 0;
            }
            if (std.mem.eql(u8, key, "h") or std.mem.eql(u8, key, "z")) {
                is_cartesian = true;
                index = 2;
            }
            if (std.mem.eql(u8, key, "y")) {
                is_cartesian = true;
                index = 1;
            }
            if (std.mem.eql(u8, key, "r")) {
                is_cartesian = false;
                index = 0;
            }
            if (std.mem.eql(u8, key, "th")) {
                is_cartesian = false;
                index = 1;
            }
            if (std.mem.eql(u8, key, "ph")) {
                is_cartesian = false;
                index = 2;
            }

            const v = try std.fmt.parseFloat(f64, value);
            out[index] = v;
            count += 1;

            if (all_cartesian) |c| {
                if (c != is_cartesian)
                    return error.CoordinateMismatch;
            } else {
                all_cartesian = is_cartesian;
            }

            if (count > 3) {
                return error.TooManyPositions;
            }
        }

        if (count < 2) {
            return error.InvalidPosition;
        }

        if (all_cartesian orelse false) {
            const r = @sqrt(out[0] * out[0] + out[1] * out[1] + out[2] * out[2]);
            const projected_r = @sqrt(out[0] * out[0] + out[1] * out[1]);
            const th = std.math.atan2(projected_r, out[2]);
            const ph = std.math.atan2(out[1], out[0]);
            return .{
                .r = r,
                .theta = th,
                .phi = ph,
            };
        } else {
            return .{
                .r = out[0],
                .theta = std.math.degreesToRadians(out[1]),
                .phi = std.math.degreesToRadians(out[2]),
            };
        }
    }
};

/// A helper for parsing initial conditions for null-geodesics.
pub const InitialParameters = struct {
    pub const Parameters = union(enum) {
        impact_parameters: struct {
            alpha: f64,
            beta: f64,
        },
        sky_angles: struct {
            theta: f64,
            phi: f64,
        },
    };

    pub fn ArgDescriptorList() []const clippy.ArgumentDescriptor {
        return &.{
            .{
                .arg = "--impact parameters",
                .display_name = "--impact alpha,beta",
                .help = "Specify the initial conditions of the null-geodesic by a comma seperated pair of impact parameters on the image plane, `alpha,beta`.",
            },
            .{
                .arg = "--angles angles",
                .display_name = "--angles theta,phi",
                .help = "Specify the initial conditions of the null-geodesic by a comma seperated pair of angles on the local sky `theta,phi`. The angles should be given in degrees. If specific units are preferred (i.e. degrees or radians), use `deg` or `rad` as `theta:rad,phi:deg`",
            },
        };
    }

    const InputType = enum {
        none,
        impact_parameters,
        sky_angles,
    };

    const AngularUnits = enum {
        radians,
        degrees,

        /// Convert the value from degrees to the specified unit.
        fn convert(self: AngularUnits, value: f64) f64 {
            return switch (self) {
                .radians => value,
                .degrees => std.math.degreesToRadians(value),
            };
        }
    };

    defaults: ?Parameters = null,

    /// Returns sky-angles in radians.
    pub fn parse(self: InitialParameters, args: anytype) !Parameters {
        var how: InputType = .none;
        var string: []const u8 = "";

        if (args.impact) |impact| {
            how = .impact_parameters;
            string = impact;
        }

        if (args.angles) |angles| {
            if (how != .none) {
                return error.MultipleInitialParameters;
            }
            how = .sky_angles;
            string = angles;
        }

        if (how == .none and self.defaults != null) {
            return self.defaults.?;
        }

        var out: [2]f64 = undefined;
        var units: [2]AngularUnits = undefined;
        var itt = std.mem.tokenizeScalar(u8, string, ',');

        for (&out, &units) |*o, *u| {
            const token = itt.next() orelse return error.TooFewParameters;

            var value_token = token;
            var unit: AngularUnits = .degrees;

            if (std.mem.indexOfScalar(u8, token, ':')) |unit_split| {
                if (how != .sky_angles) {
                    return error.InvalidUnit;
                }

                value_token = token[0..unit_split];

                if (token.len == unit_split + 1) {
                    return error.InvalidUnit;
                }

                const unit_token = token[unit_split + 1 ..];

                if (std.mem.eql(u8, unit_token, "deg")) {
                    unit = .degrees;
                } else if (std.mem.eql(u8, unit_token, "rad")) {
                    unit = .radians;
                } else {
                    return error.InvalidUnit;
                }
            }

            o.* = try std.fmt.parseFloat(f64, value_token);
            u.* = unit;
        }

        if (itt.next() != null) return error.TooManyParameters;

        switch (how) {
            .impact_parameters => return .{ .impact_parameters = .{
                .alpha = out[0],
                .beta = out[1],
            } },
            .sky_angles => return .{ .sky_angles = .{
                .theta = units[0].convert(out[0]),
                .phi = units[1].convert(out[1]),
            } },
            else => unreachable,
        }
    }
};

const TestInitialParametersArgs = struct {
    angles: ?[]const u8 = null,
    impact: ?[]const u8 = null,
};

test "parameter parsing" {
    const ip: InitialParameters = .{};
    {
        const values = try ip.parse(TestInitialParametersArgs{ .angles = "90,45" });
        try std.testing.expectEqual(std.math.degreesToRadians(90), values.sky_angles.theta);
        try std.testing.expectEqual(std.math.degreesToRadians(45), values.sky_angles.phi);
    }
    {
        const values = try ip.parse(TestInitialParametersArgs{ .angles = "1:rad,45" });
        try std.testing.expectEqual(1.0, values.sky_angles.theta);
        try std.testing.expectEqual(std.math.degreesToRadians(45), values.sky_angles.phi);
    }
}

/// Parse the number of threads from the argument, else return whatever the
/// system can support.
pub fn getNumThreads(n: ?usize) usize {
    return n orelse std.Thread.getCpuCount() catch 1;
}

const PropertiesIterator = struct {
    string: []const u8,
    index: usize = 0,

    pub fn init(string: []const u8) PropertiesIterator {
        return .{ .string = string };
    }

    pub fn next(self: *PropertiesIterator) !?struct { key: []const u8, value: []const u8 } {
        if (self.index >= self.string.len) return null;

        const token: []const u8 = b: {
            const comma = std.mem.indexOfScalarPos(u8, self.string, self.index, ',');
            if (comma) |end| {
                const token = self.string[self.index..end];
                self.index = end + 1;
                break :b token;
            }
            const token = self.string[self.index..self.string.len];
            self.index = self.string.len;
            break :b token;
        };

        const kv_split = std.mem.indexOfScalar(u8, token, ':') orelse
            return error.MalformedProperty;

        const key = token[0..kv_split];
        const value = token[kv_split + 1 .. token.len];

        return .{
            .key = key,
            .value = value,
        };
    }
};

/// An argument for parsing corona values.
pub const CoronalParameters = union(enum) {
    lamppost: struct {
        height: f64,
        radial_velocity: f64,
    },
    ring: struct {
        height: f64,
        radius: f64,
    },
    disc: struct {
        height: f64,
        inner_radius: f64,
        outer_radius: f64,
        n_rings: usize,
    },
    umbrella: struct {
        offset_radius: f64,
        inner_opening_angle: f64,
        outer_opening_angle: f64,
        n_rings: usize,
    },

    /// The helper used to construct the CLI.
    pub fn ArgDescriptorList() []const clippy.ArgumentDescriptor {
        return &.{
            .{
                .arg = "--lamppost properties",
                .help = "The lamppost corona, where the properties are height `h` and radial velocity `vr`. If no corona is specified, this is the default with `h:5,vr:0`.",
            },
            .{
                .arg = "--ring-like properties",
                .help = "The ring-like coronal model, where the properties are the height `h` and the radius `x`.",
            },
            .{
                .arg = "--disc-like properties",
                .help = "A disc-like coronal model, where the properties are the height `h`, inner and router radius `rin` and `rout` respectively, and the number of rings to compose the corona out of `nr`.",
            },
            .{
                .arg = "--umbrella properties",
                .help = "An umbrella coronal model. The properies are the offset radius from the black hole, `r`, along with the inner (`inner`) and outer opening angles (`angle`), specified in degrees. If not given, `inner` is set to zero. The number of rings used to compose the umbrella, `nr`, may also be specified.",
            },
        };
    }

    /// Returns the parsed coronal model properties.
    pub fn parse(args: anytype) !CoronalParameters {
        var corona_count: usize = 0;
        if (args.lamppost != null) corona_count += 1;
        if (args.@"ring-like" != null) corona_count += 1;
        if (args.@"disc-like" != null) corona_count += 1;
        if (args.umbrella != null) corona_count += 1;
        if (corona_count > 1) {
            return clippy.ParseError.TooManyArguments;
        }

        if (args.lamppost) |lamppost| {
            var height: ?f64 = null;
            var radial_velocity: ?f64 = null;
            var itt = PropertiesIterator.init(lamppost);
            while (try itt.next()) |kv| {
                if (kv.key.len == 1 and kv.key[0] == 'h') {
                    height = try std.fmt.parseFloat(f64, kv.value);
                } else if (std.mem.eql(u8, kv.key, "vr")) {
                    radial_velocity = try std.fmt.parseFloat(f64, kv.value);
                } else {
                    return error.InvalidProperty;
                }
            }
            if (height == null) return error.MissingProperty;
            return .{ .lamppost = .{
                .height = height.?,
                .radial_velocity = radial_velocity orelse 0,
            } };
        }

        if (args.@"ring-like") |ring_like| {
            var height: ?f64 = null;
            var radius: ?f64 = null;
            var offset: ?f64 = null;
            var angle: ?f64 = null;
            var itt = PropertiesIterator.init(ring_like);
            while (try itt.next()) |kv| {
                if (kv.key.len > 1) {
                    if (std.mem.eql(u8, kv.key, "th") or
                        (std.mem.eql(u8, kv.key, "angle")) or
                        std.mem.eql(u8, kv.key, "theta"))
                    {
                        angle = try std.fmt.parseFloat(f64, kv.value);
                    }
                } else {
                    switch (kv.key[0]) {
                        'h' => {
                            height = try std.fmt.parseFloat(f64, kv.value);
                        },
                        'x' => {
                            radius = try std.fmt.parseFloat(f64, kv.value);
                        },
                        'r' => {
                            offset = try std.fmt.parseFloat(f64, kv.value);
                        },
                        else => {
                            return error.InvalidProperty;
                        },
                    }
                }
            }
            const is_cart = (height != null or radius != null);
            const is_spher = (offset != null or angle != null);
            if (is_cart and is_spher) {
                return error.CoordinateMismatch;
            }

            if (is_spher) {
                if (offset == null) return error.MissingProperty;
                if (angle == null) return error.MissingProperty;
                height = offset.? * @cos(std.math.degreesToRadians(angle.?));
                radius = offset.? * @sin(std.math.degreesToRadians(angle.?));
            }

            if (height == null) return error.MissingProperty;
            if (radius == null) return error.MissingProperty;

            return .{ .ring = .{ .radius = radius.?, .height = height.? } };
        }

        if (args.@"disc-like") |ring_like| {
            var height: ?f64 = null;
            var inner_radius: f64 = 0;
            var outer_radius: ?f64 = null;
            var n_rings: usize = 20;
            var itt = PropertiesIterator.init(ring_like);
            while (try itt.next()) |kv| {
                if (kv.key.len == 1 and kv.key[0] == 'h') {
                    height = try std.fmt.parseFloat(f64, kv.value);
                } else if (std.mem.eql(u8, kv.key, "rin")) {
                    inner_radius = try std.fmt.parseFloat(f64, kv.value);
                } else if (std.mem.eql(u8, kv.key, "rout")) {
                    outer_radius = try std.fmt.parseFloat(f64, kv.value);
                } else if (std.mem.eql(u8, kv.key, "nr")) {
                    n_rings = try std.fmt.parseInt(usize, kv.value, 10);
                } else {
                    return error.InvalidProperty;
                }
            }
            if (height == null) {
                return error.MissingProperty;
            }
            if (outer_radius == null) {
                return error.MissingProperty;
            }
            return .{ .disc = .{
                .inner_radius = if (inner_radius == 0) inner_radius + 1e-3 else inner_radius,
                .outer_radius = outer_radius.?,
                .height = height.?,
                .n_rings = n_rings,
            } };
        }

        if (args.umbrella) |umbrella| {
            var offset_radius: ?f64 = null;
            var inner_opening_angle: f64 = 0;
            var outer_opening_angle: ?f64 = null;
            var n_rings: usize = 20;
            var itt = PropertiesIterator.init(umbrella);
            while (try itt.next()) |kv| {
                if (kv.key.len == 1 and kv.key[0] == 'r') {
                    offset_radius = try std.fmt.parseFloat(f64, kv.value);
                } else if (std.mem.eql(u8, kv.key, "inner")) {
                    inner_opening_angle = try std.fmt.parseFloat(f64, kv.value);
                } else if (std.mem.eql(u8, kv.key, "angle")) {
                    outer_opening_angle = try std.fmt.parseFloat(f64, kv.value);
                } else if (std.mem.eql(u8, kv.key, "nr")) {
                    n_rings = try std.fmt.parseInt(usize, kv.value, 10);
                } else {
                    return error.InvalidProperty;
                }
            }
            if (offset_radius == null) {
                return error.MissingProperty;
            }
            if (outer_opening_angle == null) {
                return error.MissingProperty;
            }
            return .{ .umbrella = .{
                .inner_opening_angle = if (inner_opening_angle == 0)
                    inner_opening_angle + 1e-3
                else
                    inner_opening_angle,
                .outer_opening_angle = outer_opening_angle.?,
                .offset_radius = offset_radius.?,
                .n_rings = n_rings,
            } };
        }

        // the default
        return .{ .lamppost = .{
            .height = 5,
            .radial_velocity = 0,
        } };
    }
};

const EmissivityProfile = kerrz.emissivity.EmissivityProfile;

pub const EmissivityFunction = struct {
    pub fn EmissivityWithParameters(comptime T: type) type {
        return struct {
            profile: kerrz.EmissivityProfile(T),
            params: ?kerrz.CoronalModelOptions(T),

            pub fn deinit(self: *const @This(), allocator: std.mem.Allocator) void {
                self.profile.deinit(allocator);
            }
        };
    }

    pub const ArgDescriptor: clippy.ArgumentDescriptor = .{
        .arg = "--emissivity-profile emissivity",
        .help = "The argument can also be a file to an emissivity table, where the first column are the disc radii and the second the emissivity values.",
    };

    /// Caller must call deinit.
    pub fn fromArgs(
        comptime T: type,
        allocator: std.mem.Allocator,
        args: anytype,
    ) !EmissivityWithParameters(T) {
        if (args.@"emissivity-profile") |ep| {
            return try fromFile(T, allocator, ep);
        }
        return .{
            .profile = .{ .powerlaw = .{} },
            .params = null,
        };
    }

    fn fromFile(
        comptime T: type,
        allocator: std.mem.Allocator,
        path: []const u8,
    ) !EmissivityWithParameters(T) {
        const dir = std.fs.cwd();

        const fits = try kerrz.FitsFile.open(
            allocator,
            .{ .dir = dir, .path = path },
        );
        defer fits.deinit();
        const hdu = &fits.hdus[1];

        const em = try kerrz.emissivity.AxisymmetricEmissivity(T).fromFITS(
            allocator,
            hdu,
        );

        const kzmodel = hdu.getRecord("KZMODEL").?;

        if (std.mem.eql(u8, kzmodel.value.string, "LAMPPOST")) {
            const height = hdu.getRecord("LPHEIGHT").?.value.float;
            const r_vel = hdu.getRecord("LPRVEL").?.value.float;
            const gamma = hdu.getRecord("LPGAMMA").?.value.float;

            return .{
                .profile = .{ .axisymmetric = em },
                .params = .{
                    .lamppost = .{
                        .height = .promote(height),
                        .radial_velocity = .promote(r_vel),
                        .photon_index = .promote(gamma),
                    },
                },
            };
        }

        if (std.mem.eql(u8, kzmodel.value.string, "RING")) {
            const radius = hdu.getRecord("RRADIUS").?.value.float;
            const height = hdu.getRecord("RHEIGHT").?.value.float;
            const gamma = hdu.getRecord("RGAMMA").?.value.float;
            const r_vel = std.meta.stringToEnum(
                kerrz.orbits.VelocityProfiles,
                hdu.getRecord("RVELPROF").?.value.string,
            ).?;

            return .{
                .profile = .{ .axisymmetric = em },
                .params = .{
                    .ring = .{
                        .radius = .promote(radius),
                        .height = .promote(height),
                        .photon_index = .promote(gamma),
                        .velocity = r_vel,
                    },
                },
            };
        }

        unreachable;
    }
};

/// An argument for parsing disc related values.
pub const DiscParameters = union(enum) {
    no_disc: void,
    thin_disc: struct {
        inner_radius: f64,
        outer_radius: f64,
    },
    datum_plane: struct {
        height: f64,
        inner_radius: f64,
        outer_radius: f64,
    },
    sphere: struct {
        radius: f64,
    },
    shakura_sunyaev: struct {
        /// Radiative efficiency.
        reff: ?f64 = null,
        /// Eddington fraction.
        edd: f64,
    },

    /// Convert to the accretion disc geometry wrapper type.
    pub fn toAccretionDisc(
        self: DiscParameters,
        comptime T: type,
        metric: kerrz.KerrMetric(T),
    ) kerrz.AccretionDisc(T) {
        return switch (self) {
            .no_disc => unreachable,
            .thin_disc => |d| .{ .thin_disc = .{
                .inner_radius = .promote(d.inner_radius),
                .outer_radius = .promote(d.outer_radius),
            } },
            .datum_plane => |d| .{ .datum_plane = .{
                .height = .promote(d.height),
                .inner_radius = .promote(d.inner_radius),
                .outer_radius = .promote(d.outer_radius),
            } },
            .sphere => |d| .{ .sphere = .{
                .radius = .promote(d.radius),
            } },
            .shakura_sunyaev => |d| {
                var ssd = kerrz.accretion_discs.ShakuraSunyaev(T){
                    .eddington_fraction = .promote(d.edd),
                    .state = .init(metric),
                };
                if (d.reff) |eff| {
                    ssd.state.radiative_efficiency = .promote(eff);
                }
                return .{ .shakura_sunyaev = ssd };
            },
        };
    }

    /// The helper used to construct the CLI.
    pub fn ArgDescriptorList() []const clippy.ArgumentDescriptor {
        return &.{
            .{
                .arg = "--thin-disc props",
                .help = "A thin disc in the equatorial plane. Properties are `rin` for inner radius and `rout` for the outer radius. If `rin` is not provided, it defaults to the ISCO.",
            },
            .{
                .arg = "--datum-plane props",
                .help = "Use a datum-plane as an accretion disc, i.e. an infinite plane that is parallel to the equatorial plane, but lifted some height `h` above the equator.",
            },
            .{
                .arg = "--sphere-disc props",
                .help = "A sphere centered at the origin. This is a utility disc and not supposed to be representative of any physical model. The only property for this disc is the radius `r`.",
            },
            .{
                .arg = "--shakura-sunyaev props",
                .help = "The Shakura-Sunyaev (1973) accretion disc, with some non-zero scale height controlled by the radiative efficiency (`reff`) and the Eddington accretion fraction (`edd`). `reff` is by default set to the radiative efficiency at the ISCO, and `edd = 0.3`",
            },
        };
    }

    pub fn parseDefault(args: anytype, default: DiscParameters) !DiscParameters {
        const parsed = try parse(args);
        if (parsed == .no_disc) return default;
        return parsed;
    }

    /// Returns the parsed coronal model properties.
    pub fn parse(args: anytype) !DiscParameters {
        var active_count: usize = 0;
        if (args.@"thin-disc" != null) active_count += 1;
        if (args.@"datum-plane" != null) active_count += 1;
        if (args.@"sphere-disc" != null) active_count += 1;
        if (args.@"shakura-sunyaev" != null) active_count += 1;

        if (active_count == 0) return .{ .no_disc = {} };

        if (active_count > 1) {
            return clippy.ParseError.TooManyArguments;
        }

        if (args.@"thin-disc") |thin_disc| {
            var rin: ?f64 = null;
            var rout: ?f64 = null;
            var itt = PropertiesIterator.init(thin_disc);
            while (try itt.next()) |kv| {
                if (std.mem.eql(u8, kv.key, "rin")) {
                    rin = try std.fmt.parseFloat(f64, kv.value);
                } else if (std.mem.eql(u8, kv.key, "rout")) {
                    rout = try std.fmt.parseFloat(f64, kv.value);
                } else {
                    return error.InvalidProperty;
                }
            }
            return .{
                .thin_disc = .{
                    .inner_radius = rin orelse -1,
                    .outer_radius = rout orelse 30,
                },
            };
        }

        if (args.@"datum-plane") |datum_plane| {
            var height: ?f64 = null;
            var rin: ?f64 = null;
            var rout: ?f64 = null;
            var itt = PropertiesIterator.init(datum_plane);
            while (try itt.next()) |kv| {
                if (std.mem.eql(u8, kv.key, "h")) {
                    height = try std.fmt.parseFloat(f64, kv.value);
                } else if (std.mem.eql(u8, kv.key, "rin")) {
                    rin = try std.fmt.parseFloat(f64, kv.value);
                } else if (std.mem.eql(u8, kv.key, "rout")) {
                    rout = try std.fmt.parseFloat(f64, kv.value);
                } else {
                    return error.InvalidProperty;
                }
            }
            return .{
                .datum_plane = .{
                    .height = height orelse 1.0,
                    .inner_radius = rin orelse 0.0,
                    .outer_radius = rout orelse std.math.floatMax(f64),
                },
            };
        }

        if (args.@"sphere-disc") |sphere_disc| {
            var radius: ?f64 = null;
            var itt = PropertiesIterator.init(sphere_disc);
            while (try itt.next()) |kv| {
                if (std.mem.eql(u8, kv.key, "r")) {
                    radius = try std.fmt.parseFloat(f64, kv.value);
                } else {
                    return error.InvalidProperty;
                }
            }
            return .{
                .sphere = .{ .radius = radius orelse 30.0 },
            };
        }

        if (args.@"sphere-disc") |sphere_disc| {
            var radius: ?f64 = null;
            var itt = PropertiesIterator.init(sphere_disc);
            while (try itt.next()) |kv| {
                if (std.mem.eql(u8, kv.key, "r")) {
                    radius = try std.fmt.parseFloat(f64, kv.value);
                } else {
                    return error.InvalidProperty;
                }
            }
            return .{
                .sphere = .{ .radius = radius orelse 30.0 },
            };
        }

        if (args.@"shakura-sunyaev") |ss_disc| {
            var reff: ?f64 = null;
            var edd: ?f64 = null;
            var itt = PropertiesIterator.init(ss_disc);
            while (try itt.next()) |kv| {
                if (std.mem.eql(u8, kv.key, "reff")) {
                    reff = try std.fmt.parseFloat(f64, kv.value);
                } else if (std.mem.eql(u8, kv.key, "edd")) {
                    edd = try std.fmt.parseFloat(f64, kv.value);
                } else {
                    return error.InvalidProperty;
                }
            }
            return .{
                .shakura_sunyaev = .{
                    .edd = edd orelse 0.3,
                    .reff = reff,
                },
            };
        }

        unreachable;
    }
};

pub const VelocityProfile = struct {
    pub const ArgDescriptor: clippy.ArgumentDescriptor = .{
        .arg = "--velocity profile",
        .help = "The velocity profile to use for the source location. May be either `stationary`, `corotate`, `lnr`, where `lnr` is the locally non-rotating frame, or `coplunge`, which also uses the plunging four-velocities.",
    };

    pub fn fromArgs(args: anytype) !kerrz.orbits.VelocityProfiles {
        const velocity = args.velocity orelse return .lnr;

        if (std.mem.eql(u8, "stationary", velocity)) {
            return .stationary;
        }

        if (std.mem.eql(u8, "corotate", velocity)) {
            return .co_rotate;
        }

        if (std.mem.eql(u8, "coplunge", velocity)) {
            return .keplerian_plunging;
        }

        if (std.mem.eql(u8, "lnr", velocity)) {
            return .lnr;
        }

        return error.UnknownVelocityProfile;
    }
};

pub const LineprofileIntegrationArguments = &[_]clippy.ArgumentDescriptor{
    .{
        .arg = "--nrsteps n",
        .argtype = usize,
        .default = "3000",
        .help = "How many radial integration steps to take when integrating the transfer function table.",
    },
    .{
        .arg = "--rstepgrid grid",
        .default = "log10",
        .help = "The grid spacing to use for the transfer function integration grid. This may take the same values as `--rgrid`.",
    },
    .{
        .arg = "--ng n",
        .argtype = usize,
        .default = "1000",
        .help = "The size of the energshift grid into which the lineprofile is calculated.",
    },
    .{
        .arg = "--ngstar n",
        .argtype = usize,
        .default = "2800",
        .help = "The fine grid g_star to integrate over.",
    },
};

/// Read a transfer function table from a file.
pub fn readOrReadAndInterpolateTable(
    comptime T: type,
    out: *std.Io.Writer,
    allocator: std.mem.Allocator,
    table_path: []const u8,
    spin: T.T,
    incl: T.T,
) !kerrz.tools.TransferFunctionTable(T).Result {
    var fits = try kerrz.FitsFile.open(allocator, .{ .path = table_path });
    defer fits.deinit();

    var table = try readOrReadAndInterpolateTableImpl(
        T,
        out,
        allocator,
        fits,
        .promote(spin),
        .promote(incl),
    );
    errdefer table.deinit(allocator);

    // Read in the radii
    const radii = try allocator.alloc(T.T, table.transfer_functions.items.len);
    errdefer allocator.free(radii);

    for (radii, table.transfer_functions.items) |*r, tf| {
        r.* = tf.target_radius;
    }

    return .{
        .radii = radii,
        .table = table,
    };
}

fn readOrReadAndInterpolateTableImpl(
    comptime T: type,
    out: *std.Io.Writer,
    allocator: std.mem.Allocator,
    fits: *kerrz.FitsFile,
    spin: T,
    incl: T,
) !kerrz.tools.TransferFunctionTable(T).Table {
    const table_type = fits.hdus[0].getRecord("KZTYPE").?;
    if (std.mem.eql(u8, table_type.value.string, "ctf")) {
        return try kerrz.transfer_tables.parseSingleFromFITS(
            T,
            allocator,
            fits,
        );
    }
    if (std.mem.eql(u8, table_type.value.string, "ctfgrid")) {
        var tot = try kerrz.transfer_tables.parseFromFITS(
            T,
            allocator,
            fits,
        );
        defer tot.deinit(allocator);

        try out.print("Interpolating from transfer function grid:\n", .{});
        try out.print(
            "  Spin: {d} [limits: {d:.5}, {d:.5}]\n",
            .{ spin.x, tot.spins[0], tot.spins[tot.spins.len - 1] },
        );
        try out.print(
            "  Incl: {d} [limits: {d:.5}, {d:.5}]\n",
            .{ incl.x, tot.observer_incls[0], tot.observer_incls[tot.observer_incls.len - 1] },
        );

        return try tot.interpolateAlloc(allocator, spin, incl);
    }

    return error.UnknownKZTYPE;
}

pub const Filters = enum {
    intersected,
};

pub fn parseFilter(string: ?[]const u8) !?Filters {
    if (string) |s| {
        return std.meta.stringToEnum(Filters, s) orelse
            return error.UnknownFilter;
    }
    return null;
}
