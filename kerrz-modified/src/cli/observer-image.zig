const std = @import("std");
const clippy = @import("clippy");
const kerrz = @import("kerrz");

const utils = @import("./utils.zig");

pub const short_description = "Render simple images of black holes.";
pub const description =
    \\Render a simple image of a black hole, simulating a distant observer's
    \\image plane. Each pixel corresponds to a single geodesic, and can be
    \\coloured using a `mapper` to represent some physical quantity.
    \\
;

pub const Args = clippy.Arguments([_]clippy.ArgumentDescriptor{
    .{
        .arg = "--spin spin",
        .argtype = f64,
        .default = "0.998",
        .help = "The black hole spin.",
    },
    .{
        .arg = "--incl inclination",
        .argtype = f64,
        .default = "80",
        .help = "The observer inclination in degrees.",
    },
    .{
        .arg = "--resolution pixels",
        .argtype = usize,
        .default = "8",
        .help = "The number of pixels to trace per step in impact parameter space.",
    },
    .{
        .arg = "--alpha alpha",
        .argtype = f64,
        .default = "40",
        .help = "The range of alpha impact parameters",
    },
    .{
        .arg = "--beta beta",
        .argtype = f64,
        .default = "40",
        .help = "The range of beta impact parameters.",
    },
    .{
        .arg = "--dist radius",
        .argtype = f64,
        .default = "1e7",
        .help = "The observer radial distance in rg.",
    },
    utils.MapperArg("redshift"),
    utils.EmissivityFunction.ArgDescriptor,
    .{
        .arg = "--no-false-image",
        .help = "Whether to also trace the false image or not.",
    },
    .{
        .arg = "--nthreads nthreads",
        .argtype = usize,
        .help = "The number of CPU threads to use.",
    },
    .{
        .arg = "--save filename",
        .help = "Save all geodesic endpoints to a CSV table.",
    },
    .{
        .arg = "--filter what",
        .help = "Filter the geodesics before writing. Only valid with the `--save` option.",
    },
    .{
        .arg = "-o/--output filename",
        .default = "output.pgm",
        .help = "The name of the file to write the image to. Currently only the `pgm` file format is supported, though this can be trivially converted to other  file formats with ImageMagick.",
    },
} ++ utils.DiscParameters.ArgDescriptorList());

const Dual0 = kerrz.DualNumber(f64, 0);
const Dual2 = kerrz.DualNumber(f64, 2);
const Range = kerrz.iterators.RangeIterator(f64);

pub fn run(out: *std.Io.Writer, allocator: std.mem.Allocator, itt: *clippy.ArgumentIterator) !void {
    const args = try Args.initParseAll(itt, .{});

    const selected_mappers = try utils.parseMapperValues(allocator, args.map);
    defer allocator.free(selected_mappers);
    const mapper = selected_mappers[0].mapper;

    const filter = try utils.parseFilter(args.filter);

    // Configure the mappers
    const mappers = try utils.configureMappers(
        allocator,
        selected_mappers,
        false,
    );
    defer allocator.free(mappers);

    // Construct the impact parameter grid
    const num_alpha = args.resolution * 2 * @as(usize, @intFromFloat(args.alpha));
    const num_beta = args.resolution * 2 * @as(usize, @intFromFloat(args.beta));

    var alpha_itt = Range.init(-args.alpha, args.alpha, num_alpha);
    const alpha = try alpha_itt.drain(allocator);
    defer allocator.free(alpha);

    var beta_itt = Range.init(-args.beta, args.beta, num_beta);
    const beta = try beta_itt.drain(allocator);
    defer allocator.free(beta);

    const num_threads = utils.getNumThreads(args.nthreads);

    const ca: CommonArguments = .{
        .out = out,
        .allocator = allocator,
        .x = alpha,
        .y = beta,
        .mappers = mappers,
        .num_threads = num_threads,
    };

    try out.print(
        "Rendering {d} x {d} on {d} thread{s}\n",
        .{ alpha.len, beta.len, num_threads, if (num_threads > 1) "s" else "" },
    );

    if (mapper.needsDerivatives()) {
        const x_obs: kerrz.FourVector(Dual2) = .{
            .t = .zero,
            .r = .promote(args.dist),
            .th = .promote(std.math.degreesToRadians(args.incl)),
            .ph = .zero,
        };
        const em_prof = (try utils.EmissivityFunction.fromArgs(
            Dual2,
            allocator,
            args,
        )).profile;
        defer em_prof.deinit(allocator);
        var result = try ca.run(Dual2, x_obs, em_prof, args, .image);
        defer result.deinit(allocator);
        try ca.postProcess(Dual2, &result, args.save, args.output, filter);
    } else {
        const x_obs: kerrz.FourVector(Dual0) = .{
            .t = .zero,
            .r = .promote(args.dist),
            .th = .promote(std.math.degreesToRadians(args.incl)),
            .ph = .zero,
        };
        const em_prof = (try utils.EmissivityFunction.fromArgs(
            Dual0,
            allocator,
            args,
        )).profile;
        defer em_prof.deinit(allocator);
        var result = try ca.run(Dual0, x_obs, em_prof, args, .image);
        defer result.deinit(allocator);
        try ca.postProcess(Dual0, &result, args.save, args.output, filter);
    }
}

pub const CommonArguments = struct {
    out: *std.Io.Writer,
    allocator: std.mem.Allocator,
    num_threads: usize,
    x: []const f64,
    y: []const f64,
    mappers: []const kerrz.Mapper,
    v_prof: kerrz.orbits.VelocityProfiles = .stationary,
    theta_0: f64 = 0,
    emitter: bool = false,

    pub fn run(
        c: CommonArguments,
        comptime T: type,
        x_obs: kerrz.FourVector(T),
        em_prof: kerrz.EmissivityProfile(T),
        args: anytype,
        comptime what: enum { image, sky },
    ) !kerrz.tools.ImageResult(T) {
        const Tool = kerrz.tools.Image(T);

        const metric: kerrz.KerrMetric(T) = .init(.one, .promote(args.spin));

        const disc = try utils.DiscParameters.parseDefault(
            args,
            .{ .thin_disc = .{ .inner_radius = 0, .outer_radius = 30.0 } },
        );

        var tool: Tool = .{
            .image_opts = .{
                .disc = disc.toAccretionDisc(T, metric),
                .em_prof = em_prof,
                .include_false_image = !args.@"no-false-image",
                .mappers = c.mappers,
                .v_source = c.v_prof,
                .theta_0 = c.theta_0,
                .emitter = c.emitter,
            },
            .metric = metric,
            .x_obs = x_obs,
        };

        var threads = try kerrz.ThreadMap.init(
            c.allocator,
            .{ .num_threads = c.num_threads },
        );
        defer threads.deinit();

        return try switch (what) {
            .image => tool.runImpactParameters(
                c.allocator,
                threads,
                c.x,
                c.y,
                .{
                    .out = c.out,
                    .thread_chunk_size = 1024,
                },
            ),
            .sky => tool.runSkyAngles(
                c.allocator,
                threads,
                c.x,
                c.y,
                .{
                    .out = c.out,
                    .thread_chunk_size = 1024,
                },
            ),
        };
    }

    pub fn postProcess(
        c: CommonArguments,
        comptime T: type,
        result: *kerrz.tools.ImageResult(T),
        output_table: ?[]const u8,
        image_filename: []const u8,
        filter: ?utils.Filters,
    ) !void {
        var max_val: ?f64 = null;
        var min_val: ?f64 = null;

        for (0..result.maxIndex()) |i| {
            const v = result.index(i)[0];
            if (v != 0) {
                if (min_val == null or v < min_val.?) {
                    min_val = v;
                }
                if (max_val == null or v > max_val.?) {
                    max_val = v;
                }
            }
        }

        if (max_val == null) {
            max_val = 0;
        }

        if (min_val == null) {
            min_val = 0;
        }

        const max_value = max_val.?;
        const min_value = min_val.?;

        try c.out.print("Min : {d:.9}\n", .{min_value});
        try c.out.print("Max : {d:.9}\n", .{max_value});

        try c.out.print("Writing image to file to '{s}'\n", .{image_filename});

        try utils.writeImageFile(f64, image_filename, result.values, .{
            .width = c.x.len,
            .height = c.y.len,
            .clip_low = min_value - (max_value - min_value) / 10,
            .clip_high = max_value,
            .stride = c.mappers.len,
        });

        if (output_table) |filename| {
            try c.out.print("Serialising table to '{s}'\n", .{filename});
            const fits = try kerrz.FitsFile.init(c.allocator);
            defer fits.deinit();
            (try fits.addHdu()).* = try result.toFITS(c.allocator, .{
                .filter_intersected = filter != null,
                .axes = .{
                    .x = c.x,
                    .y = c.y,
                },
            });
            try fits.save(.{ .path = filename });
        }
    }
};
