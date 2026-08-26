const std = @import("std");
const clippy = @import("clippy");
const kerrz = @import("kerrz");

const utils = @import("./utils.zig");

const CommonArguments = @import("observer-image.zig").CommonArguments;

pub const short_description = "Render the sky of a point in the spacetime.";
pub const description = short_description;

pub const Args = clippy.Arguments([_]clippy.ArgumentDescriptor{
    .{
        .arg = "--spin spin",
        .argtype = f64,
        .default = "0.998",
        .help = "The black hole spin.",
    },
    .{
        .arg = "--resolution pixels",
        .argtype = usize,
        .default = "8",
        .help = "The number of pixels to trace per 5 degree interval.",
    },
    utils.Position.ArgDescriptorDefaults(10, 45, 0),
    utils.VelocityProfile.ArgDescriptor,
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
        .arg = "--reorient",
        .help = "Reorient the image so that the black hole is always at the center.",
    },
    .{
        .arg = "--half",
        .help = "Calculate only the upper hemisphere of the sky.",
    },
    .{
        .arg = "--emitter",
        .help = "Reverse the time direction from the usual calculation to act as if the point is an emitter instead of an observer. This will change, e.g., the lensing effect of the velocity.",
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

/// How many degrees the resolution parameter refines.
const THETA_SCALE = 40;
/// How many degrees the resolution parameter refines. This is twice the
/// THETA_SCALE scale, as phi has twice the range.
const PHI_SCALE = 2 * THETA_SCALE;

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

    const v_prof = try utils.VelocityProfile.fromArgs(args);

    // Parse the potision
    const pos = try utils.Position.fromString(args.position);

    // Construct the sky angle grid
    const grid_height = args.resolution * THETA_SCALE;
    const adjusted_grid_height = @divFloor(grid_height, @as(usize, if (args.half) 2 else 1));
    const grid_width = args.resolution * PHI_SCALE;

    const theta_upper_limit: f64 = if (args.half) std.math.pi / 2.0 else std.math.pi;

    // Run in reverse so we have the correct orientation for the image.
    var theta_itt = Range.init(theta_upper_limit, 0, adjusted_grid_height);
    const theta = try theta_itt.drain(allocator);
    defer allocator.free(theta);

    var phi_itt = Range.init(0, 2 * std.math.pi, grid_width);
    const phi = try phi_itt.drain(allocator);
    defer allocator.free(phi);

    // Reorientation
    const theta_0: f64 = if (args.reorient)
        pos.theta - std.math.pi / 2.0
    else
        0;

    const num_threads = utils.getNumThreads(args.nthreads);

    const ca: CommonArguments = .{
        .out = out,
        .allocator = allocator,
        .x = phi,
        .y = theta,
        .mappers = mappers,
        .num_threads = num_threads,
        .theta_0 = theta_0,
        .v_prof = v_prof,
        .emitter = args.emitter,
    };

    try out.print(
        "Parsed position as x = (0, {d:.2}, {d:.2}, {d:.2}) for '{s}'\n",
        .{
            pos.r,
            std.math.radiansToDegrees(pos.theta),
            std.math.radiansToDegrees(pos.phi),
            if (args.emitter) "emitter" else "observer",
        },
    );

    try out.print(
        "Rendering {d} x {d} on {d} thread{s}\n",
        .{ phi.len, theta.len, num_threads, if (num_threads > 1) "s" else "" },
    );

    if (mapper.needsDerivatives()) {
        const x_obs: kerrz.FourVector(Dual2) = .{
            .t = .zero,
            .r = .promote(pos.r),
            .th = .promote(pos.theta),
            .ph = .promote(pos.phi),
        };
        const em_prof = (try utils.EmissivityFunction.fromArgs(
            Dual2,
            allocator,
            args,
        )).profile;
        defer em_prof.deinit(allocator);
        var result = try ca.run(Dual2, x_obs, em_prof, args, .sky);
        defer result.deinit(allocator);
        try ca.postProcess(Dual2, &result, args.save, args.output, filter);
    } else {
        const x_obs: kerrz.FourVector(Dual0) = .{
            .t = .zero,
            .r = .promote(pos.r),
            .th = .promote(pos.theta),
            .ph = .promote(pos.phi),
        };
        const em_prof = (try utils.EmissivityFunction.fromArgs(
            Dual0,
            allocator,
            args,
        )).profile;
        defer em_prof.deinit(allocator);
        var result = try ca.run(Dual0, x_obs, em_prof, args, .sky);
        defer result.deinit(allocator);
        try ca.postProcess(Dual0, &result, args.save, args.output, filter);
    }
}
