const std = @import("std");
const clippy = @import("clippy");
const kerrz = @import("kerrz");

const utils = @import("./utils.zig");

pub const short_description = "Calculate emissivity profiles for a corona model.";
pub const description = short_description ++ "\n";

pub const Args = clippy.Arguments([_]clippy.ArgumentDescriptor{
    .{
        .arg = "--spin spin",
        .argtype = f64,
        .default = "0.998",
        .help = "The black hole spin.",
    },
    .{
        .arg = "--photon-index photon-index",
        .argtype = f64,
        .default = "2.0",
        .help = "The photon index for the powerlaw intensity of the coronal emissision, `E^(1 - photon-index)`.",
    },
    .{
        .arg = "-o/--output filepath",
        .default = "emissivity.fits",
        .help = "The filepath for the serialised output file.",
    },
    .{
        .arg = "--nphotons nphotons",
        .argtype = usize,
        .default = "3000",
        .help = "The number of photons to trace to calculate the emissivity.",
    },
    .{
        .arg = "--nthreads nthreads",
        .argtype = usize,
        .help = "The number of CPU threads to use.",
    },
    .{
        .arg = "--full-table",
        .help = "Save the full time-dependent emissivity table to the output FITS file. This is not done by default to save diskspace.",
    },
    .{
        .arg = "--nphi n",
        .argtype = usize,
        .default = "500",
        .help = "The number of azimuthal bins to use when averaging or serialising the full emissivity table.",
    },
    .{
        .arg = "--nradii n",
        .argtype = usize,
        .default = "500",
        .help = "The number of radial bins to compute the emissivity profile over.",
    },
    utils.VelocityProfile.ArgDescriptor,
} ++ utils.CoronalParameters.ArgDescriptorList());

const Dual = kerrz.DualNumber(f64, 0);
const Range = kerrz.iterators.RangeIterator(Dual.T);

// A wrapper that adds the `format` and `shouldPrint` functions.
const WriteWrapper = struct {
    eps: []const kerrz.EmissivityPoint(f64),

    pub fn writeAll(
        self: @This(),
        writer: *std.Io.Writer,
    ) std.Io.Writer.Error!void {
        for (self.eps) |ep| {
            if (std.math.isNan(ep.em) or ep.r == 0) {
                continue;
            }
            try writer.print(
                "{d}, {d}, {d}, {d}\n",
                .{ ep.r, ep.em, ep.phi, ep.t },
            );
        }
    }
};

const Tool = kerrz.tools.Emissivity(Dual);

pub fn run(
    out: *std.Io.Writer,
    allocator: std.mem.Allocator,
    itt: *clippy.ArgumentIterator,
) !void {
    const args = try Args.initParseAll(itt, .{});
    const num_threads = utils.getNumThreads(args.nthreads);

    const metric: kerrz.KerrMetric(Dual) = .init(.one, .promote(args.spin));
    const coronal_parameters = try utils.CoronalParameters.parse(args);
    const v_prof = try utils.VelocityProfile.fromArgs(args);

    switch (coronal_parameters) {
        .lamppost => |p| {
            try out.print("Lamppost model h:{d}, vr:{d}\n", .{
                p.height,
                p.radial_velocity,
            });
        },
        .ring => |p| {
            try out.print("Ring model h:{d}, x:{d}\n", .{
                p.height,
                p.radius,
            });
        },
        .disc => |p| {
            try out.print("Disc model h:{d}, rin:{d}, rout:{d}, nr:{d}\n", .{
                p.height,
                p.inner_radius,
                p.outer_radius,
                p.n_rings,
            });
        },
        .umbrella => |p| {
            try out.print("Umbrella model r:{d}, θin:{d}°, θout:{d}°, nr:{d}\n", .{
                p.offset_radius,
                p.inner_opening_angle,
                p.outer_opening_angle,
                p.n_rings,
            });
        },
    }

    const corona: kerrz.CoronalModelOptions(Dual) = switch (coronal_parameters) {
        .lamppost => |p| .{
            .lamppost = .{
                .height = .promote(p.height),
                .radial_velocity = .promote(p.radial_velocity),
                .photon_index = .promote(args.@"photon-index"),
            },
        },
        .ring => |p| .{
            .ring = .{
                .height = .promote(p.height),
                .radius = .promote(p.radius),
                .velocity = v_prof,
                .photon_index = .promote(args.@"photon-index"),
            },
        },
        .disc => |p| .{
            .disc = .{
                .height = .promote(p.height),
                .inner_radius = .promote(p.inner_radius),
                .outer_radius = .promote(p.outer_radius),
                .n_rings = p.n_rings,
                .velocity = v_prof,
                .photon_index = .promote(args.@"photon-index"),
            },
        },
        .umbrella => |p| .{
            .umbrella = .{
                .offset_radius = .promote(p.offset_radius),
                .inner_opening_angle = .promote(
                    std.math.degreesToRadians(p.inner_opening_angle),
                ),
                .outer_opening_angle = .promote(
                    std.math.degreesToRadians(p.outer_opening_angle),
                ),
                .n_rings = p.n_rings,
                .velocity = v_prof,
                .photon_index = .promote(args.@"photon-index"),
            },
        },
    };

    switch (corona) {
        .lamppost => {},
        else => {
            if (args.nphotons < 500_000) {
                try out.writeAll("Warning: extended corona may not converge when `--nphotons` is less than 500,000.\n");
            }
        },
    }

    var threads = try kerrz.ThreadMap.init(
        allocator,
        .{ .num_threads = num_threads },
    );
    defer threads.deinit();

    switch (corona) {
        inline .disc, .umbrella => |c| {
            try out.print(
                "Tracing {d} rings with {d} photons ({d} photons per ring) on {d} thread{s}\n",
                .{
                    c.n_rings,
                    c.n_rings * args.nphotons,
                    args.nphotons,
                    num_threads,
                    if (num_threads > 1) "s" else "",
                },
            );
        },
        else => {
            try out.print(
                "Tracing {d} photons on {d} thread{s}\n",
                .{ args.nphotons, num_threads, if (num_threads > 1) "s" else "" },
            );
        },
    }

    const tool: Tool = .{
        .metric = metric,
        .em_opts = .{
            .corona = corona,
            .num_photons = args.nphotons,
            .num_phi = args.nphi,
            .num_radii = args.nradii,
        },
    };

    var result = try tool.run(
        allocator,
        threads,
        .{
            .out = out,
        },
    );
    defer result.deinit(allocator);

    const emissivity_table = b: switch (result) {
        .single => |single_result| {
            break :b single_result.emissivity_table;
        },
        .sum => |sum_result| {
            try out.print("Total time elapsed: {D}\n", .{
                sum_result.runtime * std.time.ns_per_ms,
            });
            break :b sum_result.emissivity_table;
        },
    };

    const photon_fractions = b: switch (result) {
        .single => |single_result| {
            break :b single_result.photon_fractions;
        },
        .sum => |sum_result| {
            break :b sum_result.photon_fractions;
        },
    };

    try out.writeAll("Photon fractions:\n");
    try out.print("      No status : {d:.5}\n", .{photon_fractions.no_status});
    try out.print("  Event horizon : {d:.5}\n", .{photon_fractions.event_horizon});
    try out.print("       Infinity : {d:.5}\n", .{photon_fractions.infinity});
    try out.print("           Disc : {d:.5}\n", .{photon_fractions.disc});
    try out.print("     Below ISCO : {d:.5}\n", .{photon_fractions.below_isco});
    try out.print("     Above ISCO : {d:.5}\n", .{photon_fractions.above_isco});

    const fits = try kerrz.FitsFile.init(allocator);
    defer fits.deinit();

    var hdu_ptr = try fits.addHdu();
    // Serialise to a FITS file.
    hdu_ptr.* = try emissivity_table.toFITS(allocator);
    try corona.addToHdu(hdu_ptr);
    try photon_fractions.addToHdu(hdu_ptr);
    try metric.addToHdu(hdu_ptr);

    // If the result is a sum of many rings, save each ring to the HDU also:
    switch (result) {
        .sum => |sum_result| {
            const hdu_ring_info_ptr = try fits.addHdu();
            hdu_ring_info_ptr.* = try kerrz.emissivity.makeRingInfoHdu(allocator);
            // Refresh the pointer, as the address may have changed:
            hdu_ptr = &fits.hdus[1];

            for (sum_result.results) |res| {
                try kerrz.emissivity.writeRingToFITS(
                    Dual,
                    fits.arena.allocator(),
                    res.corona.ring,
                    res.photon_fractions,
                    res.emissivity_table,
                    hdu_ptr,
                    hdu_ring_info_ptr,
                );
            }
        },
        else => {},
    }

    if (args.@"full-table") {
        try out.print("Serialing full time-dependent emissivity table\n", .{});
        switch (result) {
            .single => |single_result| {
                switch (corona) {
                    .ring => {
                        const alloc = fits.arena.allocator();

                        var time_dependent = try single_result.cache.rebinTimeDependent(allocator, .{
                            .n_phi = args.nphi,
                            .n_radii = args.nradii,
                        });
                        defer time_dependent.deinit(allocator);
                        time_dependent.interpolateMissing();

                        const table_hdu = try time_dependent.toFITS(alloc);
                        // Do not need to defer deinit since using the arena
                        // allocator of the FITS file.
                        (try fits.addHdu()).* = table_hdu;
                    },
                    .lamppost => {},
                    else => unreachable,
                }
            },
            .sum => {},
        }
    }

    const output_file = args.output;
    try out.print(
        "Writing emissivity data to '{s}'\n",
        .{output_file},
    );
    try fits.save(.{ .path = output_file });
}
