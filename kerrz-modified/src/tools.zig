const std = @import("std");
const zfits = @import("zfits");
const geometry = @import("geometry.zig");
const geodesic = @import("geodesic.zig");
const transfer_tables = @import("transfer-tables.zig");
const iterators = @import("iterators.zig");
const emissivity = @import("emissivity.zig");
const interpolations = @import("interpolations.zig");
const accretion_discs = @import("accretion-discs.zig");
const orbits = @import("orbits.zig");
const utils = @import("utils.zig");

const KerrMetric = geometry.KerrMetric;
const FourVector = geometry.FourVector;
const NullGeodesic = geodesic.NullGeodesic;
const Heuristic = transfer_tables.Heuristic;
const Refine = transfer_tables.Refine;
const Grid = iterators.Grid;
const EmissivityProfile = emissivity.EmissivityProfile;
const ThreadMap = @import("threads.zig").ThreadMap;
const AccretionDisc = accretion_discs.AccretionDisc;
const Mapper = @import("Mapper.zig");
const Matrix = @import("matrix.zig").Matrix;

const logger = std.log.scoped(.tools);

/// Options that control how a tool is run.
pub const RunToolOptions = struct {
    /// The output stream. If none given, will not print any output.
    out: ?*std.Io.Writer = null,
    /// Should a progress bar be displayed. The progress bar is written with
    /// the `out` writer.
    show_progress: bool = true,
    /// How the threads should be chunked.
    thread_chunk_size: usize = 4096,
    /// Whether to add a new line at the end of the time elapsed printing.
    new_line: bool = true,
};

fn _print(opts: RunToolOptions, comptime fmt: []const u8, args: anytype) !void {
    if (opts.out) |out| {
        try out.print(fmt, args);
    }
}

/// A utility method for running threaded workloads with optional progress
/// indicators and runtime reporting.
pub fn runThreads(
    threads: *ThreadMap,
    output: anytype,
    ctx: anytype,
    func: anytype,
    options: RunToolOptions,
) !void {
    const time_now = std.time.milliTimestamp();
    // Set the thread pool going and block until it is done.
    try threads.map(@typeInfo(@TypeOf(output)).pointer.child, output, ctx, func, .{
        .chunk_size = options.thread_chunk_size,
    });
    if (options.show_progress) {
        try threads.blockWithProgress(.{});
    } else {
        threads.blockUntilDone();
    }
    const duration = std.time.milliTimestamp() - time_now;

    try _print(options, "Time elapsed: {D}{s}", .{
        duration * std.time.ns_per_ms,
        if (options.new_line) "\n" else "",
    });
}

/// This is the `tf` command in the CLI.
pub fn TransferFunction(comptime T: type) type {
    return struct {
        pub const Table = transfer_tables.CunninghamTransferFunctionTable(T);
        pub const Options = Table.Options;
        const Self = @This();

        metric: KerrMetric(T),
        x_obs: FourVector(T),
        disc: accretion_discs.AccretionDisc(T) = .equatorial_plane,

        r_target: T.T,

        opts: Options,

        pub const Result = struct {
            ctf: Table.TransferFunction,

            pub fn deinit(self: *Result, allocator: std.mem.Allocator) void {
                self.ctf.deinit(allocator);
                self.* = undefined;
            }
        };

        pub fn run(self: *const Self, allocator: std.mem.Allocator) !Result {
            const table = Table.init(self.metric, self.x_obs, self.disc);
            return .{
                .ctf = try table.calculateRadius(allocator, self.r_target, self.opts),
            };
        }
    };
}

/// This is the `table` command in the CLI.
pub fn TransferFunctionTable(comptime T: type) type {
    return struct {
        pub const Table = transfer_tables.CunninghamTransferFunctionTable(T);
        pub const TFOptions = Table.Options;
        const Self = @This();

        pub const Options = struct {
            num_radii: usize = 100,
            r_in: ?T.T = null,
            r_out: T.T = 1000,
            r_grid: Grid(T.T) = .log10,
        };

        metric: KerrMetric(T),
        x_obs: FourVector(T),
        disc: AccretionDisc(T) = .equatorial_plane,

        table_opts: Options = .{},
        tf_opts: TFOptions,

        pub const Result = struct {
            radii: []const T.T,
            table: Table,

            pub fn deinit(self: *Result, allocator: std.mem.Allocator) void {
                allocator.free(self.radii);
                self.table.deinit(allocator);
                self.* = undefined;
            }
        };

        /// Load a single table from a file instead of computing it.
        ///
        /// TODO: maybe remove this and prefer `transfer_tables.readFromFile`
        /// for all file operations?
        pub fn readFromFile(
            allocator: std.mem.Allocator,
            table_path: []const u8,
        ) !Result {
            var table = try transfer_tables.readSingleFromFile(T, allocator, table_path);
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

        pub fn run(
            self: *const Self,
            allocator: std.mem.Allocator,
            threads: *ThreadMap,
            run_opts: RunToolOptions,
        ) !Result {
            var table = Table.init(self.metric, self.x_obs, self.disc);
            errdefer table.deinit(allocator);

            // Allocate caches
            const radii = try self.table_opts.r_grid.fill(
                allocator,
                self.table_opts.r_in orelse self.metric.isco.x,
                self.table_opts.r_out,
                self.table_opts.num_radii,
            );
            errdefer allocator.free(radii);

            const transfer_functions = try allocator.alloc(Table.CTF, self.table_opts.num_radii);
            errdefer allocator.free(transfer_functions);

            const thread_ctx: ThreadContext = .{
                .table = &table,
                .radii = radii,
                .tf_opts = self.tf_opts,
                .allocator = allocator,
            };

            if (run_opts.out) |out| {
                const n_threads = threads.opts.num_threads;
                try out.print(
                    "Calculating {d} transfer functions on {d} thread{s}\n",
                    .{
                        radii.len,
                        n_threads,
                        if (n_threads == 1) "" else "s",
                    },
                );
            }

            try runThreads(
                threads,
                transfer_functions,
                thread_ctx,
                ThreadContext.work,
                run_opts,
            );

            table.transfer_functions = .fromOwnedSlice(transfer_functions);

            return .{
                .radii = radii,
                .table = table,
            };
        }

        const ThreadContext = struct {
            table: *const Table,
            radii: []const T.T,
            tf_opts: TFOptions,
            allocator: std.mem.Allocator,

            fn work(
                self: ThreadContext,
                out_ctf: *Table.CTF,
                index: usize,
                id: ThreadMap.ThreadId,
            ) void {
                const target_radius = self.radii[index];
                const ctf = self.table.calculateRadius(
                    self.allocator,
                    target_radius,
                    self.tf_opts,
                ) catch |err| {
                    std.debug.print(
                        "Error on thread {d}: {t} on index {d} for radius {d}",
                        .{
                            id.@"0",
                            err,
                            index,
                            target_radius,
                        },
                    );
                    return;
                };
                out_ctf.* = ctf;
            }
        };
    };
}

/// This is the `lineprof` command in the CLI.
pub fn Lineprofile(comptime T: type) type {
    return struct {
        pub const Table = transfer_tables.CunninghamTransferFunctionTable(T);
        pub const Integrator = transfer_tables.TableIntegrator(T);
        pub const IntegrationOptions = Integrator.IntegrationOptions;
        const Self = @This();

        /// The transfer function table to integrate.
        table: Table,

        num_g: usize = 1000,
        num_fine_g: usize = 2800,

        lp_opts: IntegrationOptions,

        pub const Result = struct {
            /// The lineprofile itself
            profile: Profile,
            /// The time spent integrating the lineprofile.
            integration_time: i64,

            pub fn deinit(self: *Result, allocator: std.mem.Allocator) void {
                self.profile.deinit(allocator);
                self.* = undefined;
            }
        };

        fn checkField(
            self: *const Self,
            field: transfer_tables.SelectedField,
            out: ?*std.Io.Writer,
        ) !void {
            if (!self.table.hasField(field)) {
                if (out) |writer| {
                    try writer.print(
                        "Missing field in the transfer table: '{s}'. The table" ++
                            " must be recomputed with the correct fields.\n",
                        .{@tagName(field)},
                    );
                }
                return error.MissingField;
            }
        }

        /// Run the lineprofile intergration but using a pre-calculated table.
        pub fn run(
            self: *const Self,
            allocator: std.mem.Allocator,
            run_opts: RunToolOptions,
        ) !Result {
            // Check the table has the necessary information.
            try self.checkField(.f, run_opts.out);

            // Initialise an energy grid.
            var g_itt = iterators.RangeIterator(T.T).init(0.0, 2.0, self.num_g);
            const g_grid = try g_itt.drain(allocator);
            errdefer allocator.free(g_grid);

            // Allocate a flux grid and zero it.
            const flux_grid = try allocator.alloc(T.T, self.num_g);
            errdefer allocator.free(flux_grid);

            var integrator = try self.table.integrator(
                allocator,
                g_grid,
                self.num_fine_g,
            );
            defer integrator.deinit(allocator);

            const time_now = std.time.milliTimestamp();
            integrator.lineprofile(flux_grid, self.lp_opts);
            const duration = std.time.milliTimestamp() - time_now;

            return .{
                .profile = .{
                    .flux = flux_grid,
                    .g_grid = g_grid,
                },
                .integration_time = duration,
            };
        }

        /// The integrated lineprofile or spectral profile.
        pub const Profile = struct {
            /// The energy-related grid.
            g_grid: []const T.T,
            /// The flux-related grid.
            flux: []const T.T,

            pub fn deinit(self: Profile, allocator: std.mem.Allocator) void {
                allocator.free(self.g_grid);
                allocator.free(self.flux);
            }
        };
    };
}

/// This is the `impulse` command in the CLI.
pub fn ImpulseResponse(comptime T: type) type {
    return struct {
        pub const Table = transfer_tables.CunninghamTransferFunctionTable(T);
        pub const Integrator = transfer_tables.TableIntegrator(T);
        pub const IntegrationOptions = Integrator.IntegrationOptions;
        const Self = @This();

        /// The transfer function table to integrate.
        table: Table,

        num_g: usize = 1000,
        num_fine_g: usize = 2800,

        /// The time grid.
        num_t: usize = 1000,
        t_min: T.T = 1.0,
        t_max: T.T = 1000,
        t_grid: Grid(T.T) = .log10,

        lp_opts: IntegrationOptions,

        pub const Result = struct {
            g_grid: []const T.T,
            t_grid: []const T.T,
            impulse_reponse: Matrix(T.T),
            integration_time: i64,

            pub fn deinit(self: *Result, allocator: std.mem.Allocator) void {
                allocator.free(self.g_grid);
                allocator.free(self.t_grid);
                self.impulse_reponse.deinit(allocator);
                self.* = undefined;
            }
        };

        pub fn run(
            self: *const Self,
            allocator: std.mem.Allocator,
            run_opts: RunToolOptions,
        ) !Result {
            _ = run_opts;

            // Initialise an energy grid.
            var g_itt = iterators.RangeIterator(T.T).init(0.0, 2.0, self.num_g);
            const g_grid = try g_itt.drain(allocator);
            errdefer allocator.free(g_grid);

            // Initialise a time grid.
            const t_grid = try self.t_grid.fill(
                allocator,
                self.t_max,
                self.t_min,
                self.num_t,
            );
            errdefer allocator.free(t_grid);
            std.mem.reverse(T.T, t_grid);

            // Allocate a flux matrix and zero it.
            var flux = try Matrix(T.T).init(allocator, t_grid.len, g_grid.len);
            errdefer flux.deinit(allocator);

            var integrator = try self.table.integrator(
                allocator,
                g_grid,
                self.num_fine_g,
            );
            defer integrator.deinit(allocator);

            const time_now = std.time.milliTimestamp();
            integrator.impulseResponse(t_grid, flux, self.lp_opts);
            const duration = std.time.milliTimestamp() - time_now;

            return .{
                .g_grid = g_grid,
                .t_grid = t_grid,
                .impulse_reponse = flux,
                .integration_time = duration,
            };
        }
    };
}

pub fn ImageResult(comptime T: type) type {
    return struct {
        const Self = @This();
        /// This holds all of the actual result values contiguously. To
        /// iterate over them, use the `index` function.
        values: []T.T,

        /// The status codes for each geodesic.
        statuses: []geodesic.Status,

        /// Whether the original `x` and `y` axes refer to impact paramters or
        /// not.
        is_impact_parameters: bool,

        /// These are the selected mapper names.
        mappers: []const Mapper,

        metric: KerrMetric(T),
        x_obs: FourVector(T),

        /// Get the total number of pixels in the image. It is to be used
        /// for iterating over all traced values.
        pub fn maxIndex(self: *const Self) usize {
            const size = self.mappers.len;
            return @divExact(self.values.len, size);
        }

        /// Get the values associated with the `i`th pixel of the image.
        pub fn index(self: *Self, i: usize) []T.T {
            const size = self.mappers.len;
            const offset = size * i;
            return self.values[offset .. offset + size];
        }

        /// Get the status associated with the `i`th pixel of the image.
        pub fn getStatus(self: *Self, i: usize) geodesic.Status {
            return self.statuses[i];
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.values);
            allocator.free(self.mappers);
            allocator.free(self.statuses);
            self.* = undefined;
        }

        pub const SerialisationOptions = struct {
            /// Only write those pixel that intersected the accretion disc.
            filter_intersected: bool = false,

            /// If not null, will also include the `x` and `y` axis as a column in the
            /// serialised fits file.
            axes: ?struct {
                x: []const T.T,
                y: []const T.T,
            } = null,
        };

        /// Serialise the image data to a FITS binary table HDU.
        pub fn toFITS(self: *Self, allocator: std.mem.Allocator, opts: SerialisationOptions) !zfits.Hdu {
            var hdu = zfits.Hdu.init(allocator, .{ .binary_table = .empty });
            errdefer hdu.deinit();

            // Setup basic information:
            try hdu.setName(
                "IMG_DATA",
                "This is data from a ray-traced image.",
            );
            try utils.addKerrzFITSInfo(&hdu);
            try utils.addObserverInformation(&hdu, self.x_obs);
            try self.metric.addToHdu(&hdu);
            // TODO: serialise information about the accretion disc.

            const table = &hdu.data.binary_table;

            var column_offset: usize = 0;

            // Setup the columns:
            if (opts.axes) |axes| {
                std.debug.assert(self.statuses.len == (axes.x.len * axes.y.len));
                try table.appendColumn(.{
                    .label = if (self.is_impact_parameters) "ax_alpha" else "ax_phi",
                    .units = if (self.is_impact_parameters) "tg" else "deg",
                    .comment = "The x-axis of the image",
                    .col_type = .float_32,
                });
                try table.appendColumn(.{
                    .label = if (self.is_impact_parameters) "ax_beta" else "ax_theta",
                    .units = if (self.is_impact_parameters) "tg" else "deg",
                    .comment = "The y-axis of the image.",
                    .col_type = .float_32,
                });
                column_offset = 2;
            }

            for (self.mappers) |mapper| {
                if (mapper.value == .time) {
                    logger.warn(
                        "Serialising `time` can lose precision. Consider using `delta_time` instead.",
                        .{},
                    );
                }

                try table.appendColumn(.{
                    .label = @tagName(mapper.value),
                    .units = mapper.getUnits(),
                    .comment = mapper.getShortDescriptor(),
                    .col_type = .float_32,
                });
            }

            // Write all of the values.
            const N = self.maxIndex();
            for (0..N) |pixel_index| {
                if (opts.filter_intersected) {
                    if (self.getStatus(pixel_index) != .intersected_disc) {
                        continue;
                    }
                }

                const pixel_values = self.index(pixel_index);

                const row = try table.addRow();

                if (opts.axes) |axes| {
                    const y_index: usize = @divFloor(pixel_index, axes.x.len);
                    const x_index: usize = @rem(pixel_index, axes.x.len);

                    var x_val = axes.x[x_index];
                    if (!self.is_impact_parameters) {
                        x_val = std.math.radiansToDegrees(x_val);
                    }
                    row.cols[0].one.float_32 = @floatCast(x_val);

                    var y_val = axes.y[y_index];
                    if (!self.is_impact_parameters) {
                        y_val = std.math.radiansToDegrees(y_val);
                    }
                    row.cols[1].one.float_32 = @floatCast(y_val);
                }

                for (column_offset.., pixel_values) |col_index, val| {
                    row.cols[col_index].one.float_32 = @floatCast(val);
                }
            }

            return hdu;
        }
    };
}

/// This is the `image` command in the CLI.
pub fn Image(comptime T: type) type {
    return struct {
        pub const Disc = AccretionDisc(T);
        const Self = @This();

        pub const Options = struct {
            include_false_image: bool = true,
            theta_0: T.T = 0,
            disc: Disc = .{ .thin_disc = .{} },
            em_prof: EmissivityProfile(T) = .{ .powerlaw = .{} },
            v_source: orbits.VelocityProfiles = .lnr,
            mappers: []const Mapper = &.{.{ .value = .redshift }},
            /// Whether this image calculation should be considered an
            /// 'emitter' or an 'observer', which is the default if this flag
            /// is false.
            emitter: bool = false,
        };

        metric: KerrMetric(T),
        x_obs: FourVector(T),

        image_opts: Options,

        pub const Result = ImageResult(T);

        /// Run the tool with an impact parameter grid.
        pub fn runSkyAngles(
            self: *const Self,
            allocator: std.mem.Allocator,
            threads: *ThreadMap,
            phi: []const T.T,
            theta: []const T.T,
            run_opts: RunToolOptions,
        ) !Result {
            return try self.runImpl(
                allocator,
                threads,
                phi,
                theta,
                run_opts,
                .sky_angles,
            );
        }

        /// Run the tool with an impact parameter grid.
        pub fn runImpactParameters(
            self: *const Self,
            allocator: std.mem.Allocator,
            threads: *ThreadMap,
            alpha: []const T.T,
            beta: []const T.T,
            run_opts: RunToolOptions,
        ) !Result {
            return try self.runImpl(
                allocator,
                threads,
                alpha,
                beta,
                run_opts,
                .impact_parameters,
            );
        }

        fn runImpl(
            self: *const Self,
            allocator: std.mem.Allocator,
            threads: *ThreadMap,
            x: []const T.T,
            y: []const T.T,
            run_opts: RunToolOptions,
            comptime what: enum { impact_parameters, sky_angles },
        ) !Result {
            const image_size = x.len * y.len;
            const num_mappers = self.image_opts.mappers.len;

            const values = try allocator.alloc(
                T.T,
                image_size * num_mappers,
            );
            errdefer allocator.free(values);

            const statuses = try allocator.alloc(geodesic.Status, image_size);
            errdefer allocator.free(statuses);

            const ts = self.metric.tangentSpace(self.x_obs);
            var v_source = self.image_opts.v_source.fourVector(
                T,
                self.metric,
                ts,
            );

            if (!self.image_opts.emitter) {
                // Negate the `phi` coordinate as for images we are backwards
                // ray-tracing for an observer image.
                v_source.ph = v_source.ph.neg();
                // Normalise again
                v_source = ts.constrainVector(v_source, 1.0);
            }

            const thread_ctx: ThreadContext = .{
                .metric = self.metric,
                .ts = ts,
                .image_opts = self.image_opts,
                .values = values,
                .statuses = statuses,
                .x = x,
                .y = y,
                .v_source = v_source,
                .max_winding = if (self.image_opts.include_false_image) 1 else 0,
            };

            try runThreads(
                threads,
                values[0..image_size],
                thread_ctx,
                switch (what) {
                    .impact_parameters => ThreadContext.workImpactParameters,
                    .sky_angles => ThreadContext.workSkyAngles,
                },
                run_opts,
            );

            return .{
                .values = values,
                .statuses = statuses,
                .mappers = try allocator.dupe(Mapper, self.image_opts.mappers),
                .metric = self.metric,
                .x_obs = self.x_obs,
                .is_impact_parameters = what == .impact_parameters,
            };
        }

        const ThreadContext = struct {
            metric: KerrMetric(T),
            ts: KerrMetric(T).TangentSpace,

            v_source: FourVector(T),

            max_winding: usize,
            image_opts: Options,

            x: []const T.T,
            y: []const T.T,
            values: []T.T,
            statuses: []geodesic.Status,

            fn trace(
                ctx: ThreadContext,
                geod: NullGeodesic(T),
                index: usize,
            ) void {
                const result = ctx.image_opts.disc.trace(
                    ctx.metric,
                    geod,
                    .{ .winding = ctx.max_winding },
                );

                var mapper_ctx = Mapper.Context(T){
                    .geometry = ctx.metric,
                    .v_source = ctx.v_source,
                    .em = ctx.image_opts.em_prof,
                };

                ctx.statuses[index] = result.status;

                // Apply each of the mappers.
                const N = ctx.image_opts.mappers.len;
                for (0..N) |i| {
                    const mapper = ctx.image_opts.mappers[i];
                    ctx.values[index * N + i] = mapper.calculateAlt(
                        T,
                        &mapper_ctx,
                        geod,
                        result,
                    );
                }
            }

            fn workSkyAngles(ctx: ThreadContext, _: *T.T, index: usize, _: ThreadMap.ThreadId) void {
                const row: usize = @divFloor(index, ctx.x.len);
                const col: usize = @rem(index, ctx.x.len);

                // TODO: is this the right way round? Also check that the
                // output write is correct.
                const phi = T.promote(ctx.x[col]);
                const theta = T.promote(ctx.y[row]);

                const geod = NullGeodesic(T).fromSkyAnglesRotatedTangentSpace(
                    ctx.metric,
                    ctx.ts,
                    ctx.v_source,
                    if (T.N != 0) theta.diff(1) else theta,
                    if (T.N != 0) phi.diff(0) else phi,
                    .{ .theta_0 = .promote(ctx.image_opts.theta_0) },
                );

                ctx.trace(geod, index);
            }

            fn workImpactParameters(ctx: ThreadContext, _: *T.T, index: usize, _: ThreadMap.ThreadId) void {
                const row: usize = @divFloor(index, ctx.x.len);
                const col: usize = @rem(index, ctx.x.len);

                // TODO: is this the right way round? Also check that the
                // output write is correct.
                const alpha = T.promote(ctx.x[col]);
                const beta = T.promote(ctx.y[row]);

                const geod = NullGeodesic(T).fromImpactParameters(
                    ctx.metric,
                    ctx.ts.x,
                    if (T.N != 0) alpha.diff(0) else alpha,
                    if (T.N != 0) beta.diff(1) else beta,
                );

                ctx.trace(geod, index);
            }
        };
    };
}

// This is the `emissivity` command.
pub fn Emissivity(comptime T: type) type {
    return struct {
        const Self = @This();

        pub const Options = struct {
            /// The maximum number of photons to trace.
            num_photons: usize = 3000,
            /// The number of radial bins.
            num_radii: usize = 500,
            /// The spacing of the radial bins.
            rgrid: Grid(T.T) = .log10,
            /// The number of phi bins.
            num_phi: usize = 1000,
            /// The maximum radius when rebinning the emissivity profile.
            max_radius: T.T = 1e4,
            /// The coronal parameters themselves.
            corona: emissivity.CoronalModelOptions(T),
        };

        metric: KerrMetric(T),
        em_opts: Options,

        pub const SingleResult = struct {
            /// The computation runtime in milliseconds.
            runtime: i64,
            /// The emissivity cache used in the computation.
            cache: emissivity.EmissivityCache(T),
            /// The (rebinned) emissivity data. Equivalently the
            /// azimuthally-averaged emissivity.
            emissivity_table: emissivity.AxisymmetricEmissivity(T),
            /// The instantiated coronal model.
            corona: emissivity.CoronalModel(T),
            /// The associated photon fractions.
            photon_fractions: emissivity.PhotonFractions(T.T),

            pub fn deinit(self: *SingleResult, allocator: std.mem.Allocator) void {
                self.emissivity_table.deinit(allocator);
                self.cache.deinit(allocator);
                self.corona.deinit(allocator);
                self.* = undefined;
            }
        };

        pub const SummedResult = struct {
            /// The computation runtime in milliseconds.
            runtime: i64,
            /// All results that are summed together.
            results: []SingleResult,
            /// The instantiated coronal model.
            corona: emissivity.CoronalModel(T),
            /// The (rebinned) and summed emissivity data across all rings.
            /// Equivalently the azimuthally-averaged emissivity.
            emissivity_table: emissivity.AxisymmetricEmissivity(T),
            /// The weighted averages of the photon fractions.
            photon_fractions: emissivity.PhotonFractions(T.T),

            pub fn deinit(self: *SummedResult, allocator: std.mem.Allocator) void {
                for (self.results) |*r| r.deinit(allocator);
                allocator.free(self.results);
                self.corona.deinit(allocator);
                self.emissivity_table.deinit(allocator);
                self.* = undefined;
            }
        };

        pub const Result = union(enum) {
            single: SingleResult,
            sum: SummedResult,

            pub fn deinit(self: *Result, allocator: std.mem.Allocator) void {
                switch (self.*) {
                    inline else => |*i| i.deinit(allocator),
                }
            }
        };

        pub fn run(
            self: *const Self,
            allocator: std.mem.Allocator,
            threads: *ThreadMap,
            run_opts: RunToolOptions,
        ) !Result {
            switch (self.em_opts.corona) {
                .lamppost => {
                    const result = try self.runLamppost(
                        allocator,
                        .{ .lamppost = .init(self.metric, self.em_opts.corona.lamppost) },
                        threads,
                        run_opts,
                    );
                    return .{ .single = result };
                },
                .ring => {
                    const result = try self.runRing(
                        allocator,
                        .{ .ring = .init(self.metric, self.em_opts.corona.ring) },
                        threads,
                        run_opts,
                    );
                    return .{ .single = result };
                },
                .disc, .umbrella => {
                    const model = try self.em_opts.corona.toCoronalModel(
                        allocator,
                        self.metric,
                    );
                    errdefer model.deinit(allocator);
                    const summed_result = try self.runRingSum(allocator, model, threads, run_opts);
                    return .{ .sum = summed_result };
                },
            }
        }

        fn runRingSum(
            self: *const Self,
            allocator: std.mem.Allocator,
            corona: emissivity.CoronalModel(T),
            threads: *ThreadMap,
            run_opts: RunToolOptions,
        ) !SummedResult {
            const A = T.Algebra;
            std.debug.assert(switch (corona) {
                .umbrella, .disc => true,
                else => false,
            });

            const rings = switch (corona) {
                inline .disc, .umbrella => |c| c.rings,
                else => unreachable,
            };

            var ring_results = try std.ArrayList(SingleResult).initCapacity(
                allocator,
                rings.len,
            );
            defer ring_results.deinit(allocator);
            // Catch for if things go wrong:
            errdefer for (ring_results.items) |*result| result.deinit(allocator);

            var ring_run_opts = run_opts;
            ring_run_opts.new_line = false;

            var total_runtime: i64 = 0;
            for (0.., rings) |i, ring| {
                if (run_opts.out) |out| {
                    switch (corona) {
                        .disc => {
                            try out.print(
                                "Tracing ring {d} of {d}: r={d:.2}rg\n",
                                .{ i + 1, rings.len, ring.getRadius().x },
                            );
                        },
                        .umbrella => {
                            try out.print(
                                "Tracing ring {d} of {d}: θ={d:.2}°\n",
                                .{ i + 1, rings.len, std.math.radiansToDegrees(ring.x.th.x) },
                            );
                        },
                        else => unreachable,
                    }
                }

                const result = ring_results.addOneAssumeCapacity();
                result.* = try self.runRing(
                    allocator,
                    .{ .ring = ring },
                    threads,
                    ring_run_opts,
                );

                // Clearup the last two lines.
                if (run_opts.out) |out| {
                    if (i == rings.len - 1) {
                        try out.print("\x1b[{d}{c}", .{ 0, 'G' });
                        try out.print("\x1b[{d}{c}", .{ 2, 'K' });
                    }
                    try out.print("\x1b[{d}{c}", .{ 1, 'A' });
                    try out.print("\x1b[{d}{c}", .{ 0, 'G' });
                    try out.print("\x1b[{d}{c}", .{ 2, 'K' });
                }

                total_runtime += result.runtime;
            }

            // Clearup the terminal ready for printing the next line.
            if (run_opts.out) |out| {
                try out.print("\x1b[{d}{c}", .{ 0, 'G' });
                try out.print("\x1b[{d}{c}", .{ 2, 'K' });
            }

            // Sum the azimuthally-averaged emissivity.
            const radii = try allocator.dupe(
                T,
                ring_results.items[0].emissivity_table.radii,
            );
            errdefer allocator.free(radii);
            const em = try allocator.dupe(
                T,
                ring_results.items[0].emissivity_table.em,
            );
            errdefer allocator.free(em);
            const time = try allocator.dupe(
                T,
                ring_results.items[0].emissivity_table.time,
            );
            errdefer allocator.free(time);
            const energyshift = try allocator.dupe(
                T,
                ring_results.items[0].emissivity_table.g,
            );
            errdefer allocator.free(energyshift);
            const local_theta = try allocator.dupe(
                T,
                ring_results.items[0].emissivity_table.g,
            );
            errdefer allocator.free(local_theta);

            // Zero the emissivity and time arrays before the summation:
            @memset(em, .zero);
            @memset(time, .zero);
            @memset(energyshift, .zero);
            @memset(local_theta, .zero);

            var total_weight: T = .zero;
            var total_photon_fractions: emissivity.PhotonFractions(T.T) = .{};

            // TODO: an intergral that takes into account the GR-corrected
            // volume element and interpolates many different rings instead of
            // using only those computed.
            for (ring_results.items, 0..) |ring, i| {
                const weight = switch (corona) {
                    inline .umbrella, .disc => |c| c.summationWeight(i),
                    else => unreachable,
                };

                for (
                    em,
                    time,
                    energyshift,
                    local_theta,
                    ring.emissivity_table.em,
                    ring.emissivity_table.time,
                    ring.emissivity_table.g,
                    ring.emissivity_table.local_theta,
                ) |*e, *t, *g, *th, r_em, r_time, r_g, r_th| {
                    e.* = A.add(e.*, A.mult(r_em, weight));
                    t.* = A.add(e.*, A.mult(r_time, weight));
                    g.* = A.add(e.*, A.mult(r_g, weight));
                    th.* = A.add(e.*, A.mult(r_th, weight));
                }

                // And do the weighted sum of the photon fractions:
                inline for (@typeInfo(@TypeOf(total_photon_fractions)).@"struct".fields) |field| {
                    @field(total_photon_fractions, field.name) += @field(
                        ring.photon_fractions,
                        field.name,
                    ) * weight.x;
                }

                // Sum the total weight for the normalisation at the end.
                total_weight = A.add(total_weight, weight);
            }

            // Normalise the emissivity and time
            if (total_weight.x > 0) {
                for (em, time, energyshift, local_theta) |*e, *t, *g, *th| {
                    e.* = A.div(e.*, total_weight);
                    t.* = A.div(t.*, total_weight);
                    g.* = A.div(g.*, total_weight);
                    th.* = A.div(th.*, total_weight);
                }

                // And do the weighted sum of the photon fractions:
                inline for (@typeInfo(@TypeOf(total_photon_fractions)).@"struct".fields) |field| {
                    @field(total_photon_fractions, field.name) /= total_weight.x;
                }
            }

            return .{
                .corona = corona,
                .results = try ring_results.toOwnedSlice(allocator),
                .runtime = total_runtime,
                .emissivity_table = .init(
                    radii,
                    em,
                    time,
                    energyshift,
                    local_theta,
                ),
                .photon_fractions = total_photon_fractions,
            };
        }

        fn runLamppost(
            self: *const Self,
            allocator: std.mem.Allocator,
            corona: emissivity.CoronalModel(T),
            threads: *ThreadMap,
            run_opts: RunToolOptions,
        ) !SingleResult {
            std.debug.assert(corona == .lamppost);
            const time_now = std.time.milliTimestamp();

            var cache = try emissivity.EmissivityCache(T).init(allocator, .{
                .max_traces = self.em_opts.num_photons,
            });
            errdefer cache.deinit(allocator);

            const thread_ctx: ThreadContext = .{
                .corona = corona,
                .cache = &cache,
            };

            try runThreads(
                threads,
                cache.traces,
                thread_ctx,
                ThreadContext.workLamppost,
                run_opts,
            );

            cache.sortByRadius();

            const duration = std.time.milliTimestamp() - time_now;

            return .{
                .runtime = duration,
                .emissivity_table = try cache.timeAveraged(allocator),
                .cache = cache,
                .corona = corona,
                .photon_fractions = cache.photonFractions(self.metric),
            };
        }

        fn runRing(
            self: *const Self,
            allocator: std.mem.Allocator,
            corona: emissivity.CoronalModel(T),
            threads: *ThreadMap,
            run_opts: RunToolOptions,
        ) !SingleResult {
            std.debug.assert(corona == .ring);
            const time_now = std.time.milliTimestamp();

            var cache = try emissivity.EmissivityCache(T).init(allocator, .{
                .max_traces = self.em_opts.num_photons,
            });
            errdefer cache.deinit(allocator);

            const thread_ctx: ThreadContext = .{
                .corona = corona,
                .cache = &cache,
            };
            try runThreads(
                threads,
                cache.traces,
                thread_ctx,
                ThreadContext.workRing,
                run_opts,
            );

            // Sort them all by radius, then can perform the 2d rebinning without
            // needing to a 2d grid.
            cache.sortByRadius();

            const time_average = try cache.rebinTimeAveraged(
                allocator,
                .{
                    .n_phi = self.em_opts.num_phi,
                    .r_max = self.em_opts.max_radius,
                    .n_radii = self.em_opts.num_radii,
                },
            );
            errdefer time_average.deinit(allocator);

            const duration = std.time.milliTimestamp() - time_now;
            return .{
                .runtime = duration,
                .emissivity_table = time_average,
                .cache = cache,
                .corona = corona,
                .photon_fractions = cache.photonFractions(self.metric),
            };
        }

        const ThreadContext = struct {
            corona: emissivity.CoronalModel(T),
            cache: *emissivity.EmissivityCache(T),

            fn workLamppost(
                ctx: ThreadContext,
                out: *emissivity.EmissivityTrace(T),
                index: usize,
                _: ThreadMap.ThreadId,
            ) void {
                out.* = ctx.cache.traceLamppost(ctx.corona.lamppost, index);
            }

            fn workRing(
                ctx: ThreadContext,
                out: *emissivity.EmissivityTrace(T),
                index: usize,
                _: ThreadMap.ThreadId,
            ) void {
                out.* = ctx.cache.traceRing(ctx.corona.ring, index);
            }
        };
    };
}
