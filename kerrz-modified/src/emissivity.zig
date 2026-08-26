const std = @import("std");
const ad = @import("zad");
const zfits = @import("zfits");
const dinterp = @import("dinterp");
const options = @import("options");

const utils = @import("utils.zig");

const accretion_discs = @import("accretion-discs.zig");
const geometry = @import("geometry.zig");
const geodesic = @import("geodesic.zig");
const redshift = @import("redshift.zig");
const orbits = @import("orbits.zig");
const interpolations = @import("interpolations.zig");
const iterators = @import("iterators.zig");

const AccretionDisc = accretion_discs.AccretionDisc;
const DualNumber = ad.DualNumber;
const KerrMetric = geometry.KerrMetric;
const FourVector = geometry.FourVector;
const NullGeodesic = geodesic.NullGeodesic;
const TraceResult = geodesic.TraceResult;
const Matrix = @import("matrix.zig").Matrix;

const TEST_TOLERANCE = options.test_numerical_tolerance;

pub fn PhotonFractions(comptime T: type) type {
    return struct {
        const Self = @This();
        disc: T = 0,
        infinity: T = 0,
        event_horizon: T = 0,
        no_status: T = 0,

        above_isco: T = 0,
        below_isco: T = 0,

        pub fn addToHdu(self: *const Self, hdu: *zfits.Hdu) !void {
            try hdu.appendHeaderRecord(
                "PF_BH",
                .{
                    .value = .{ .float = @floatCast(self.event_horizon) },
                    .comment = "Photon fraction event horizon",
                },
            );
            try hdu.appendHeaderRecord(
                "PF_INF",
                .{
                    .value = .{ .float = @floatCast(self.infinity) },
                    .comment = "Photon fraction infinity",
                },
            );
            try hdu.appendHeaderRecord(
                "PF_NONE",
                .{
                    .value = .{ .float = @floatCast(self.no_status) },
                    .comment = "Photon fraction with unknown status",
                },
            );
            try hdu.appendHeaderRecord(
                "PF_DISC",
                .{
                    .value = .{ .float = @floatCast(self.disc) },
                    .comment = "Photon fraction accretion disc",
                },
            );
            try hdu.appendHeaderRecord(
                "PF_AISCO",
                .{
                    .value = .{ .float = @floatCast(self.above_isco) },
                    .comment = "Photon fraction above the ISCO",
                },
            );
            try hdu.appendHeaderRecord(
                "PF_BISCO",
                .{
                    .value = .{ .float = @floatCast(self.below_isco) },
                    .comment = "Photon fraction below the ISCO",
                },
            );
        }
    };
}

/// Values associated with calculating emissivity terms.
///
/// This structure is intended for internal calculations only, and is not
/// supposed to be part of any public interface.
fn EmissivityValues(comptime T: type) type {
    const A = T.Algebra;
    return struct {
        const Self = @This();
        /// The energyshift
        g: T,
        /// The Lorentz factor
        lorentz_factor: T,
        /// The incident intensity
        intensity: T,
        /// The area term on the disc
        area: T,
        /// The emissivity itself
        emissivity: T,
        /// The jacobian term (d(theta, phi) / d(r, phi))
        jacobian: T,
        /// The incident angle, measured from the disc normal.
        local_theta: T,

        /// Calculate the emissivity values. The `jacobian` should be the
        /// d(theta, phi) / d(r, phi) term for integrating over the disc.
        /// `theta` is here the elevation angle in the local sky of the corona.
        ///
        /// Finally, `v_src` is the velocity of the emitter, required for the
        /// energyshift calculations.
        pub fn calculate(
            metric: KerrMetric(T),
            geod: NullGeodesic(T),
            result: TraceResult(T),
            photon_index: T,
            jacobian: T,
            theta: T,
            v_src: FourVector(T),
        ) Self {
            const ts = metric.tangentSpaceAlt(result.r, result.theta);
            const lorentz_factor = orbits.lorentzFactorKeplerian(T, metric, ts);
            // Calculate the energyshift. The inverse must be taken, as we are
            // interested in the energyshift in the frame of the endpoint.
            const g = A.div(.one, redshift.keplerianRedshiftResult(
                T,
                metric,
                geod,
                result,
                v_src,
            ));
            const _I = A.pow(g, photon_index);

            const m = ts.metric_components;

            // Determine the general relativistically corrected area, modulo
            // the special relativistic correction
            const area = A.sqrt(A.mult(m.rr, m.phph));

            // Baker PhD thesis, 2026, Equation (113)
            // Calculate the emissivity itself
            const em = A.div(
                A.mult(A.mult(A.sin(theta), jacobian), _I),
                A.mult(lorentz_factor, area),
            );

            const local_angles = result.localAnglesReverse(
                metric,
                geod,
                orbits.keplerianPlunging(T, metric, ts),
            );

            return .{
                .g = g,
                .lorentz_factor = lorentz_factor,
                .jacobian = jacobian,
                .intensity = _I,
                .area = area,
                .emissivity = em,
                .local_theta = local_angles.theta,
            };
        }
    };
}

pub fn LamppostOptions(comptime T: type) type {
    return struct {
        photon_index: T = .promote(2.0),
        /// The radial out-/inflowing velocity.
        radial_velocity: T = .zero,
        height: T,
    };
}

/// The lamppost coronal model. A single point-like source on the spin axis of
/// the black hole, geometrically parameterised by the height along the `z`
/// axis.
///
/// All methods implemented for the lamppost corona exploit the complete
/// axis-symmetry of the model.
pub fn Lamppost(comptime T: type) type {
    const A = T.Algebra;

    return struct {
        pub const Options = LamppostOptions(T);
        const Self = @This();

        x: FourVector(T),
        v: FourVector(T),
        metric: KerrMetric(T),
        ts: KerrMetric(T).TangentSpace,

        /// The photon index of the powerlaw emission.
        photon_index: T,

        pub fn init(metric: KerrMetric(T), opts: Options) Self {
            const x: FourVector(T) = .{
                .t = .zero,
                .r = opts.height,
                .th = .promote(1e-7), // TODO: allow complete on-axis
                .ph = .zero,
            };
            const ts = metric.tangentSpace(x);

            // Calculate the four-velocity of the source.
            const v = if (opts.radial_velocity.x == 0)
                orbits.stationary(T, ts)
            else
                orbits.radialMotion(T, ts, opts.radial_velocity);

            return .{
                .x = x,
                .v = v,
                .metric = metric,
                .ts = ts,
                .photon_index = opts.photon_index,
            };
        }

        /// Adapt to a new dual number type.
        pub fn adapt(self: Self, comptime NewT: type) Lamppost(NewT) {
            return .{
                .x = self.x.adapt(NewT),
                .v = self.v.adapt(NewT),
                .metric = self.metric.adapt(NewT),
                .ts = self.ts.adapt(NewT),
                .photon_index = self.photon_index.pushSlot(),
            };
        }

        /// Trace a single photon and calculate its contribution to the
        /// emissivity. The photon is parameterised by the elevation angle in
        /// the local frame of the corona, denoted `theta`.
        pub fn traceEmissivity(
            self: Self,
            disc: AccretionDisc(T),
            theta: T,
            opts: geodesic.TracingConfig,
        ) EmissivityTrace(T) {
            const Dual = T.PushSlot();
            const geod = NullGeodesic(Dual).fromSkyAnglesTangentSpace(
                self.metric.adapt(Dual),
                self.ts.adapt(Dual),
                self.v.adapt(Dual),
                // Set the derivative to be calculate with respect to theta
                theta.pushSlot().diff(0),
                .zero,
            );

            const result = geod.traceDisc(self.metric.adapt(Dual), disc.adapt(Dual), opts);

            if (result.status != .intersected_disc) {
                return .zero(result.status);
            }

            const dr_dtheta: T = .promote(@abs(result.r.dx[0]));
            const jacobian = A.div(.one, dr_dtheta);

            std.debug.assert(std.math.isFinite(jacobian.x));

            const values = EmissivityValues(T).calculate(
                self.metric,
                geod.adapt(T),
                result.adapt(T),
                self.photon_index,
                jacobian,
                theta,
                self.v,
            );

            const total = result.adapt(T).totalAntiderivatives(
                self.metric,
                geod.adapt(T),
            );

            return .{
                .r = .adaptFrom(result.r),
                .phi = total.coordinateAzimuth(self.metric, geod.adapt(T)),
                .t = total.coordinateTime(self.metric, geod.adapt(T)),
                .em = .adaptFrom(values.emissivity),
                .g = .adaptFrom(values.g),
                .local_theta = .adaptFrom(values.local_theta),
                .status = result.status,
            };
        }
    };
}

test "lamppost emissivity" {
    const Dual = DualNumber(f64, 1);
    const Lamp = Lamppost(Dual);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    const disc: accretion_discs.ThinDisc(Dual) = .{
        .inner_radius = .zero,
        .outer_radius = .promote(10000000),
    };
    const lp = Lamp.init(metric, .{ .height = .promote(4.0) });

    const em = lp.traceEmissivity(
        .{ .thin_disc = disc },
        .promote(std.math.degreesToRadians(82.0)),
        .{},
    );

    // TODO: this is a regression test, it has not yet been verified.
    try std.testing.expectApproxEqAbs(
        0.0026342757362146406,
        em.em.x,
        TEST_TOLERANCE,
    );
}

pub fn RingOptions(comptime T: type) type {
    return struct {
        photon_index: T = .promote(2.0),
        velocity: orbits.VelocityProfiles = .co_rotate,
        height: T,
        radius: T,
    };
}

/// The ring-like coronal model. This model is parameterised by a radius `x`
/// and a height along the spin axis `z`, such that the ring is always parallel
/// to the equatorial plane.
///
/// For the purposes of calculations, the axis-symmetry of the ring is
/// exploited, and a single off-axis point can be averaged over the azimuthal
/// coordinate of the disc to obtain the axisymmetric emissivity profile.
pub fn Ring(comptime T: type) type {
    const A = T.Algebra;

    return struct {
        pub const Options = RingOptions(T);
        const Self = @This();

        // The coordinate of the off-axis point of the ring
        x: FourVector(T),
        v: FourVector(T),
        metric: KerrMetric(T),
        ts: KerrMetric(T).TangentSpace,

        /// The photon index of the powerlaw emission.
        photon_index: T,

        // TODO: provide alt constructor for (R, theta) coords

        pub fn init(metric: KerrMetric(T), opts: Options) Self {
            // Convert the coordinates
            const r = A.sqrt(A.add(A.powi(opts.radius, 2), A.powi(opts.height, 2)));
            const th = A.atan2(opts.radius, opts.height);
            const x: FourVector(T) = .{
                .t = .zero,
                .r = r,
                .th = th,
                .ph = .zero,
            };
            const ts = metric.tangentSpace(x);
            // TODO: don't use the velocity here, but calculate the frame and
            // store that, since it will not change.
            const v = opts.velocity.fourVector(T, metric, ts);
            return .{
                .x = x,
                .v = v,
                .metric = metric,
                .ts = ts,
                .photon_index = opts.photon_index,
            };
        }

        /// Trace a single photon and calculate its contribution to the
        /// emissivity. The photon is parameterised by the elevation and
        /// azimuthal angle in the local frame of the corona, denoted `theta`
        /// and `phi` respectively.
        pub fn traceEmissivity(
            self: Self,
            disc: AccretionDisc(T),
            theta: T,
            phi: T,
            opts: geodesic.TracingConfig,
        ) EmissivityTrace(T) {
            const Dual = T.PushSlot().PushSlot();

            const geod = NullGeodesic(Dual).fromSkyAnglesRotatedTangentSpace(
                self.metric.adapt(Dual),
                self.ts.adapt(Dual),
                self.v.adapt(Dual),
                // Set the derivative to be calculate with respect to theta and
                // phi
                theta.pushSlot().pushSlot().diff(0),
                phi.pushSlot().pushSlot().diff(1),
                .{ .theta_0 = .promote(std.math.pi / 2.0) },
            );

            const result = geod.traceDisc(
                self.metric.adapt(Dual),
                disc.adapt(Dual),
                opts,
            );

            if (result.status != .intersected_disc) {
                return .zero(result.status);
            }

            const total = result.totalAntiderivatives(self.metric.adapt(Dual), geod);
            const azimuth = total.coordinateAzimuth(self.metric.adapt(Dual), geod);

            const det_jac: T = .promote(@abs(ad.jacobianDeterminant(
                Dual,
                result.r,
                azimuth,
            )));

            const jacobian = A.div(.one, det_jac);

            const values = EmissivityValues(T).calculate(
                self.metric,
                geod.adapt(T),
                result.adapt(T),
                self.photon_index,
                jacobian,
                theta,
                self.v,
            );

            return .{
                .r = .adaptFrom(result.r),
                .phi = .adaptFrom(azimuth),
                .t = .adaptFrom(total.coordinateTime(self.metric.adapt(Dual), geod)),
                .em = .adaptFrom(values.emissivity),
                .g = .adaptFrom(values.g),
                .local_theta = .adaptFrom(values.local_theta),
                .status = result.status,
            };
        }

        /// Utility function to obtain the height of the ring.
        pub fn getHeight(self: Self) T {
            return A.mult(self.x.r, A.cos(self.x.th));
        }

        /// Utility function to obtain the radius of the ring.
        pub fn getRadius(self: Self) T {
            return A.mult(self.x.r, A.sin(self.x.th));
        }
    };
}

test "ring emissivity" {
    const Dual = DualNumber(f64, 0);
    const RingModel = Ring(Dual);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    const disc: accretion_discs.ThinDisc(Dual) = .{
        .inner_radius = .zero,
        .outer_radius = .promote(10000000),
    };
    const ring = RingModel.init(metric, .{
        .height = .promote(4.0),
        .radius = .promote(4.0),
    });

    const em = ring.traceEmissivity(
        .{ .thin_disc = disc },
        .promote(std.math.degreesToRadians(82.0)),
        .promote(std.math.degreesToRadians(2.0)),
        .{},
    );

    // TODO: this is a regression test, it has not yet been verified.
    try std.testing.expectApproxEqAbs(
        0.06159082905384861,
        em.em.x,
        TEST_TOLERANCE,
    );
}

/// Options for a disc-like corona, summed together from many ring-like
/// elements.
pub fn DiscOptions(comptime T: type) type {
    return struct {
        photon_index: T = .promote(2.0),
        velocity: orbits.VelocityProfiles = .co_rotate,
        height: T,
        inner_radius: T,
        outer_radius: T,
        n_rings: usize = 20,
    };
}

pub fn Disc(comptime T: type) type {
    return struct {
        pub const Options = DiscOptions(T);
        const Self = @This();

        /// The rings this disc is made out of
        rings: []Ring(T),
        /// The radii associated with the rings. This is `rings.len + 1`, so
        /// the lower and upper radii associated with each ring.
        radii: []const T.T,

        pub fn init(
            allocator: std.mem.Allocator,
            metric: KerrMetric(T),
            opts: Options,
        ) !Self {
            const rings = try allocator.alloc(Ring(T), opts.n_rings);
            errdefer allocator.free(rings);
            // TODO: use a dual number iterator
            const radii = try iterators.Grid(T.T).linear.fill(
                allocator,
                opts.inner_radius.x,
                opts.outer_radius.x,
                rings.len + 1,
            );
            errdefer allocator.free(radii);

            for (rings, radii[0..rings.len]) |*ring, r| {
                ring.* = .init(metric, .{
                    .photon_index = opts.photon_index,
                    .velocity = opts.velocity,
                    .height = opts.height,
                    .radius = .promote(r),
                });
            }

            return .{
                .rings = rings,
                .radii = radii,
            };
        }

        pub fn deinit(self: *const Self, allocator: std.mem.Allocator) void {
            allocator.free(self.rings);
            allocator.free(self.radii);
        }

        /// The summation weight of the `ith` ring when combining the rings
        /// together into the full disc. This includes all relativistic
        /// correction.
        pub fn summationWeight(self: *const Self, i: usize) T {
            const A = T.Algebra;

            const ring = self.rings[i];

            const r_low = self.radii[i];
            const r_high = self.radii[i + 1];

            const lorentz_factor = orbits.lorentzFactor(
                T,
                ring.ts,
                ring.v,
            );

            const g = ring.ts.metric_components;

            // The weight of the annulus, 2 pi r dr, along with the
            // relativistically corrected area element along the surface of
            // `z=h`.
            const weight = A.mult(
                A.mult(lorentz_factor, A.sqrt(g.rr)),
                .promote((r_high - r_low) * r_low * 2 * std.math.pi),
            );

            return weight;
        }
    };
}

/// Options for an umbrella-like corona, summed together from many ring-like
/// elements.
pub fn UmbrellaOptions(comptime T: type) type {
    return struct {
        photon_index: T = .promote(2.0),
        velocity: orbits.VelocityProfiles = .co_rotate,
        offset_radius: T,
        inner_opening_angle: T,
        outer_opening_angle: T,
        n_rings: usize = 20,
    };
}

pub fn Umbrella(comptime T: type) type {
    return struct {
        pub const Options = UmbrellaOptions(T);
        const Self = @This();

        /// The rings this umbrella is made of.
        rings: []Ring(T),
        /// The polar angles associated with the rings. This is `rings.len +
        /// 1`, so a lower and upper angle are associated with each ring.
        angles: []const T.T,

        pub fn init(
            allocator: std.mem.Allocator,
            metric: KerrMetric(T),
            opts: Options,
        ) !Self {
            const A = T.Algebra;

            const rings = try allocator.alloc(Ring(T), opts.n_rings);
            errdefer allocator.free(rings);

            // TODO: use a dual number iterator
            const angles = try iterators.Grid(T.T).linear.fill(
                allocator,
                opts.inner_opening_angle.x,
                opts.outer_opening_angle.x,
                rings.len + 1,
            );
            errdefer allocator.free(angles);

            for (rings, angles[0..rings.len]) |*ring, theta| {
                const h = A.mult(
                    opts.offset_radius,
                    A.cos(.promote(theta)),
                );
                const x = A.mult(
                    opts.offset_radius,
                    A.sin(.promote(theta)),
                );

                ring.* = .init(metric, .{
                    .photon_index = opts.photon_index,
                    .velocity = opts.velocity,
                    .height = h,
                    .radius = x,
                });
            }

            return .{
                .rings = rings,
                .angles = angles,
            };
        }

        pub fn deinit(self: *const Self, allocator: std.mem.Allocator) void {
            allocator.free(self.rings);
            allocator.free(self.angles);
        }

        /// The summation weight of the `ith` ring when combining the rings
        /// together into the full disc. This includes all relativistic
        /// correction.
        pub fn summationWeight(self: *const Self, i: usize) T {
            const A = T.Algebra;

            const ring = self.rings[i];

            const th_low = self.angles[i];
            const th_high = self.angles[i + 1];

            const lorentz_factor = orbits.lorentzFactor(
                T,
                ring.ts,
                ring.v,
            );

            const g = ring.ts.metric_components;

            // The weight here is from considering a patch of the sphere:
            //
            //     dA = r^2 sin(θ) dθ dφ
            //
            // Since `dr = 0` there is no `g_rr` term, and the azimuthal
            // coordinate can be integrated out:
            //
            //     dA = 2π r^2 sin(θ) dθ
            //
            const dth_sinth = (th_high - th_low) * @sin(th_low);
            const r2 = A.powi(ring.ts.x.r, 2);
            const weight = A.mult(
                A.mult(lorentz_factor, A.sqrt(g.thth)),
                A.mult(r2, .promote(dth_sinth * 2 * std.math.pi)),
            );

            return weight;
        }
    };
}

/// An emissivity trace.
///
/// This is what the `EmissivityValues` are unpacked to as part of the public
/// interface.
///
/// The `status` field is augmented from the `EmissivityValues`, and therefore
/// the `EmissivityTrace` can hold information also about photons that did not
/// hit the accretion disc. The `status` field should always be checked before
/// accessing any of the values.
pub fn EmissivityTrace(comptime T: type) type {
    return struct {
        const Self = @This();
        /// The radius on the accretion disc.
        r: T,
        /// The azimuthal coordinate on the accretion disc.
        phi: T,
        /// The corona-to-disc light travel time.
        t: T,
        /// The corona-to-disc energyshift.
        g: T,
        /// The emissivity itself.
        em: T,
        /// The incident angle, measured from the disc normal.
        local_theta: T,

        /// The status of the geodesic. If this is not `intersected_disc`, then
        /// all the above values will be zeroed.
        status: geodesic.Status,

        fn zero(status: geodesic.Status) Self {
            return .{
                .r = .zero,
                .phi = .zero,
                .t = .zero,
                .g = .zero,
                .em = .zero,
                .local_theta = .zero,
                .status = status,
            };
        }

        fn sortRadius(_: void, left: Self, right: Self) bool {
            return left.r.x < right.r.x;
        }
    };
}

/// The cache used to calculate emissivity profiles.
pub fn EmissivityCache(comptime T: type) type {
    return struct {
        /// Options for controlling the emissivity profile calculation.
        pub const Options = struct {
            /// The maximum number of traces to include in the calculation.
            max_traces: usize = 3000,
            /// The accretion disc to calculate the emissivity profile over.
            disc: AccretionDisc(T) = .{
                .thin_disc = .{ .inner_radius = .zero, .outer_radius = .promote(std.math.inf(T.T)) },
            },
            /// The angular offset around the poles to avoid singularities.
            eps: T.T = 1e-6,

            /// Radial grid size for rebinning the emissivity profile, if needed.
            radial_grid_size: usize = 1000,
        };

        const Self = @This();

        /// The trace buffer.
        traces: []EmissivityTrace(T),
        /// Options passed during initialisation.
        opts: Options,

        /// Caller must manage memory. Initialises without pre-allocating any
        /// cache.
        pub fn initEmpty(opts: Options) Self {
            return .{
                .traces = &.{},
                .opts = opts,
            };
        }

        /// Caller owns memory and must call deinit.
        pub fn init(allocator: std.mem.Allocator, opts: Options) !Self {
            return .{
                .traces = try allocator.alloc(EmissivityTrace(T), opts.max_traces),
                .opts = opts,
            };
        }

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            allocator.free(self.traces);
            self.* = undefined;
        }

        /// Compute the photon fractions for the emissivity point traces in
        /// this cache.
        pub fn photonFractions(self: *const Self, metric: KerrMetric(T)) PhotonFractions(T.T) {
            var fractions: PhotonFractions(usize) = .{};
            for (self.traces) |trace| {
                switch (trace.status) {
                    .event_horizon => fractions.event_horizon += 1,
                    .infinity => fractions.infinity += 1,
                    .intersected_disc => {
                        fractions.disc += 1;
                        if (trace.r.x > metric.isco.x) {
                            fractions.above_isco += 1;
                        } else {
                            fractions.below_isco += 1;
                        }
                    },
                    .no_status => fractions.no_status += 1,
                }
            }

            const total: T.T = @floatFromInt(
                fractions.disc + fractions.event_horizon +
                    fractions.infinity + fractions.no_status,
            );

            return .{
                .disc = @as(T.T, @floatFromInt(fractions.disc)) / total,
                .above_isco = @as(T.T, @floatFromInt(fractions.above_isco)) / total,
                .below_isco = @as(T.T, @floatFromInt(fractions.below_isco)) / total,
                .event_horizon = @as(T.T, @floatFromInt(fractions.event_horizon)) / total,
                .infinity = @as(T.T, @floatFromInt(fractions.infinity)) / total,
                .no_status = @as(T.T, @floatFromInt(fractions.no_status)) / total,
            };
        }

        fn convertIndex(self: Self, index: usize) T.T {
            const index_f: T.T = @floatFromInt(index);
            const max_f: T.T = @floatFromInt(self.opts.max_traces);
            // Avoid the poles.
            return ((index_f / max_f) + self.opts.eps) / (1 - self.opts.eps * 2);
        }

        /// Used to trace a single photon for a lamppost coronal model.
        ///
        /// The index parameter should be a monotonically increasing value
        /// between 0 and `self.opts.max_traces`.
        pub fn traceLamppost(self: Self, lp: Lamppost(T), index: usize) EmissivityTrace(T) {
            const x = self.convertIndex(index);
            const theta = x * std.math.pi;
            return lp.traceEmissivity(
                self.opts.disc,
                .promote(theta),
                .{},
            );
        }

        /// Used to trace a single photon for a ring-like coronal model.
        ///
        /// The index parameter should be a monotonically increasing value
        /// between 0 and `self.opts.max_traces`.
        pub fn traceRing(self: Self, ring: Ring(T), index: usize) EmissivityTrace(T) {
            const x = self.convertIndex(index);
            const index_f = x * @as(T.T, @floatFromInt(self.opts.max_traces));

            // This is the Golden spiral method for generating evenly distributed point
            // on the surface of a sphere as the limit of `n` goes to infinity.
            //
            // TODO: this algorithm should be a choice, and appropriate
            // weightings used in computing e.g. photon fractions or emissivity
            // profiles.
            const phi = @mod(std.math.pi * (1 + @sqrt(5.0)) * index_f, std.math.pi * 2.0);
            const theta = x * std.math.pi;

            return ring.traceEmissivity(
                self.opts.disc,
                .promote(theta),
                .promote(phi),
                .{},
            );
        }

        /// Sort all of the traces by radius.
        pub fn sortByRadius(self: *Self) void {
            std.sort.heap(
                EmissivityTrace(T),
                self.traces,
                {},
                EmissivityTrace(T).sortRadius,
            );
        }

        /// Options for controlling how the emissivity profile is rebinned.
        pub const BinningOptions = struct {
            /// The innermost radius on the disc for the radial grid.
            r_min: T.T = 1.0,
            /// The outermost radius on the disc for the radial grid.
            r_max: T.T = 1e4,
            /// The number of radial bins (logarithmically spaced).
            n_radii: usize = 500,
            /// The number of phi bins, used internally for interpolating
            /// values.
            n_phi: usize = 1000,

            /// The minimal value of the time axis in time-dependent rebinning.
            t_min: T.T = 0.0,
            /// The maximal value of the time axis in time-dependent rebinning.
            t_max: T.T = 1e4,
            /// The number of time binds.
            n_time: usize = 500,
        };

        fn predicatePhi(lhs: T.T, rhs: T.T) bool {
            return lhs > rhs;
        }

        fn findPhiIndex(phis: []const T.T, phi: T.T) usize {
            return std.sort.partitionPoint(
                T.T,
                phis[0 .. phis.len - 1],
                @mod(phi, std.math.pi * 2.0),
                predicatePhi,
            );
        }

        /// Used for time-dependent emissivities.
        ///
        /// The temporary and output arrays needed for rebinning the emissivity
        /// profile in various ways.
        pub const BinningCache = struct {
            // These booleans are used to track what needs to be freed.
            owns_radii: bool = true,
            owns_emissivity: bool = true,
            owns_time: bool = true,
            owns_energyshift: bool = true,
            owns_local_theta: bool = true,
            owns_phi: bool = true,

            radii: []T.T,
            phis: []T.T,
            em: Matrix(T.T),
            time: Matrix(T.T),
            g: Matrix(T.T),
            local_theta: Matrix(T.T),
            counts: Matrix(usize),

            pub fn deinit(
                self: BinningCache,
                allocator: std.mem.Allocator,
            ) void {
                if (self.owns_emissivity) {
                    self.em.deinit(allocator);
                }
                if (self.owns_time) {
                    self.time.deinit(allocator);
                }
                if (self.owns_energyshift) {
                    self.g.deinit(allocator);
                }
                if (self.owns_radii) {
                    allocator.free(self.radii);
                }
                if (self.owns_local_theta) {
                    self.local_theta.deinit(allocator);
                }
                if (self.owns_phi) {
                    allocator.free(self.phis);
                }
                self.counts.deinit(allocator);
            }

            fn toOwnedRadii(self: *BinningCache) []T.T {
                self.owns_radii = false;
                return self.radii;
            }

            fn toOwnedPhi(self: *BinningCache) []T.T {
                self.owns_phi = false;
                return self.phis;
            }

            fn toOwnedEmissivity(self: *BinningCache) Matrix(T.T) {
                self.owns_emissivity = false;
                return self.em;
            }

            fn toOwnedTime(self: *BinningCache) Matrix(T.T) {
                self.owns_time = false;
                return self.time;
            }

            fn toOwnedEnergyshift(self: *BinningCache) Matrix(T.T) {
                self.owns_energyshift = false;
                return self.g;
            }

            fn toOwnedLocalTheta(self: *BinningCache) Matrix(T.T) {
                self.owns_local_theta = false;
                return self.local_theta;
            }

            fn reset(self: *BinningCache) void {
                @memset(self.em.values, 0);
                @memset(self.time.values, 0);
                @memset(self.g.values, 0);
                @memset(self.local_theta.values, 0);
                @memset(self.counts.values, 0);
            }

            pub fn init(
                allocator: std.mem.Allocator,
                opts: BinningOptions,
            ) !BinningCache {
                const radii = try iterators.Grid(T.T).log10.fill(
                    allocator,
                    opts.r_min,
                    opts.r_max,
                    opts.n_radii,
                );
                errdefer allocator.free(radii);

                // TODO: should this be to (1 - 1/n) * 2 pi?
                const phis = try iterators.Grid(T.T).linear.fill(
                    allocator,
                    0,
                    2 * std.math.pi,
                    opts.n_phi,
                );
                errdefer allocator.free(phis);

                const em = try Matrix(T.T).init(allocator, radii.len, phis.len);
                errdefer em.deinit(allocator);
                @memset(em.values, 0);

                const time = try Matrix(T.T).init(allocator, radii.len, phis.len);
                errdefer time.deinit(allocator);
                @memset(time.values, 0);

                const g = try Matrix(T.T).init(allocator, radii.len, phis.len);
                errdefer g.deinit(allocator);
                @memset(g.values, 0);

                const local_theta = try Matrix(T.T).init(allocator, radii.len, phis.len);
                errdefer local_theta.deinit(allocator);
                @memset(local_theta.values, 0);

                const counts = try Matrix(usize).init(allocator, radii.len, phis.len);
                errdefer counts.deinit(allocator);
                @memset(counts.values, 0);

                return .{
                    .radii = radii,
                    .em = em,
                    .time = time,
                    .g = g,
                    .local_theta = local_theta,
                    .phis = phis,
                    .counts = counts,
                };
            }

            /// Equivalently obtain the time- and radially-averaged emissivity
            /// profile. Bin the emissivity data by radial coordinate on the
            /// accretion disc.
            fn binTraces(
                cache: *BinningCache,
                traces: []const EmissivityTrace(T),
            ) void {
                // Rebin by running over all points and the grid, advancing one bin
                // at a time.
                // The idea is to take everything between some r and r + dr and bin
                // it azimuthally, and then to interpolate the values around the
                // annulus and put them into the output array.
                var r_index: usize = 0;
                outer: for (traces) |point| {
                    if (point.em.x == 0) continue;

                    while (point.r.x > cache.radii[r_index + 1]) {
                        r_index += 1;
                        // Are we done?
                        if (r_index == cache.radii.len - 1) break :outer;
                    }

                    const index = findPhiIndex(cache.phis, point.phi.x);
                    cache.counts.getPtr(index, r_index).* += 1;
                    cache.em.getPtr(index, r_index).* += point.em.x;
                    cache.time.getPtr(index, r_index).* += point.t.x;
                    cache.g.getPtr(index, r_index).* += point.g.x;
                    cache.local_theta.getPtr(index, r_index).* += point.local_theta.x;
                }

                // Normalise by the number of counts in each bin.
                for (
                    cache.em.values,
                    cache.time.values,
                    cache.g.values,
                    cache.local_theta.values,
                    cache.counts.values,
                ) |*em, *t, *g, *th, count| {
                    if (count > 0) {
                        const count_f: T.T = @floatFromInt(count);
                        t.* /= count_f;
                        em.* /= count_f;
                        g.* /= count_f;
                        th.* /= count_f;
                    }
                }
            }
        };

        /// The temporary and output arrays needed for rebinning the emissivity
        /// profile in various ways. Used for time-averages.
        pub const AveragedBinningCache = struct {
            // These booleans are used to track what needs to be freed.
            owns_radii: bool = true,
            owns_emissivity: bool = true,
            owns_energyshift: bool = true,
            owns_local_theta: bool = true,
            owns_time: bool = true,

            radii: []T,
            em: []T,
            time: []T,
            g: []T,
            local_theta: []T,
            phis: []T.T,
            tmp_em: []T.T,
            tmp_time: []T.T,
            tmp_g: []T.T,
            tmp_local_theta: []T.T,
            counts: []usize,

            pub fn deinit(
                self: AveragedBinningCache,
                allocator: std.mem.Allocator,
            ) void {
                if (self.owns_emissivity) {
                    allocator.free(self.em);
                }
                if (self.owns_time) {
                    allocator.free(self.time);
                }
                if (self.owns_energyshift) {
                    allocator.free(self.g);
                }
                if (self.owns_local_theta) {
                    allocator.free(self.local_theta);
                }
                if (self.owns_radii) {
                    allocator.free(self.radii);
                }
                allocator.free(self.phis);
                allocator.free(self.tmp_em);
                allocator.free(self.tmp_time);
                allocator.free(self.tmp_g);
                allocator.free(self.tmp_local_theta);
                allocator.free(self.counts);
            }

            fn toOwnedRadii(self: *AveragedBinningCache) []T {
                self.owns_radii = false;
                return self.radii;
            }

            fn toOwnedEmissivity(self: *AveragedBinningCache) []T {
                self.owns_emissivity = false;
                return self.em;
            }

            fn toOwnedTime(self: *AveragedBinningCache) []T {
                self.owns_time = false;
                return self.time;
            }

            fn toOwnedLocalTheta(self: *AveragedBinningCache) []T {
                self.owns_local_theta = false;
                return self.local_theta;
            }

            fn toOwnedEnergyshift(self: *AveragedBinningCache) []T {
                self.owns_energyshift = false;
                return self.g;
            }

            fn reset(self: *AveragedBinningCache) void {
                @memset(self.em, .zero);
                @memset(self.time, .zero);
                @memset(self.g, .zero);
                @memset(self.local_theta, .zero);
                @memset(self.tmp_time, 0);
                @memset(self.tmp_em, 0);
                @memset(self.tmp_g, 0);
                @memset(self.tmp_local_theta, 0);
                @memset(self.counts, 0);
            }

            pub fn init(
                allocator: std.mem.Allocator,
                opts: BinningOptions,
            ) !AveragedBinningCache {
                var r_itt = iterators.Grid(T.T).log10.iterator(
                    opts.r_min,
                    opts.r_max,
                    opts.n_radii,
                );
                const radii = try allocator.alloc(T, opts.n_radii);
                errdefer allocator.free(radii);

                for (radii) |*r| r.* = .promote(r_itt.next().?);

                const em = try allocator.alloc(T, opts.n_radii);
                errdefer allocator.free(em);
                @memset(em, .zero);

                const time = try allocator.alloc(T, opts.n_radii);
                errdefer allocator.free(time);
                @memset(time, .zero);

                const g = try allocator.alloc(T, opts.n_radii);
                errdefer allocator.free(g);
                @memset(g, .zero);

                const local_theta = try allocator.alloc(T, opts.n_radii);
                errdefer allocator.free(local_theta);
                @memset(local_theta, .zero);

                // TODO: should this be to (1 - 1/n) * 2 pi?
                const phis = try iterators.Grid(T.T).linear.fill(
                    allocator,
                    0,
                    2 * std.math.pi,
                    opts.n_phi,
                );
                errdefer allocator.free(phis);

                const tmp_em = try allocator.alloc(T.T, phis.len);
                errdefer allocator.free(tmp_em);
                @memset(tmp_em, 0);

                const tmp_time = try allocator.alloc(T.T, phis.len);
                errdefer allocator.free(tmp_time);
                @memset(tmp_time, 0);

                const tmp_g = try allocator.alloc(T.T, phis.len);
                errdefer allocator.free(tmp_g);
                @memset(tmp_g, 0);

                const tmp_local_theta = try allocator.alloc(T.T, phis.len);
                errdefer allocator.free(tmp_local_theta);
                @memset(tmp_local_theta, 0);

                const counts = try allocator.alloc(usize, phis.len);
                errdefer allocator.free(counts);
                @memset(counts, 0);

                return .{
                    .radii = radii,
                    .em = em,
                    .time = time,
                    .phis = phis,
                    .g = g,
                    .local_theta = local_theta,
                    .tmp_em = tmp_em,
                    .tmp_time = tmp_time,
                    .tmp_g = tmp_g,
                    .tmp_local_theta = tmp_local_theta,
                    .counts = counts,
                };
            }

            /// Equivalently obtain the time- and radially-averaged emissivity
            /// profile. Bin the emissivity data by radial coordinate on the
            /// accretion disc.
            fn binTraces(
                cache: *AveragedBinningCache,
                traces: []const EmissivityTrace(T),
            ) void {
                const A = T.Algebra;
                // Rebin by running over all points and the grid, advancing one bin
                // at a time.
                // The idea is to take everything between some r and r + dr and bin
                // it azimuthally, and then to interpolate the values around the
                // annulus and put them into the output array.
                var r_index: usize = 0;
                for (traces) |point| {
                    if (point.em.x == 0) continue;

                    if (point.r.x > cache.radii[r_index + 1].x) {
                        // Apply what we have so far:
                        for (cache.counts, cache.tmp_em) |count, *em| {
                            const weight: T.T = @floatFromInt(@max(1, count));
                            em.* = em.* / weight;
                        }
                        interpolations.interpolateZeroes(T.T, cache.tmp_em);
                        interpolations.interpolateZeroes(T.T, cache.tmp_time);

                        const len_phi_f: T.T = @floatFromInt(cache.phis.len);
                        for (cache.tmp_em) |em| {
                            cache.em[r_index + 1] = A.add(
                                cache.em[r_index + 1],
                                .promote(em / len_phi_f),
                            );
                        }
                        for (cache.tmp_time) |t| {
                            cache.time[r_index + 1] = A.add(
                                cache.time[r_index + 1],
                                .promote(t / len_phi_f),
                            );
                        }
                        for (cache.tmp_g) |g| {
                            cache.g[r_index + 1] = A.add(
                                cache.g[r_index + 1],
                                .promote(g / len_phi_f),
                            );
                        }
                        for (cache.tmp_local_theta) |local_theta| {
                            cache.local_theta[r_index + 1] = A.add(
                                cache.local_theta[r_index + 1],
                                .promote(local_theta / len_phi_f),
                            );
                        }

                        // Reset the buffers.
                        @memset(cache.counts, 0);
                        @memset(cache.tmp_em, 0);
                        @memset(cache.tmp_time, 0);
                        @memset(cache.tmp_g, 0);
                        @memset(cache.tmp_local_theta, 0);

                        while (point.r.x > cache.radii[r_index + 1].x) {
                            r_index += 1;
                            // Are we done?
                            if (r_index == cache.radii.len - 1) return;
                        }
                    }

                    const index = findPhiIndex(cache.phis, point.phi.x);
                    cache.counts[index] += 1;
                    cache.tmp_em[index] += point.em.x;
                    cache.tmp_time[index] += point.t.x;
                    cache.tmp_g[index] += point.g.x;
                    cache.tmp_local_theta[index] += point.local_theta.x;
                }
            }
        };

        /// Interprets the traces directly as complete axisymmetric emissivity
        /// profile without any rebinning on the disc. This is only valid for
        /// on-axis and axisymmetric sources, such as the lamppost model.
        pub fn timeAveraged(self: Self, allocator: std.mem.Allocator) !AxisymmetricEmissivity(T) {

            // Count how many are non-zero radius
            var start: usize = 0;
            for (self.traces, 0..) |t, i| {
                if (t.r.x != 0) {
                    start = i;
                    break;
                }
            }

            const N = self.traces.len - start;

            const radii = try allocator.alloc(T, N);
            errdefer allocator.free(radii);

            const em = try allocator.alloc(T, N);
            errdefer allocator.free(em);

            const time = try allocator.alloc(T, N);
            errdefer allocator.free(time);

            const energyshift = try allocator.alloc(T, N);
            errdefer allocator.free(energyshift);

            const local_theta = try allocator.alloc(T, N);
            errdefer allocator.free(local_theta);

            for (
                radii,
                em,
                time,
                energyshift,
                local_theta,
                self.traces[start..],
            ) |*r, *e, *t, *g, *th, trace| {
                r.* = trace.r;
                e.* = trace.em;
                t.* = trace.t;
                g.* = trace.g;
                th.* = trace.local_theta;
            }

            return .{
                .em = em,
                .radii = radii,
                .time = time,
                .g = energyshift,
                .local_theta = local_theta,
            };
        }

        /// Used to indicate whether to keep ownership of the memory or
        /// transfer it to the caller.
        const TransferOrKeep = enum { keep, transfer };

        /// Rebin the emissivity profile into the time-averaged emissivity as a
        /// function of radius.
        ///
        /// Caller owns the memory and must free the resulting table.
        pub fn rebinTimeAveraged(
            self: Self,
            allocator: std.mem.Allocator,
            opts: BinningOptions,
        ) !AxisymmetricEmissivity(T) {
            var cache = try AveragedBinningCache.init(allocator, opts);
            defer cache.deinit(allocator);
            return self.rebinTimeAveragedInplace(&cache, .transfer);
        }

        /// Same as `rebinTimeAveraged` but with a pre-allocated cache.
        pub fn rebinTimeAveragedInplace(
            self: Self,
            cache: *AveragedBinningCache,
            keep: TransferOrKeep,
        ) AxisymmetricEmissivity(T) {
            cache.reset();
            // TODO: assert that the radii are sorted.
            cache.binTraces(self.traces);
            return switch (keep) {
                .transfer => .{
                    .em = cache.toOwnedEmissivity(),
                    .radii = cache.toOwnedRadii(),
                    .time = cache.toOwnedTime(),
                    .g = cache.toOwnedEnergyshift(),
                    .local_theta = cache.toOwnedLocalTheta(),
                },
                .keep => .{
                    .em = cache.em,
                    .radii = cache.radii,
                    .time = cache.time,
                    .g = cache.g,
                    .local_theta = cache.local_theta,
                },
            };
        }

        /// Rebin the emissivity profile into a time-dependent grid.
        ///
        /// Caller owns the memory and must free the resulting table.
        pub fn rebinTimeDependent(
            self: Self,
            allocator: std.mem.Allocator,
            opts: BinningOptions,
        ) !TableEmissivity(T.T) {
            var cache = try BinningCache.init(allocator, opts);
            defer cache.deinit(allocator);
            return self.rebinTimeDependentInplace(&cache, .transfer);
        }

        /// Same as `rebinTimeDependent` but with a pre-allocated cache.
        pub fn rebinTimeDependentInplace(
            self: Self,
            cache: *BinningCache,
            keep: TransferOrKeep,
        ) TableEmissivity(T.T) {
            cache.reset();
            // TODO: assert that the radii are sorted.
            cache.binTraces(self.traces);
            return switch (keep) {
                .transfer => .{
                    .r_grid = cache.toOwnedRadii(),
                    .phi_grid = cache.toOwnedPhi(),
                    .em = cache.toOwnedEmissivity(),
                    .time = cache.toOwnedTime(),
                    .g = cache.toOwnedEnergyshift(),
                    .local_theta = cache.toOwnedLocalTheta(),
                },
                .keep => .{
                    .r_grid = cache.radii,
                    .phi_grid = cache.phis,
                    .em = cache.em,
                    .time = cache.time,
                    .g = cache.g,
                    .local_theta = cache.local_theta,
                },
            };
        }
    };
}

/// The emissivity row is what is serialised to the various FITS tables.
pub fn EmissivityRow(comptime T: type) type {
    return struct {
        /// The emissivity itself.
        em: T,
        /// The corona-to-disc light travel time.
        t: T,
        /// The corona-to-disc energyshift.
        g: T,
        /// The incident angle, measured from the disc normal.
        local_theta: T,
    };
}

/// A constant powerlaw index emissivity function.
pub fn PowerLaw(comptime T: type) type {
    const _default_alpha = if (@typeInfo(T) == .float) -3 else T.promote(-3);
    return struct {
        const Self = @This();
        // TODO: some kind of normalisation?
        alpha: T = _default_alpha,

        /// Calculate the powerlaw emissivity at a particular radius.
        pub fn radialValues(self: Self, r: T) EmissivityRow(T) {
            const A = T.Algebra;
            return .{
                .em = A.pow(r, self.alpha),
                .t = .zero,
                .g = .zero,
                .local_theta = .zero,
            };
        }
    };
}

/// An emissivity table, that contains the emissivity as a function of disc
/// parameters. The table does not own any memory it refers to, and the table
/// data must outlive the lifetime of the AxisymmetricEmissivity.
pub fn AxisymmetricEmissivity(comptime T: type) type {
    if (@typeInfo(T) == .float) {
        @compileError("In AxisymmetricEmissivity, T must be a dual number type");
    }
    return struct {
        const Self = @This();
        /// The radii on the disc the each of the emissivity values corresponds to.
        radii: []const T,
        /// The emissivity slice itself.
        em: []const T,
        /// The corona-to-disc time.
        time: []const T,
        /// The corona-to-disc energyshift.
        g: []const T,
        /// The incident angle, measured from the disc normal.
        local_theta: []const T,

        /// Calculate the table emissivity at a particular radius. Returns zero
        /// if out of bounds for the table.
        pub fn radialValues(self: Self, r: T) EmissivityRow(T) {
            if (r.x < self.radii[0].x or r.x > self.radii[self.radii.len - 1].x) {
                return .{
                    .em = .zero,
                    .t = .zero,
                    .g = .zero,
                    .local_theta = .zero,
                };
            }

            const index_of_first_greater = std.sort.partitionPoint(
                T,
                self.radii,
                r,
                PartitionSorted(T).f,
            );

            const r1 = self.radii[index_of_first_greater -| 1];
            const r2 = self.radii[index_of_first_greater];

            const em1 = self.em[index_of_first_greater -| 1];
            const em2 = self.em[index_of_first_greater];

            const t1 = self.time[index_of_first_greater -| 1];
            const t2 = self.time[index_of_first_greater];

            const g1 = self.g[index_of_first_greater -| 1];
            const g2 = self.g[index_of_first_greater];

            const ang1 = self.local_theta[index_of_first_greater -| 1];
            const ang2 = self.local_theta[index_of_first_greater];

            const w = interpolations.lerpWeight(T, r, r1, r2);
            const interp_em = interpolations.lerpValue(T, w, em1, em2);
            const interp_time = interpolations.lerpValue(T, w, t1, t2);
            const interp_g = interpolations.lerpValue(T, w, g1, g2);
            const interp_local_theta = interpolations.lerpValue(T, w, ang1, ang2);

            return .{
                .em = interp_em,
                .t = interp_time,
                .g = interp_g,
                .local_theta = interp_local_theta,
            };
        }

        pub fn init(
            radii: []const T,
            em: []const T,
            time: []const T,
            energyshift: []const T,
            local_theta: []const T,
        ) Self {
            return .{
                .radii = radii,
                .em = em,
                .time = time,
                .g = energyshift,
                .local_theta = local_theta,
            };
        }

        pub fn deinit(self: Self, allocator: std.mem.Allocator) void {
            allocator.free(self.time);
            allocator.free(self.radii);
            allocator.free(self.em);
            allocator.free(self.g);
            allocator.free(self.local_theta);
        }

        /// Serialise to a FITS table.
        pub fn toFITS(self: *const Self, allocator: std.mem.Allocator) !zfits.Hdu {
            var hdu = zfits.Hdu.init(allocator, .{ .binary_table = .empty });
            errdefer hdu.deinit();

            try hdu.setName(
                "EMISSIVITY",
                "This is an axisymmetric emissivity table",
            );
            try utils.addKerrzFITSInfo(&hdu);

            // Setup the columns
            try hdu.data.binary_table.appendColumn(.{
                .label = "radius",
                .comment = "The radial coordinate on the accretion disc",
                .units = "rg",
                .units_comment = "Gravitational radii rg = GM/c^2",
            });
            try hdu.data.binary_table.appendColumn(.{
                .label = "em",
                .comment = "The emissivity value itself",
                .units = "arbitrary",
            });
            try hdu.data.binary_table.appendColumn(.{
                .label = "time",
                .comment = "Corona-to-disc light-travel time",
                .units = "tg",
                .units_comment = "Light-crossing time tg = GM/c^3",
            });
            try hdu.data.binary_table.appendColumn(.{
                .label = "energyshift",
                .comment = "Corona-to-disc energyshift",
            });
            try hdu.data.binary_table.appendColumn(
                .{
                    .label = "local_theta",
                    .comment = "The incident angle off disc surface normal.",
                    .units = "rad",
                },
            );

            for (
                self.radii,
                self.em,
                self.time,
                self.g,
                self.local_theta,
            ) |r, em, t, g, local_theta| {
                const row = try hdu.data.binary_table.addRow();
                row.cols[0].one.float_32 = @floatCast(r.x);
                row.cols[1].one.float_32 = @floatCast(em.x);
                row.cols[2].one.float_32 = @floatCast(t.x);
                row.cols[3].one.float_32 = @floatCast(g.x);
                row.cols[4].one.float_32 = @floatCast(local_theta.x);
            }

            return hdu;
        }

        /// Load from a FITS table. Copies the memory in the table and casts to
        /// the appropriate number type.
        pub fn fromFITS(allocator: std.mem.Allocator, hdu: *zfits.Hdu) !Self {
            const table = &hdu.data.binary_table;

            // TODO: assert that this is actually an emissivity table.
            const radii = try allocator.alloc(T, table.num_rows);
            errdefer allocator.free(radii);

            const em = try allocator.alloc(T, table.num_rows);
            errdefer allocator.free(em);

            const time = try allocator.alloc(T, table.num_rows);
            errdefer allocator.free(time);

            const energyshift = try allocator.alloc(T, table.num_rows);
            errdefer allocator.free(energyshift);

            const local_theta = try allocator.alloc(T, table.num_rows);
            errdefer allocator.free(local_theta);

            for (0..table.num_rows) |i| {
                const row = try table.getOrParseRow(allocator, i);

                radii[i] = .promote(@floatCast(row.cols[0].one.float_32));
                em[i] = .promote(@floatCast(row.cols[1].one.float_32));
                time[i] = .promote(@floatCast(row.cols[2].one.float_32));
                energyshift[i] = .promote(@floatCast(row.cols[3].one.float_32));
                local_theta[i] = .promote(@floatCast(row.cols[4].one.float_32));
            }

            return .{
                .radii = radii,
                .em = em,
                .time = time,
                .g = energyshift,
                .local_theta = local_theta,
            };
        }
    };
}

fn PartitionSorted(comptime T: type) type {
    return struct {
        pub fn f(self: T, other: T) bool {
            return self.x > other.x;
        }
    };
}

pub fn TableEmissivity(comptime T: type) type {
    return struct {
        const Self = @This();

        /// The azimuthal coordinate grid. These are the columns of the
        /// emissivity table.
        phi_grid: []const T,
        /// The radial coordinate grid. These are the rows of the emissivity
        /// table.
        r_grid: []const T,

        /// The emissivity values themselves.
        em: Matrix(T),

        /// The corona-to-disc energyshifts.
        g: Matrix(T),

        /// The corona-to-disc time.
        time: Matrix(T),

        /// The incident angle, measured from the disc normal.
        local_theta: Matrix(T),

        /// Obtain the emissivity averaged over a particular radius.
        pub fn radialValues(self: *const Self, r: T) EmissivityRow(T) {
            _ = self;
            _ = r;
            unreachable; // TODO: implement me
        }

        const Interp = dinterp.HypercubeInterpolator(2, T);

        /// Obtain the emissivity values at a particular point on the disc.
        pub fn values(self: *const Self, r: T, phi: T) EmissivityRow(T) {
            const interp = Interp.interpolate(
                .fromSlices(.{ self.r_grid, self.phi_grid }),
                .{ r, phi },
            );
            var out: EmissivityRow(T) = .{ .em = 0, .t = 0, .g = 0, .local_theta = 0 };
            interp.applyContext(EmissivityRow(T), &out, self, applyInterpolation);

            // Positivity preserving:
            out.em = @max(out.em, 0);
            out.t = @max(out.t, 0);
            out.g = @max(out.g, 0);
            out.local_theta = @max(out.local_theta, 0);
            return out;
        }

        fn applyInterpolation(out: *EmissivityRow(T), op: dinterp.Operation(T), self: *const Self) void {
            out.em += self.em.values[op.index] * op.weight;
            out.t += self.time.values[op.index] * op.weight;
            out.g += self.g.values[op.index] * op.weight;
            out.local_theta += self.local_theta.values[op.index] * op.weight;
        }

        pub fn deinit(self: *const Self, allocator: std.mem.Allocator) void {
            allocator.free(self.phi_grid);
            allocator.free(self.r_grid);
            self.em.deinit(allocator);
            self.time.deinit(allocator);
            self.g.deinit(allocator);
            self.local_theta.deinit(allocator);
        }

        /// Interpolate all missing values along the azimuthal coordinate:
        pub fn interpolateMissing(self: *Self) void {
            for (0..self.r_grid.len) |i| {
                interpolations.interpolateZeroes(T, self.em.getColumn(i));
            }
            for (0..self.r_grid.len) |i| {
                interpolations.interpolateZeroes(T, self.time.getColumn(i));
            }
            for (0..self.r_grid.len) |i| {
                interpolations.interpolateZeroes(T, self.g.getColumn(i));
            }
            for (0..self.r_grid.len) |i| {
                interpolations.interpolateZeroes(T, self.local_theta.getColumn(i));
            }
        }

        /// Write the emissivity and time tables to a FITS HDUs. This creates
        /// three new HDUs. The first is for the r and phi axes. The second is
        /// the emissivity over r and phi, where the columns correpond to
        /// different `r` and the rows to different `phi`. The third is the
        /// same but for the time.
        pub fn toFITS(self: *const Self, allocator: std.mem.Allocator) !zfits.Hdu {
            // Create the HDU for the emissivity and time data:
            var data_hdu = zfits.Hdu.init(allocator, .{ .binary_table = .empty });
            errdefer data_hdu.deinit();
            try data_hdu.setName(
                "EM_TABLE",
                "A disc emissivity table",
            );
            try utils.addKerrzFITSInfo(&data_hdu);

            try data_hdu.data.binary_table.appendColumn(
                .{
                    .label = "radius",
                    .comment = "The radial coordinate on the accretion disc",
                    .units = "rg",
                    .units_comment = "Gravitational radius rg = GM/c^2",
                },
            );
            try data_hdu.data.binary_table.appendColumn(
                .{
                    .label = "em",
                    .col_type = .{ .repeat = self.phi_grid.len },
                    .comment = "The emissivity value itself",
                    .units = "arbitrary",
                },
            );
            try data_hdu.data.binary_table.appendColumn(
                .{
                    .label = "time",
                    .col_type = .{ .repeat = self.phi_grid.len },
                    .comment = "Corona-to-disc light-travel time",
                    .units = "tg",
                    .units_comment = "Light-crossing time tg = GM/c^3",
                },
            );
            try data_hdu.data.binary_table.appendColumn(
                .{
                    .label = "energyshift",
                    .col_type = .{ .repeat = self.phi_grid.len },
                    .comment = "Corona-to-disc energyshift",
                },
            );
            try data_hdu.data.binary_table.appendColumn(
                .{
                    .label = "local_theta",
                    .comment = "The incident angle off disc surface normal.",
                    .units = "rad",
                },
            );

            for (0..self.r_grid.len) |i| {
                const row = try data_hdu.data.binary_table.addRow();

                row.cols[0].one.float_32 = @floatCast(self.r_grid[i]);
                for (row.cols[1].many, 0..self.phi_grid.len) |*value, j| {
                    value.float_32 = @floatCast(self.em.get(j, i));
                }
                for (row.cols[2].many, 0..self.phi_grid.len) |*value, j| {
                    value.float_32 = @floatCast(self.time.get(j, i));
                }
                for (row.cols[3].many, 0..self.phi_grid.len) |*value, j| {
                    value.float_32 = @floatCast(self.g.get(j, i));
                }
                for (row.cols[4].many, 0..self.phi_grid.len) |*value, j| {
                    value.float_32 = @floatCast(self.local_theta.get(j, i));
                }
            }

            return data_hdu;
        }
    };
}

/// Uniform interface to many different emissivity profiles.
pub fn EmissivityProfile(comptime T: type) type {
    return union(enum) {
        const Self = @This();

        powerlaw: PowerLaw(T),
        axisymmetric: AxisymmetricEmissivity(T),
        table: TableEmissivity(T),

        /// Get the emissivity at a paricular annular radius on the accretion
        /// disc, averaging over the azimuthal coordinate if needed.
        pub fn radialValues(self: Self, r: T) EmissivityRow(T) {
            return switch (self) {
                inline else => |p| p.radialValues(r),
            };
        }

        pub fn deinit(self: *const Self, allocator: std.mem.Allocator) void {
            switch (self.*) {
                .powerlaw => {},
                .axisymmetric => self.axisymmetric.deinit(allocator),
                .table => self.table.deinit(allocator),
            }
        }
    };
}

/// A union type for the parameters of all coronal models.
pub fn CoronalModelOptions(comptime T: type) type {
    return union(enum) {
        const Self = @This();
        lamppost: LamppostOptions(T),
        ring: RingOptions(T),
        disc: DiscOptions(T),
        umbrella: UmbrellaOptions(T),

        /// Initialise a `CoronaModel` from the parameters.
        pub fn toCoronalModel(
            self: Self,
            allocator: std.mem.Allocator,
            metric: KerrMetric(T),
        ) !CoronalModel(T) {
            return switch (self) {
                .lamppost => |opts| .{ .lamppost = Lamppost(T).init(metric, opts) },
                .ring => |opts| .{ .ring = Ring(T).init(metric, opts) },
                .disc => |opts| .{ .disc = try Disc(T).init(allocator, metric, opts) },
                .umbrella => |opts| .{ .umbrella = try Umbrella(T).init(allocator, metric, opts) },
            };
        }

        /// Add metadata about the coronal model to a FITS HDU header.
        pub fn addToHdu(self: Self, hdu: *zfits.Hdu) !void {
            switch (self) {
                .lamppost => |lp| {
                    try hdu.appendHeaderRecord(
                        "KZMODEL",
                        .{
                            .value = .{ .string = "LAMPPOST" },
                            .comment = "This is a kerrz lamppost corona table",
                        },
                    );
                    try hdu.appendHeaderRecord("LPHEIGHT", .{
                        .value = .{
                            .float = @floatCast(lp.height.x),
                        },
                        .comment = "Lamppost height in rg",
                    });
                    try hdu.appendHeaderRecord("LPRVEL", .{
                        .value = .{
                            .float = @floatCast(lp.radial_velocity.x),
                        },
                        .comment = "Lamppost radial velocity in units of c",
                    });
                    try hdu.appendHeaderRecord("LPGAMMA", .{
                        .value = .{
                            .float = @floatCast(lp.photon_index.x),
                        },
                        .comment = "Lamppost powerlaw photon index",
                    });
                },
                .ring => |ring| {
                    try hdu.appendHeaderRecord(
                        "KZMODEL",
                        .{
                            .value = .{ .string = "RING" },
                            .comment = "This is a kerrz ring corona table",
                        },
                    );
                    try hdu.appendHeaderRecord("RRADIUS", .{
                        .value = .{
                            .float = @floatCast(ring.radius.x),
                        },
                        .comment = "Ring radius in rg",
                    });
                    try hdu.appendHeaderRecord("RHEIGHT", .{
                        .value = .{
                            .float = @floatCast(ring.height.x),
                        },
                        .comment = "Ring height in rg",
                    });
                    try hdu.appendHeaderRecord("RGAMMA", .{
                        .value = .{
                            .float = @floatCast(ring.photon_index.x),
                        },
                        .comment = "Ring powerlaw photon index",
                    });
                    try hdu.appendHeaderRecord("RVELPROF", .{
                        .value = .{
                            .string = @tagName(ring.velocity),
                        },
                        .comment = "The velocity profile of the ring corona",
                    });
                },
                .disc => |disc| {
                    try hdu.appendHeaderRecord(
                        "KZMODEL",
                        .{
                            .value = .{ .string = "DISC" },
                            .comment = "This is a disc-like coronal model.",
                        },
                    );
                    try hdu.appendHeaderRecord("DHEIGHT", .{
                        .value = .{
                            .float = @floatCast(disc.height.x),
                        },
                        .comment = "The height of the disc (rg).",
                    });
                    try hdu.appendHeaderRecord("DRIN", .{
                        .value = .{
                            .float = @floatCast(disc.height.x),
                        },
                        .comment = "The inner radius of the disc (rg).",
                    });
                    try hdu.appendHeaderRecord("DROUT", .{
                        .value = .{
                            .float = @floatCast(disc.height.x),
                        },
                        .comment = "The outer radius of the disc (rg).",
                    });
                    try hdu.appendHeaderRecord("DGAMMA", .{
                        .value = .{
                            .float = @floatCast(disc.photon_index.x),
                        },
                        .comment = "The disc powerlaw photon index",
                    });
                    try hdu.appendHeaderRecord("DVELPROF", .{
                        .value = .{
                            .string = @tagName(disc.velocity),
                        },
                        .comment = "The velocity profile of the ring corona",
                    });
                },
                .umbrella => |umbrella| {
                    try hdu.appendHeaderRecord(
                        "KZMODEL",
                        .{
                            .value = .{ .string = "UMBRELLA" },
                            .comment = "This is an umbrella coronal model",
                        },
                    );
                    try hdu.appendHeaderRecord("UOFFR", .{
                        .value = .{
                            .float = @floatCast(umbrella.offset_radius.x),
                        },
                        .comment = "The offset radius of the corona (rg)",
                    });
                    try hdu.appendHeaderRecord("UINNER", .{
                        .value = .{
                            .float = @floatCast(umbrella.inner_opening_angle.x),
                        },
                        .comment = "The inner opening angle in radians.",
                    });
                    try hdu.appendHeaderRecord("UANGLE", .{
                        .value = .{
                            .float = @floatCast(umbrella.outer_opening_angle.x),
                        },
                        .comment = "The outer opening angle in radians.",
                    });
                    try hdu.appendHeaderRecord("UGAMMA", .{
                        .value = .{
                            .float = @floatCast(umbrella.photon_index.x),
                        },
                        .comment = "The umbrella powerlaw photon index",
                    });
                    try hdu.appendHeaderRecord("UVELPROF", .{
                        .value = .{
                            .string = @tagName(umbrella.velocity),
                        },
                        .comment = "The velocity profile of the umbrella",
                    });
                },
            }
        }
    };
}

/// A union type for all coronal models.
pub fn CoronalModel(comptime T: type) type {
    return union(enum) {
        const Self = @This();
        lamppost: Lamppost(T),
        ring: Ring(T),
        disc: Disc(T),
        umbrella: Umbrella(T),

        pub fn deinit(self: *const Self, allocator: std.mem.Allocator) void {
            switch (self.*) {
                .disc => |*disc| disc.deinit(allocator),
                .umbrella => |*umbrella| umbrella.deinit(allocator),
                .lamppost, .ring => {},
            }
        }
    };
}

/// A utility method for serialising many ring-like coronae to a FITS HDUs.
///
/// The first HDU is where the emissivity data will be written. This does not
/// include the radius on the disc, as it is assumed that is a repeated column,
/// and so the caller must write the radii themselves.
///
/// The second HDU is where the parameters of the ring-like corona will be
/// written. It must be created using `makeRingInfoHdu()`.
///
/// The allocator must be an arena'd allocator, as this function is leaky when
/// generating the labels for the FITS table.
pub fn writeRingToFITS(
    comptime T: type,
    leaky_allocator: std.mem.Allocator,
    corona: Ring(T),
    photon_fractions: PhotonFractions(T.T),
    em_table: AxisymmetricEmissivity(T),
    emissivity_hdu: *zfits.Hdu,
    ring_info_hdu: *zfits.Hdu,
) !void {
    // The index for this ring in the serialisation:
    const index = ring_info_hdu.data.binary_table.num_rows;

    const info_row = try ring_info_hdu.data.binary_table.addRow();
    info_row.cols[0].one.float_32 = @floatCast(corona.getRadius().x);
    info_row.cols[1].one.float_32 = @floatCast(corona.getHeight().x);
    info_row.cols[2].one.float_32 = @floatCast(corona.photon_index.x);
    // Photon fractions:
    info_row.cols[3].one.float_32 = @floatCast(photon_fractions.infinity);
    info_row.cols[4].one.float_32 = @floatCast(photon_fractions.event_horizon);
    info_row.cols[5].one.float_32 = @floatCast(photon_fractions.disc);
    info_row.cols[6].one.float_32 = @floatCast(photon_fractions.above_isco);
    info_row.cols[7].one.float_32 = @floatCast(photon_fractions.below_isco);
    info_row.cols[8].one.float_32 = @floatCast(photon_fractions.no_status);

    const num_cols = emissivity_hdu.data.binary_table.column_headers.len;

    const em_label = try std.fmt.allocPrint(
        leaky_allocator,
        "r{d}em",
        .{index},
    );

    try emissivity_hdu.data.binary_table.appendColumn(.{
        .label = em_label,
        .comment = "Radial coord on disc for nth ring",
    });

    const t_label = try std.fmt.allocPrint(
        leaky_allocator,
        "r{d}t",
        .{index},
    );
    try emissivity_hdu.data.binary_table.appendColumn(.{
        .label = t_label,
        .comment = "Time coord to disc for nth ring",
    });

    // Write the emissivity data itself.
    var row_itt = emissivity_hdu.data.binary_table.cached_rows.iterator();
    for (
        em_table.em,
        em_table.time,
    ) |e, t| {
        const row = row_itt.next().?.value_ptr;
        row.cols[num_cols].one.float_32 = @floatCast(e.x);
        row.cols[num_cols + 1].one.float_32 = @floatCast(t.x);
    }
}

/// Create and format a HDU that can be used as the `ring_info_hdu` argument in
/// `writeRingToFITS`.
pub fn makeRingInfoHdu(allocator: std.mem.Allocator) !zfits.Hdu {
    var hdu = zfits.Hdu.init(allocator, .{ .binary_table = .empty });
    try hdu.setName("RINGDATA", "Parameters of the constituent ring coronae.");

    const table = &hdu.data.binary_table;

    try table.appendColumn(.{
        .label = "radius",
        .comment = "The radius on the accretion disc",
        .units = "rg",
        .units_comment = "In units of rg = GM/c^2",
    });
    try table.appendColumn(.{
        .label = "height",
        .comment = "The height of the ring",
        .units = "rg",
    });
    try table.appendColumn(.{
        .label = "gamma",
        .comment = "Photon index E^-Gamma",
    });
    try table.appendColumn(.{
        .label = "pf_inf",
        .comment = "Photon fraction at infinity",
    });
    try table.appendColumn(.{
        .label = "pf_bh",
        .comment = "Photon fraction in black hole",
    });
    try table.appendColumn(.{
        .label = "pf_disc",
        .comment = "Photon fraction on disc",
    });
    try table.appendColumn(.{
        .label = "pf_aisco",
        .comment = "Photon fraction above isco on disc",
    });
    try table.appendColumn(.{
        .label = "pf_bisco",
        .comment = "Photon fraction below isco on disc",
    });
    try table.appendColumn(.{
        .label = "pf_none",
        .comment = "Unknown status photon fraction",
    });

    return hdu;
}

test "problematic cases" {
    const Dual = DualNumber(f64, 1);
    const Lamp = Lamppost(Dual);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    {
        const lp = Lamp.init(metric, .{ .height = .promote(4.0) });
        const em = lp.traceEmissivity(
            .equatorial_plane,
            .promote(0.000003141598936787666),
            .{},
        );
        try std.testing.expectEqual(
            geodesic.Status.event_horizon,
            em.status,
        );
    }

    {
        const lp = Lamp.init(metric, .{ .height = .promote(19.7) });
        const em = lp.traceEmissivity(
            .equatorial_plane,
            .promote(0.000003141598936787666),
            .{},
        );
        try std.testing.expectEqual(
            geodesic.Status.event_horizon,
            em.status,
        );
    }
}
