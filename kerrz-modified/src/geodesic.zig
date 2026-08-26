/// References:
/// - G&L: Gralla and Lupsasca, 2020
/// - B72: Bardeen et al, 1972
const std = @import("std");
const ad = @import("zad");
const dinterp = @import("dinterp");
const options = @import("options");

const geometry = @import("geometry.zig");
const potentials = @import("potentials.zig");
const antiderivatives = @import("antiderivatives.zig");
const accretion_discs = @import("accretion-discs.zig");
const orbits = @import("orbits.zig");
const polarisation = @import("polarisation.zig");
const ComplexNumber = @import("complex.zig").ComplexNumber;
const integration_state = @import("integration-state.zig");

const IntegrationState = integration_state.IntegrationState;
const KerrMetric = geometry.KerrMetric;
const AccretionDisc = accretion_discs.AccretionDisc;
const FourVector = geometry.FourVector;
const ThreeVector = geometry.ThreeVector;
const DualNumber = ad.DualNumber;
const AnglePair = geometry.AnglePair;

const TEST_TOLERANCE = options.test_numerical_tolerance;

test "sample geodesic" {
    const Dual = ad.DualNumber(f64, 1);
    const geom: KerrMetric(Dual) = .init(.one, .promote(0.998));
    const x: FourVector(Dual) = .{
        .t = .zero,
        .r = .promote(1e8),
        .th = .promote(std.math.degreesToRadians(50.0)),
        .ph = .zero,
    };

    const geod = NullGeodesic(Dual).fromImpactParameters(geom, x, .promote(1.1), .promote(0.5));
    try std.testing.expectApproxEqAbs(0.3384179922859939, geod.eta.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-0.8426488874308758, geod.lambda.x, TEST_TOLERANCE);

    const roots = potentials.rootsOfRadialPotential(Dual, geom, geod.eta, geod.lambda);

    try std.testing.expect(roots.r1.isReal());
    try std.testing.expectApproxEqAbs(-1.9769897749725187, roots.r1.real().x, TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(0.04524048760408694, roots.r2.real().x, TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(0.9658746436842156, roots.r3.real().x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-1.6839567463476899, roots.r3.imag().x, TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(0.9658746436842156, roots.r4.real().x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(1.6839567463476899, roots.r4.imag().x, TEST_TOLERANCE);
}

pub fn ImpactParameters(comptime T: type) type {
    return struct {
        alpha: T,
        beta: T,
    };
}

/// Map the angles on the local sky of a point `x_src` moving with velocity
/// `v_src` onto an intial velocity vector for a null-geodesic. When
/// `sky_theta` is zero, the vector is pointing towards the zenith, and when
/// `sky_phi` is `pi`, the vector is pointing towards the symmetry-axis,
/// equivalently the z-axis.
pub fn skyAnglesToVelocity(
    comptime T: type,
    ts: KerrMetric(T).TangentSpace,
    v_src: FourVector(T),
    sky_theta: T,
    sky_phi: T,
) FourVector(T) {
    // TODO: this is an expensive calculation that doesn't need to be done
    // again and again for the same v_src
    return skyAnglesToVelocityFrame(
        T,
        ts.x,
        ts.localFrame(v_src),
        sky_theta,
        sky_phi,
    );
}

/// Alternative signature for `skyAnglesToVelocity`.
pub fn skyAnglesToVelocityFrame(
    comptime T: type,
    x: FourVector(T),
    frame: geometry.TetradFrame(T),
    sky_theta: T,
    sky_phi: T,
) FourVector(T) {
    const transform = geometry.cartesianFromSpherical(
        T,
        x.th,
        x.ph,
    );
    // This implementation differs from Gradus.jl in that it does not have the
    // -1 multiplication to be consistent with other parts of the code.
    const v_xfm = transform.apply(.init(.one, sky_theta, sky_phi));
    // Promote to a FourVector with a dummy time variable
    const v_fourvec: FourVector(T) = .{
        .t = .one,
        .r = v_xfm.v[0],
        .th = v_xfm.v[1],
        .ph = v_xfm.v[2],
    };

    return frame.apply(v_fourvec);
}

test "skyAnglesToVelocity" {
    const Dual = ad.DualNumber(f64, 0);
    const FV = FourVector(Dual);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    const x: FV = .{
        .t = .zero,
        .r = .promote(10.0),
        .th = .promote(std.math.degreesToRadians(20)),
        .ph = .zero,
    };

    const v_stationary: FV = .{
        .t = .promote(1.1168175629588084),
        .r = .zero,
        .th = .zero,
        .ph = .zero,
    };

    const ts = metric.tangentSpace(x);

    const result_s = skyAnglesToVelocity(
        Dual,
        ts,
        v_stationary,
        .promote(std.math.degreesToRadians(20)),
        .promote(std.math.degreesToRadians(30)),
    );

    try std.testing.expectApproxEqAbs(1.115381473629473, result_s.t.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(0.8820032115685071, result_s.r.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-0.004287049627681557, result_s.th.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(0.04974574883067166, result_s.ph.x, TEST_TOLERANCE);

    const v: FV = .{
        .t = .promote(0.3),
        .r = .promote(0.1),
        .th = .zero,
        .ph = .promote(-0.3),
    };

    const result = skyAnglesToVelocity(
        Dual,
        ts,
        v,
        .promote(std.math.degreesToRadians(20)),
        .promote(std.math.degreesToRadians(30)),
    );

    try std.testing.expectApproxEqAbs(1.5717228271574222, result.t.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(1.4348252841784677, result.r.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-0.004287049627681566, result.th.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-0.3356958408658669, result.ph.x, TEST_TOLERANCE);

    // Make sure all the constraints work too
    const constrained = ts.constrainVector(result, 0);
    try std.testing.expectApproxEqAbs(2.214584272022076, constrained.t.x, TEST_TOLERANCE);

    const lowered = ts.lowerIndices(constrained);
    try std.testing.expectApproxEqAbs(-1.76775913324798, lowered.t.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(1.7870566573299376, lowered.r.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-0.4324753957562234, lowered.th.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-4.018171479439075, lowered.ph.x, TEST_TOLERANCE);
}

/// Rotate a vector `v` by an `angle` about the axis `k`. Uses Rodrigues'
/// rotation formula.
fn rotateAboutVector(
    comptime T: type,
    v: ThreeVector(T),
    angle: T,
    k: ThreeVector(T),
) ThreeVector(T) {
    const A = T.Algebra;
    const cos_angle = A.cos(angle);
    const sin_angle = A.sin(angle);

    const term_1 = v.scalarMult(cos_angle);
    const term_2 = k.cross(v).scalarMult(sin_angle);
    const term_3 = k.scalarMult(A.mult(A.sub(.one, cos_angle), k.dot(v)));

    return term_1.add(term_2).add(term_3);
}

test "rotateAboutVector" {
    const Dual = ad.DualNumber(f64, 0);
    const v1: ThreeVector(Dual) = .init(.promote(2.0), .promote(-3.0), .promote(0.5));
    const v2: ThreeVector(Dual) = .init(.promote(-9.3), .promote(1.0), .promote(2.5));

    const new = rotateAboutVector(
        Dual,
        v1,
        .promote(std.math.degreesToRadians(30)),
        v2,
    );

    try std.testing.expectApproxEqAbs(31.087413014344932, new.v[0].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-0.499459244339989, new.v[1].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(6.567055119425537, new.v[2].x, TEST_TOLERANCE);
}

pub fn SkyAnglesRotatedOptions(comptime T: type) type {
    return struct {
        theta_0: T = .zero,
        phi_0: T = .zero,
    };
}

/// The same as `skyAnglesToVelocity`, but offsets the origin of the sky
/// coordinates by angles `theta_0, phi_0`. That is, the entire coordinate
/// system of the sky is reorientated so that `theta = 0, phi = 0` points in
/// the direction `theta_0, phi_0`.
pub fn skyAnglesToVelocityRotated(
    comptime T: type,
    ts: KerrMetric(T).TangentSpace,
    v_src: FourVector(T),
    sky_theta: T,
    sky_phi: T,
    opts: SkyAnglesRotatedOptions(T),
) FourVector(T) {
    const A = T.Algebra;
    const Xfm = geometry.CartesianFromSpherical(T);
    const InvXfm = geometry.SphericalFromCartesian(T);

    // The axis about which to rotate.
    const k = Xfm.transformAlt(
        .init(.one, opts.theta_0, opts.phi_0),
    );

    // The offset vector corresponding to the given elevation.
    const q = Xfm.transformAlt(
        .init(.one, A.add(sky_theta, opts.theta_0), opts.phi_0),
    );

    const rotated = rotateAboutVector(T, q, sky_phi, k);
    // Transform back to spherical vector to get the new angles.
    const spherical = InvXfm.transformAlt(rotated);
    return skyAnglesToVelocity(T, ts, v_src, spherical.v[1], spherical.v[2]);
}

/// Translate a velocity vector at some point `x_src` that is moving with
/// velocity `v_src` onto the local sky at `x_src`, returning the elevation and
/// azimuthal angles `theta, phi`. This mapping performs the special
/// relativistic Lorentz transformation to account for boosting of the local
/// frame. It is the inverse of `skyAnglesToVelocity`.
pub fn velocityToSkyAngles(
    comptime T: type,
    ts: KerrMetric(T).TangentSpace,
    v_src: FourVector(T),
    velocity: FourVector(T),
) AnglePair(T) {
    const basis = ts.localBasis(v_src);
    const local_velocity = basis.apply(velocity);
    const transform = geometry.sphericalFromCartesian(T, ts.x.th, ts.x.ph);
    // drop the time component
    const cartesian_velocity = transform.apply(.{ .v = [_]T{
        local_velocity.r,
        local_velocity.th,
        local_velocity.ph,
    } });
    // negate the phi component for consistency with `skyAnglesToVelocity`.
    return .{ .theta = cartesian_velocity.v[1], .phi = cartesian_velocity.v[2].neg() };
}

test velocityToSkyAngles {
    const Dual = ad.DualNumber(f64, 0);
    const FV = FourVector(Dual);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    const x: FV = .{
        .t = .zero,
        .r = .promote(10.0),
        .th = .promote(std.math.degreesToRadians(20)),
        .ph = .zero,
    };

    const v_stationary: FV = .{
        .t = .promote(1.1168175629588084),
        .r = .zero,
        .th = .zero,
        .ph = .zero,
    };

    const v: FV = .{
        .t = .promote(1.115381473629473),
        .r = .promote(0.8820032115685071),
        .th = .promote(-0.004287049627681557),
        .ph = .promote(0.04974574883067166),
    };

    const ts = metric.tangentSpace(x);

    const angles = velocityToSkyAngles(
        Dual,
        ts,
        v_stationary,
        v,
    );

    try std.testing.expectApproxEqAbs(std.math.degreesToRadians(20), angles.theta.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(std.math.degreesToRadians(30), angles.phi.x, TEST_TOLERANCE);
}

fn reverseVelocity(
    comptime T: type,
    ts: KerrMetric(T).TangentSpace,
    v: FourVector(T),
    magnitude: T.T,
) FourVector(T) {
    const A = T.Algebra;
    // The reversal strategy is simply to rotate the angles of the velocity
    // vector:
    var out = v;
    out.th = A.sub(.promote(std.math.pi), v.th);
    out.ph = A.sub(v.ph, .promote(std.math.pi));
    return ts.constrainVector(out, magnitude);
}

test "velocity reversal" {
    const Dual = ad.DualNumber(f64, 0);
    const FV = FourVector(Dual);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    const x: FV = .{
        .t = .zero,
        .r = .promote(10.0),
        .th = .promote(std.math.degreesToRadians(20)),
        .ph = .zero,
    };

    // Choosing some arbitrary vector:
    const v_unconstrained: FV = .{
        .t = .zero,
        .r = .promote(0.5),
        .th = .promote(0.2),
        .ph = .promote(0.1),
    };

    const ts = metric.tangentSpace(x);
    const v_null = ts.constrainVector(v_unconstrained, 0);

    const v_rev = reverseVelocity(Dual, ts, v_null, 0);

    // Check the normalisation
    try std.testing.expectApproxEqAbs(
        v_null.properNorm(ts).x,
        v_rev.properNorm(ts).x,
        TEST_TOLERANCE,
    );

    // Then compute the same vector the long-way-round:
    const cartFromSpher = geometry.cartesianFromSpherical(Dual, x.th, x.ph);
    const spherFromCart = geometry.sphericalFromCartesian(Dual, x.th, x.ph);

    const v_null_cart = cartFromSpher.apply(v_null.toThreeVector());

    // Reversing the direction in cartesian is really just inverting all of the
    // coordinates:
    const v_null_cart_reverse = v_null_cart.scalarMult(.promote(-1));

    const v_null_spher = spherFromCart.apply(v_null_cart_reverse);

    const w_unconstrained = v_null_spher.toFourVector(.zero);
    const w = ts.constrainVector(w_unconstrained, 0);

    try std.testing.expectApproxEqAbs(
        w.r.x,
        v_rev.r.x,
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        w.th.x,
        v_rev.th.x,
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        w.ph.x,
        v_rev.ph.x,
        TEST_TOLERANCE,
    );
}

pub const TracingConfig = struct {
    /// The number of photon windings this geodesic should follow.
    winding: usize = 0,

    /// The radial windings.
    r_winding: usize = 0,

    /// For algorithms that need a step size on the Mino time, this controls
    /// what that step size should be. Currently, this is only used for the
    /// geometrically thick disc intersection calculations.
    ///
    /// This has to currently be quite small, as the algorithms can go awry if
    /// the root solving bracketing interval goes over the black hole. Within
    /// that interval, there are invalid points which cannot currently be
    /// walked over.
    mino_step: f64 = 0.01,

    /// Controls what an acceptable error tolerance for this geodesic is. This
    /// is also currently only used in the (absolute) position error when
    /// tracing geometrically thick disc intersections.
    error_tolerance: f64 = 1e-5,
};

/// The ultimate status associated with a geodesic.
pub const Status = enum(u8) {
    /// May be in free space or otherwise some unknown status.
    no_status,
    /// The geodesic connects to the event horizon.
    event_horizon,
    /// The geodesic goes to infinity.
    infinity,
    /// The geodesic intersected with the accretion disc.
    intersected_disc,

    /// Update the status field of the result.
    /// Determine the status of the geodesic based on information about
    /// the trace and potentials.
    fn determine(
        comptime T: type,
        result: TraceResult(T),
        geod: NullGeodesic(T),
        metric: KerrMetric(T),
    ) Status {
        const r_sign: T.T = @floatFromInt(geod.radial_sign);

        if (result.mino_time.x <= 0 and result.r.x >= 0) {
            if (result.state.angular_case == .vortical) {
                return .event_horizon;
            }
            return .infinity;
        }

        if (result.mino_time.x > 0 and result.r.x <= 0) {
            return .infinity;
        }

        if (result.r.x < metric.horizon_radius.x * 1.01) {
            return .event_horizon;
        }

        switch (result.state.radial_case) {
            .case_II => {
                const integrals = result.rv.integrals;
                // This protects against erroneous false-images.
                if (std.math.sign(integrals.cn.x) < 0 and std.math.sign(integrals.sn.x) > 0) {
                    return .infinity;
                }
            },
            .case_III, .case_IV => {
                if (r_sign < 0 and result.rv.X.x > 0) {
                    if (result.state.angular_case != .vortical) {
                        return .event_horizon;
                    }
                }

                if (r_sign > 0) {
                    const half_mino = @min(1, result.mino_time.x / 2.0);
                    const half_rv = result.state.radiusAtMinoTime(.promote(half_mino), r_sign);

                    if (half_rv.r.x < 0) {
                        return .infinity;
                    }
                }
            },
            else => {},
        }

        return .no_status;
    }
};

pub fn TraceResult(comptime T: type) type {
    const A = T.Algebra;

    return struct {
        const Geodesic = NullGeodesic(T);
        const Result = @This();

        /// The determined Mino time.
        mino_time: T,
        /// The radial coordinate this geodesic connects to.
        r: T,
        /// The angular coordinate this geodesic connects to.
        theta: T,
        /// The status of the geodesic.
        status: Status,
        /// The state of the integration. This is a value cache for further
        /// computations, but also stores things related to e.g. the radial
        /// case.
        state: IntegrationState(T),
        /// The values of the angular potentials.
        rv: antiderivatives.PrincipalRadialValues(T),
        /// The winding number
        winding: usize,
        /// The final value of the X potential parameter.
        potential_X: T,
        /// The sign of the poloidal momentum at this point in the geodesic.
        final_theta_sign: T.T,

        pub fn adapt(self: Result, comptime NewT: type) TraceResult(NewT) {
            return .{
                .mino_time = .adaptFrom(self.mino_time),
                .r = .adaptFrom(self.r),
                .theta = .adaptFrom(self.theta),
                .status = self.status,
                .state = self.state.adapt(NewT),
                .rv = self.rv.adapt(NewT),
                .winding = self.winding,
                .potential_X = .adaptFrom(self.potential_X),
                .final_theta_sign = self.final_theta_sign,
            };
        }

        /// Compute the total anti-derivatives given the trace result.
        pub fn totalAntiderivatives(
            self: Result,
            metric: KerrMetric(T),
            geod: Geodesic,
        ) antiderivatives.TotalAntiderivatives(T) {
            return antiderivatives.totalAntiderivatives(
                T,
                metric,
                self.state.radial_case,
                self.state.angular_case,
                self.state.radial_cache,
                self.state.radial_roots,
                self.state.angular_cache,
                self.state.angular_roots,
                self.rv,
                geod.x_init.r,
                self.mino_time,
                geod.x_init.th,
                self.theta,
                // These signs are always the initial signs, and all subsequent
                // sign changes are inferred from the winding number.
                @floatFromInt(geod.radial_sign),
                @floatFromInt(geod.theta_sign),
                self.final_theta_sign,
                self.winding,
            );
        }

        /// Calculates the four-position of the photon at the endpoint of the
        /// geodesic in the global coordinates.
        ///
        /// This function is inefficient if other result quantities that depend
        /// on the antiderivatives are needed, but is a useful short-hand for
        /// quickly obtaining the position.
        pub fn position(
            self: Result,
            metric: KerrMetric(T),
            geod: Geodesic,
        ) FourVector(T) {
            const total = self.totalAntiderivatives(metric, geod);
            const t = total.coordinateTime(metric, geod);
            const ph = total.coordinateAzimuth(metric, geod);
            return .{
                .t = t,
                .r = self.r,
                .th = self.theta,
                .ph = ph,
            };
        }

        /// Calculates the four-velocity of the photon at the endpoint of
        /// the geodesic in the global coordinates.
        pub fn velocity(
            self: Result,
            metric: KerrMetric(T),
            geod: Geodesic,
        ) FourVector(T) {
            var vel = fourVelocity(T, metric, geod.eta, geod.lambda, self.r, self.theta);
            // get the turning point information to correct the sign of the velocity
            const integrals = self.rv.integrals;
            vel.r = A.mult(.promote(std.math.sign(integrals.sn.x)), vel.r);
            vel.th = A.mult(.promote(std.math.sign(self.final_theta_sign)), vel.th);
            return vel;
        }

        /// Map the endpoint of the geodesic to the local angles of a
        /// medium at `x_medium == (0, result.r, result.theta, 0)` with
        /// velocity `v_medium`.
        pub fn localAngles(
            self: Result,
            metric: KerrMetric(T),
            geod: Geodesic,
            v_medium: FourVector(T),
        ) AnglePair(T) {
            // TODO: this function should use pre-calculated versions of
            // the geodesic velocity and the tangent space to avoid
            // exessive recalculations.
            const ts = metric.tangentSpaceAlt(self.r, self.theta);
            return velocityToSkyAngles(T, ts, v_medium, self.velocity(metric, geod));
        }

        /// Same as `localAngles` but makes the angles point in the antipodal
        /// direction.
        pub fn localAnglesReverse(
            self: Result,
            metric: KerrMetric(T),
            geod: Geodesic,
            v_medium: FourVector(T),
        ) AnglePair(T) {
            // TODO: this function should use pre-calculated versions of
            // the geodesic velocity and the tangent space to avoid
            // exessive recalculations.
            var angles = self.localAngles(metric, geod, v_medium);
            angles.theta.x = @abs(std.math.pi - angles.theta.x);
            angles.phi.x = @abs(angles.phi.x - std.math.pi);
            return angles;
        }

        /// Calculate the polarisation degree at infinity with impact
        /// parameters `alpha` and `beta`.
        pub fn polarisationDegree(
            self: Result,
            metric: KerrMetric(T),
            geod: Geodesic,
            v_medium: FourVector(T),
        ) T {
            return self.polarisationXYAtInfinity(metric, geod, v_medium).degree();
        }

        /// Calculate the polarisation angle at infinity with impact
        /// parameters `alpha` and `beta`.
        pub fn polarisationAngle(
            self: Result,
            metric: KerrMetric(T),
            geod: Geodesic,
            v_medium: FourVector(T),
        ) T {
            return self.polarisationXYAtInfinity(metric, geod, v_medium).angle();
        }

        /// Compute the X and Y normalised Stokes parameters as seen by an
        /// observer at infinity.
        pub fn polarisationXYAtInfinity(
            self: Result,
            metric: KerrMetric(T),
            geod: Geodesic,
            v_medium: FourVector(T),
        ) polarisation.PolarisationConstant(T).XY {
            const pc = self.polarisationConstant(metric, geod, v_medium);
            const ip = geod.calculateImpactParameters(metric, geod.x_init.th);
            return pc.polarisationXYAtInfinityAlt(
                metric,
                geod.x_init.th,
                ip.alpha,
                ip.beta,
            );
        }

        /// Calculate the polarisation constant from the result.
        pub fn polarisationConstant(
            self: Result,
            metric: KerrMetric(T),
            geod: Geodesic,
            v_medium: FourVector(T),
        ) polarisation.PolarisationConstant(T) {
            const photon_velocity = self.velocity(metric, geod);
            const ts = metric.tangentSpaceAlt(self.r, self.theta);
            return calculatePolarisationConstant(
                T,
                metric,
                ts,
                v_medium,
                photon_velocity,
                self.r,
                self.theta,
            );
        }

        /// Calculate the polarisation vector along the geodesic.
        ///
        /// This function calculates a polarisation vector `f` at the endpoint
        /// of the geodesic, along with the polarisation constant. It then uses
        /// the polarisation constant to calculate the polarisation vector at
        /// the start point of the geodesic. This 'translated' vector is then
        /// returned.
        pub fn polarisationVector(
            self: Result,
            metric: KerrMetric(T),
            geod: Geodesic,
            v_medium: FourVector(T),
        ) FourVector(T) {
            // evaluate endpoint
            const photon_velocity = self.velocity(metric, geod);
            const ts = metric.tangentSpaceAlt(self.r, self.theta);
            const pc = calculatePolarisationConstant(
                T,
                metric,
                ts,
                v_medium,
                photon_velocity,
                self.r,
                self.theta,
            );
            // evaluate startpoint
            const initial_photon_velocity = fourVelocity(
                T,
                metric,
                geod.eta,
                geod.lambda,
                geod.x_init.r,
                geod.x_init.th,
            );
            return polarisation.vectorFromPolarisationConstant(
                T,
                metric,
                geod.x_init.r,
                geod.x_init.th,
                initial_photon_velocity,
                pc,
            );
        }

        /// Extrapolate the solution (assuming asymptotically flat space) to
        /// infinity, and return the `theta` and `phi` coordinates.
        ///
        /// The `phi` argument is the azimuthal coordinate of the result. It
        /// must be calculated by first calling `.totalAntiderivatives()`, and
        /// is passed in so that the result may be reused.
        pub fn extrapolateAsymptoticRay(
            self: Result,
            metric: KerrMetric(T),
            geod: Geodesic,
            phi: T,
        ) AnglePair(T) {
            // Convert the tangent velocity vector to cartesian coordinates,
            // and then project out to infinity
            return geometry.projectRayToInfinitySpherical(
                T,
                .init(.one, self.theta, phi),
                self.velocity(metric, geod).toThreeVector(),
            );
        }
    };
}

// T is the number type of the geodesic
pub fn NullGeodesic(comptime T: type) type {
    const A = T.Algebra;

    return struct {
        const Result = TraceResult(T);
        const Self = @This();

        // constants of motion
        E: T,
        L: T,
        Q: T,

        // This is L / E
        lambda: T,
        // This is Q / E^2
        eta: T,

        // Initial position of the geodesic
        x_init: FourVector(T),

        theta_sign: i32 = 1,
        radial_sign: i32 = -1,
        windings: usize = 0,

        /// Adapt the dual number type
        pub fn adapt(self: Self, comptime NewT: type) NullGeodesic(NewT) {
            return .{
                .E = .adaptFrom(self.E),
                .L = .adaptFrom(self.L),
                .Q = .adaptFrom(self.Q),
                .lambda = .adaptFrom(self.lambda),
                .eta = .adaptFrom(self.eta),
                .x_init = self.x_init.adapt(NewT),
                .theta_sign = self.theta_sign,
                .radial_sign = self.radial_sign,
                .windings = self.windings,
            };
        }

        /// Initialise a geodesic from the constants of motion. Note that this
        /// implicitly sets the sign of the poloidal momentum to +ive. Use the
        /// `theta_sign` field to adjust as needed.
        pub fn fromConstantsOfMotion(x: FourVector(T), E: T, L: T, Q: T) Self {
            return .{
                .E = E,
                .L = L,
                .Q = Q,
                .lambda = A.div(L, E),
                .eta = A.div(Q, A.powi(E, 2)),
                .x_init = x,
            };
        }

        /// Initialise the geodesic at a point `x` with velocity `v`. The
        /// velocity must already be constrained to be a null-geodesic (use
        /// `KerrMetric.constrainVector`). This is not checked by this
        /// initialiser. Use `fromVelocity` if the vector is not guaranteed to
        /// be null-like.
        pub fn fromNullVelocity(
            metric: KerrMetric(T),
            x: FourVector(T),
            v: FourVector(T),
        ) Self {
            return fromNullVelocityTangentSpace(metric, metric.tangentSpace(x), v);
        }

        /// Same as `fromNullVelocity` but with a pre-computed tangent space at
        /// `x`.
        pub fn fromNullVelocityTangentSpace(
            metric: KerrMetric(T),
            ts: KerrMetric(T).TangentSpace,
            v: FourVector(T),
        ) Self {
            return fromMomenta(metric, ts.x, ts.lowerIndices(v));
        }

        /// Initialise the geodesic at a point `x` with an (unconstrained)
        /// velocity `v`. If the velocity is already constrained to be
        /// null-like, use `fromNullVelocity` instead.
        pub fn fromVelocity(
            metric: KerrMetric(T),
            x: FourVector(T),
            v: FourVector(T),
        ) Self {
            return fromVelocityTangentSpace(metric, metric.tangentSpace(x), v);
        }

        /// Same as `fromVelocity` but with a pre-computed tangent space at
        /// `x`.
        pub fn fromVelocityTangentSpace(
            metric: KerrMetric(T),
            ts: KerrMetric(T).TangentSpace,
            v: FourVector(T),
        ) Self {
            return fromNullVelocityTangentSpace(metric, ts, ts.constrainVector(v, 0));
        }

        /// Initialise a geodesic from angles on the local sky of a point `x`
        /// assuming it is stationary in the spacetime, `v = (1, 0, 0, 0)` (not
        /// LNRF!). Calculates the effects of the Lorentz boosting implicitly
        /// into the constants of motion.
        pub fn fromStationarySkyAngles(
            metric: KerrMetric(T),
            x: FourVector(T),
            sky_theta: T,
            sky_phi: T,
        ) Self {
            // TODO: there's a lot of room for optimisation here. For example,
            // the metric components are calculated a number of times, when
            // they only need to be calculated once. That is, they are
            // calculated for the orthonomalisation of the local frame, again
            // for constraining the magnitude of the velocity vector, and then
            // again for lowering the indices of the velocity vector to a
            // momentum vector. That's a lot of wasted computation that could
            // build up on a big render.
            const ts = metric.tangentSpace(x);
            return fromStationarySkyAnglesTangentSpace(metric, ts, sky_theta, sky_phi);
        }

        pub fn fromStationarySkyAnglesTangentSpace(
            metric: KerrMetric(T),
            ts: KerrMetric(T).TangentSpace,
            sky_theta: T,
            sky_phi: T,
        ) Self {
            return fromSkyAnglesTangentSpace(
                metric,
                ts,
                orbits.stationary(T, ts),
                sky_theta,
                sky_phi,
            );
        }

        /// Initialise a geodesic from angles on the local sky of a point `x`
        /// that is moving with velocity `v`.  Calculates the effects of the
        /// Lorentz boosting implicitly into the constants of motion. See also
        /// `fromStationarySkyAngles`.
        pub fn fromSkyAngles(
            metric: KerrMetric(T),
            x: FourVector(T),
            v: FourVector(T),
            sky_theta: T,
            sky_phi: T,
        ) Self {
            const ts = metric.tangentSpace(x);
            return fromSkyAnglesTangentSpace(metric, ts, v, sky_theta, sky_phi);
        }

        pub fn fromSkyAnglesTangentSpace(
            metric: KerrMetric(T),
            ts: KerrMetric(T).TangentSpace,
            v: FourVector(T),
            sky_theta: T,
            sky_phi: T,
        ) Self {
            const unconstrained_velocity = skyAnglesToVelocity(
                T,
                ts,
                v,
                sky_theta,
                sky_phi,
            );
            return fromVelocityTangentSpace(metric, ts, unconstrained_velocity);
        }

        pub fn fromSkyAnglesRotatedTangentSpace(
            metric: KerrMetric(T),
            ts: KerrMetric(T).TangentSpace,
            v: FourVector(T),
            sky_theta: T,
            sky_phi: T,
            opts: SkyAnglesRotatedOptions(T),
        ) Self {
            const unconstrained_velocity = skyAnglesToVelocityRotated(
                T,
                ts,
                v,
                sky_theta,
                sky_phi,
                opts,
            );
            return fromVelocityTangentSpace(metric, ts, unconstrained_velocity);
        }

        /// Calculate the initial constants of motions from a four-momentum
        /// vector at a particular point in the spacetime.
        pub fn fromMomenta(m: KerrMetric(T), x: FourVector(T), p: FourVector(T)) Self {
            var geod = Self.fromConstantsOfMotion(
                x,
                p.t.neg(),
                p.ph,
                m.carterConstant(x, p),
            );
            geod.theta_sign = @intFromFloat(std.math.sign(p.th.x));
            geod.radial_sign = @intFromFloat(std.math.sign(p.r.x));

            // Avoid zeroing a bunch of calculations later down the line.
            if (geod.theta_sign == 0) geod.theta_sign = 1;
            if (geod.radial_sign == 0) geod.radial_sign = 1;

            return geod;
        }

        /// Calculate the initial constant of motions from a four-position and
        /// a set of impact parameters. Note: this currently only uses the
        /// formula for impact parameters at infinity, and not for finite
        /// distance.
        pub fn fromImpactParameters(m: KerrMetric(T), x: FourVector(T), alpha: T, beta: T) Self {
            // TODO: this assumes the observer is at infinity, but they are not by `x`!
            const L = A.mult(alpha.neg(), A.sin(x.th));
            const Q = A.add(
                A.powi(beta, 2),
                A.mult(
                    A.powi(A.cos(x.th), 2),
                    A.sub(
                        A.powi(alpha, 2),
                        A.powi(m.a, 2),
                    ),
                ),
            );
            var geod = Self.fromConstantsOfMotion(x, .one, L, Q);

            geod.theta_sign = @intFromFloat(-std.math.sign(beta.x));
            if (geod.theta_sign == 0) {
                geod.theta_sign = 1;
            }

            return geod;
        }

        /// From the constants of motion, calculate the equivalent impact
        /// parameters at some (infinite distance) observer inclination
        /// `theta`.
        pub fn calculateImpactParameters(self: Self, m: KerrMetric(T), theta: T) ImpactParameters(T) {
            const alpha = A.div(self.lambda, A.sin(theta)).neg();
            const beta_squared = A.sub(
                self.eta,
                A.mult(
                    A.powi(A.cos(theta), 2),
                    A.sub(A.powi(alpha, 2), A.powi(m.a, 2)),
                ),
            );
            return .{
                .alpha = alpha,
                .beta = A.mult(
                    .promote(@floatFromInt(-self.theta_sign)),
                    A.sqrt(beta_squared),
                ),
            };
        }

        /// Trace the null-geodesic conditioned by some accretion disc,
        /// calculating all and any intersection points.
        pub fn traceDisc(
            self: Self,
            m: KerrMetric(T),
            disc: AccretionDisc(T),
            config: TracingConfig,
        ) Result {
            return disc.trace(m, self, config);
        }

        /// Trace the null-geodesic to a a particular Mino time.
        pub fn atMinoTime(
            self: Self,
            m: KerrMetric(T),
            mino_time: T,
            config: TracingConfig,
        ) Result {
            const state = IntegrationState(T).init(m, self);
            const radial_values = state.radiusAtMinoTime(
                mino_time,
                @floatFromInt(self.radial_sign),
            );

            const angle = state.angleAtMinoTime(
                mino_time,
                @floatFromInt(self.theta_sign),
            );

            const theta_sign: T.T = @floatFromInt(self.theta_sign);

            var result: Result = .{
                .mino_time = mino_time,
                .r = radial_values.r,
                .theta = angle.theta,
                .status = .no_status,
                .state = state,
                .rv = radial_values,
                .winding = config.winding,
                .potential_X = radial_values.X,
                .final_theta_sign = theta_sign,
            };
            result.status = Status.determine(T, result, self, m);
            return result;
        }

        /// Trace the null-geodesic to a given radial coordinate specified by
        /// `radius`.
        pub fn traceToRadius(
            self: Self,
            m: KerrMetric(T),
            radius: T,
            config: TracingConfig,
        ) Result {
            var state = IntegrationState(T).init(m, self);
            if (config.r_winding == 1) {
                state.r_winding_sign = -1;
            }

            const mino_time = state.minoTimeToRadius(
                radius,
                @floatFromInt(self.radial_sign),
            );

            const radial_values = state.radiusAtMinoTime(
                mino_time,
                @floatFromInt(self.radial_sign),
            );

            const angle = state.angleAtMinoTime(
                mino_time,
                @floatFromInt(self.theta_sign),
            );

            const theta_sign: T.T = @floatFromInt(self.theta_sign);

            var result: Result = .{
                .mino_time = mino_time,
                .r = radial_values.r,
                .theta = angle.theta,
                .status = .no_status,
                .state = state,
                .rv = radial_values,
                .winding = config.winding,
                .potential_X = radial_values.X,
                .final_theta_sign = theta_sign,
            };
            result.status = Status.determine(T, result, self, m);
            return result;
        }

        /// Trace the null-geodesic to a given theta coordinate specified by
        /// `angle`.
        pub fn traceToAngle(
            self: Self,
            m: KerrMetric(T),
            angle: T,
            config: TracingConfig,
        ) Result {
            const state = IntegrationState(T).init(m, self);

            const mino_time = state.minoTimeToAngle(
                angle,
                config.winding,
                @floatFromInt(self.theta_sign),
            );

            const radial_values = state.radiusAtMinoTime(
                mino_time,
                @floatFromInt(self.radial_sign),
            );

            const theta_sign: T.T = @floatFromInt(self.theta_sign);
            var winding = config.winding;

            if (theta_sign < 0) {
                winding += 1;
            }

            var result: Result = .{
                .mino_time = mino_time,
                .r = radial_values.r,
                .theta = angle,
                .status = .no_status,
                .state = state,
                .rv = radial_values,
                .winding = winding,
                .potential_X = radial_values.X,
                .final_theta_sign = if (@mod(winding, 2) == 1) -theta_sign else theta_sign,
            };
            result.status = Status.determine(T, result, self, m);
            return result;
        }

        pub const PathBuilder = struct {
            state: IntegrationState(T),
            geod: Self,
            metric: KerrMetric(T),
            config: TracingConfig,

            /// The Mino time until the first theta turning point is reached.
            theta_0_time: T,
            /// The Mino time until the second theta turning point is reached.
            theta_1_time: T,
            /// The Mino time until the first r turning point is reached.
            r_0_time: T,
            /// The Mino time until the second r turning point is reached.
            r_1_time: T,

            /// Trace to a particular angle. This is the same as
            /// `NullGeodesic.traceToAngle`.
            pub fn traceToAngle(self: PathBuilder, angle: T) TraceResult(T) {
                const time = self.state.minoTimeToAngle(
                    angle,
                    0,
                    @floatFromInt(self.geod.theta_sign),
                );
                return self.atMinoTime(time);
            }

            /// Trace to a particular radius. This is the same as
            /// `NullGeodesic.traceToRadius`.
            pub fn traceToRadius(self: PathBuilder, radius: T) TraceResult(T) {
                const time = self.state.minoTimeToRadius(
                    radius,
                    @floatFromInt(self.geod.radial_sign),
                );
                return self.atMinoTime(time);
            }

            /// Obtain a `TraceResult` at a particular Mino time.
            pub fn atMinoTime(self: PathBuilder, mino_time: T) Result {
                const radial_values = self.state.radiusAtMinoTime(
                    mino_time,
                    @floatFromInt(self.geod.radial_sign),
                );

                const theta_vals = self.state.angleAtMinoTime(
                    mino_time,
                    @floatFromInt(self.geod.theta_sign),
                );

                // Handle any possible turning points. This doesn't so much
                // matter if all we're interested in is the start and the
                // endpoint, as the approach manifestly seems to flip the
                // position of the observer. But when the full path needs to be
                // reconstructed, it is important the turning points are
                // handled correctly.
                var winding = self.config.winding;

                var theta_sign: T.T = @floatFromInt(self.geod.theta_sign);

                const _K = self.state.angular_cache.G_theta_half_libration.x;

                const mino_time_now = mino_time.x;
                var winding_count: usize = 0;

                if (self.theta_0_time.x > 0 and mino_time_now > self.theta_0_time.x) {
                    winding_count += 1;
                }

                winding_count += @intFromFloat(
                    @floor(@max(0, (mino_time_now - self.theta_0_time.x) / _K)),
                );

                if (self.theta_0_time.x > 0 and self.theta_1_time.x > self.theta_0_time.x) {
                    winding_count += @intFromFloat(
                        @floor(@max(0, (mino_time_now - self.theta_1_time.x) / _K)),
                    );
                }

                winding += winding_count;

                if (@mod(winding, 2) == 1) {
                    theta_sign = -theta_sign;
                }

                winding = @mod(winding, 2);

                var result: Result = .{
                    .mino_time = mino_time,
                    .r = radial_values.r,
                    .theta = theta_vals.theta,
                    .state = self.state,
                    .status = .no_status,
                    .rv = radial_values,
                    .winding = winding,
                    .potential_X = radial_values.X,
                    .final_theta_sign = theta_sign,
                };

                // Do not determine the status if we've not gone anywhere.
                if (mino_time.x != 0) {
                    result.status = Status.determine(T, result, self.geod, self.metric);
                }

                return result;
            }

            fn fromState(
                geod: Self,
                metric: KerrMetric(T),
                state: IntegrationState(T),
                config: TracingConfig,
            ) PathBuilder {
                const angular_turns = state.minoTimeToAngularTurns(
                    config.winding,
                    @floatFromInt(geod.theta_sign),
                );
                const radial_turns = state.minoTimeToRadialTurns(
                    @floatFromInt(geod.radial_sign),
                );

                return .{
                    .state = state,
                    .geod = geod,
                    .metric = metric,
                    .config = config,
                    .theta_0_time = angular_turns.tau_0,
                    .theta_1_time = angular_turns.tau_1,
                    .r_0_time = radial_turns.r_0 orelse .zero,
                    .r_1_time = radial_turns.r_1 orelse .zero,
                };
            }

            /// Same as `NullGeodesic(T).traceBuilder` but uses a pre-computed
            /// integration state from the result of another trace.
            pub fn fromResult(
                geod: Self,
                metric: KerrMetric(T),
                result: Result,
            ) PathBuilder {
                return .fromState(
                    geod,
                    metric,
                    result.state,
                    .{ .winding = result.winding },
                );
            }
        };

        /// Returns a tracer that can be used to evaluate points along the
        /// geodesic as a function of Mino time.
        pub fn traceBuilder(
            self: Self,
            m: KerrMetric(T),
            config: TracingConfig,
        ) PathBuilder {
            const state = IntegrationState(T).init(m, self);
            return .fromState(self, m, state, config);
        }

        /// Calculates the four-velocity of the photon at the startpoint of the
        /// geodesic in the global coordinates.
        pub fn initialVelocity(
            self: Self,
            metric: KerrMetric(T),
        ) FourVector(T) {
            return fourVelocity(
                T,
                metric,
                self.eta,
                self.lambda,
                self.x_init.r,
                self.x_init.th,
            );
        }
    };
}

/// A helper function for calculating the polarisation constant at the endpoint
/// of a null geodesic.
fn calculatePolarisationConstant(
    comptime T: type,
    metric: KerrMetric(T),
    ts: KerrMetric(T).TangentSpace,
    v_medium_backward: FourVector(T),
    photon_velocity: FourVector(T),
    r: T,
    theta: T,
) polarisation.PolarisationConstant(T) {
    // Need to transform the velocity into the 'forward' direction.
    var v_medium = v_medium_backward;
    v_medium.r = v_medium.r.neg();
    const disc_norm = polarisation.normalVectorEquatorial(T, ts, v_medium);

    if (@import("builtin").mode == .Debug) {
        // TODO: this should be removed once the sub-isco velocities are
        // implemented.
        // Check the assumptions:
        const n_magnitude = disc_norm.dot(ts, disc_norm);
        std.debug.assert(std.math.approxEqAbs(T.T, 1.0, n_magnitude.x, TEST_TOLERANCE));
        const k_magnitude = photon_velocity.dot(ts, photon_velocity);
        std.debug.assert(std.math.approxEqAbs(T.T, 0.0, k_magnitude.x, TEST_TOLERANCE));
        const v_magnitude = v_medium.dot(ts, v_medium);
        std.debug.assert(std.math.approxEqAbs(T.T, -1.0, v_magnitude.x, TEST_TOLERANCE));
    }

    const initial_vector = polarisation.initialPolarisationVector(
        T,
        ts,
        disc_norm,
        v_medium,
        photon_velocity,
    );

    if (@import("builtin").mode == .Debug) {
        // verify the properties of the polarisation vector
        const f_magnitude = initial_vector.dot(ts, initial_vector);
        std.debug.assert(std.math.approxEqAbs(T.T, 1.0, f_magnitude.x, TEST_TOLERANCE));
        const f_perp = initial_vector.dot(ts, photon_velocity);
        std.debug.assert(std.math.approxEqAbs(T.T, 0.0, f_perp.x, TEST_TOLERANCE));
        const k_magnitude = photon_velocity.dot(ts, photon_velocity);
        std.debug.assert(std.math.approxEqAbs(T.T, 0.0, k_magnitude.x, TEST_TOLERANCE));
    }

    return polarisation.polarisationConstant(
        T,
        metric,
        r,
        theta,
        photon_velocity,
        initial_vector,
    );
}

test "basic traces" {
    const Dual = ad.DualNumber(f64, 1);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    const x: FourVector(Dual) = .{
        .t = .zero,
        .r = .promote(1e8),
        .th = .promote(std.math.degreesToRadians(50.0)),
        .ph = .zero,
    };
    {
        const geod = NullGeodesic(Dual).fromImpactParameters(
            metric,
            x,
            .promote(-2.16),
            .promote(2.16),
        );
        const result = geod.traceToAngle(metric, .promote(std.math.pi / 2.0), .{});
        try std.testing.expectApproxEqAbs(0.7568040584050688, result.mino_time.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(1.7894030555257678, result.r.x, TEST_TOLERANCE);
    }

    {
        const geod = NullGeodesic(Dual).fromImpactParameters(
            metric,
            x,
            .promote(1.0),
            .promote(-0.5),
        );
        const result = geod.traceToAngle(metric, .promote(std.math.pi / 2.0), .{});
        try std.testing.expectEqual(potentials.RadialCase.case_III, result.state.radial_case);
        try std.testing.expectEqual(Status.event_horizon, result.status);
        try std.testing.expectApproxEqAbs(1.3218288537122953, result.mino_time.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.07184383752334776, result.r.x, TEST_TOLERANCE);
    }
}

test "ray extrapolation" {
    const Dual = ad.DualNumber(f64, 0);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));

    {
        const geod = NullGeodesic(Dual).fromStationarySkyAngles(
            metric,
            .{
                .t = .zero,
                .r = .promote(8.0),
                .th = .promote(1e-5),
                .ph = .zero,
            },
            .promote(std.math.degreesToRadians(20.0)),
            .zero,
        );

        const result = geod.traceToRadius(metric, .promote(1e8), .{});

        const projected = result.extrapolateAsymptoticRay(metric, geod, .zero);
        // At this distance, everything should be consistent with simply using the
        // coordinates.
        try std.testing.expectApproxEqAbs(result.theta.x, projected.theta.x, TEST_TOLERANCE);
    }

    {
        const geod = NullGeodesic(Dual).fromStationarySkyAngles(
            metric,
            .{
                .t = .zero,
                .r = .promote(2.0),
                .th = .promote(1e-5),
                .ph = .zero,
            },
            .promote(std.math.degreesToRadians(20.0)),
            .zero,
        );

        const result = geod.traceToRadius(metric, .promote(1e1), .{});

        const projected = result.extrapolateAsymptoticRay(metric, geod, .zero);
        try std.testing.expectApproxEqAbs(0.6403394476178499, result.theta.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.6243923820022169, projected.theta.x, TEST_TOLERANCE);
    }
}

test "impact parameters back and forth" {
    const Dual = ad.DualNumber(f64, 1);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    const x: FourVector(Dual) = .{
        .t = .zero,
        .r = .promote(1e8),
        .th = .promote(std.math.degreesToRadians(50.0)),
        .ph = .zero,
    };
    {
        const geod = NullGeodesic(Dual).fromImpactParameters(
            metric,
            x,
            .promote(-12.16),
            .promote(2.16),
        );
        const ip = geod.calculateImpactParameters(metric, x.th);
        try std.testing.expectApproxEqAbs(-12.16, ip.alpha.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(2.16, ip.beta.x, TEST_TOLERANCE);
    }
    {
        const geod = NullGeodesic(Dual).fromImpactParameters(
            metric,
            x,
            .promote(0.16),
            .promote(-18.16),
        );
        const ip = geod.calculateImpactParameters(metric, x.th);
        try std.testing.expectApproxEqAbs(0.16, ip.alpha.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-18.16, ip.beta.x, TEST_TOLERANCE);
    }
}

test "traces by mino time" {
    const Dual = ad.DualNumber(f64, 1);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    const x: FourVector(Dual) = .{
        .t = .zero,
        .r = .promote(1e8),
        .th = .promote(std.math.degreesToRadians(50.0)),
        .ph = .zero,
    };

    {
        // Recreate tracing to a equatorial plane.
        const geod = NullGeodesic(Dual).fromImpactParameters(
            metric,
            x,
            .promote(-2.16),
            .promote(2.16),
        );
        const path_builder = geod.traceBuilder(metric, .{});
        const result = path_builder.atMinoTime(.promote(0.7568040584050688));
        try std.testing.expectApproxEqAbs(0.7568040584050688, result.mino_time.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(1.7894030555257678, result.r.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(std.math.pi / 2.0, result.theta.x, TEST_TOLERANCE);
    }
}

test "constants of motion" {
    const Dual = ad.DualNumber(f64, 0);
    const FV = FourVector(Dual);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    {
        const x: FV = .{
            .t = .zero,
            .r = .promote(10.0),
            .th = .promote(std.math.degreesToRadians(20)),
            .ph = .zero,
        };
        const geod = NullGeodesic(Dual).fromStationarySkyAngles(
            metric,
            x,
            .promote(std.math.degreesToRadians(20)),
            .promote(std.math.degreesToRadians(30)),
        );
        try std.testing.expectApproxEqAbs(0.8954013915671949, geod.E.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.5620292966418767, geod.L.x, TEST_TOLERANCE);
    }
    {
        const x: FV = .{
            .t = .zero,
            .r = .promote(10.0),
            .th = .promote(std.math.degreesToRadians(1)),
            .ph = .zero,
        };
        const geod = NullGeodesic(Dual).fromStationarySkyAngles(
            metric,
            x,
            .promote(std.math.degreesToRadians(70)),
            .promote(std.math.degreesToRadians(120)),
        );
        try std.testing.expectApproxEqAbs(0.8955287646699938, geod.E.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.14266565959566255, geod.L.x, TEST_TOLERANCE);
    }
    {
        const x: FV = .{
            .t = .zero,
            .r = .promote(1e4),
            .th = .promote(std.math.degreesToRadians(1)),
            .ph = .zero,
        };
        const geod = NullGeodesic(Dual).fromImpactParameters(
            metric,
            x,
            .zero,
            .promote(-2.0),
        );
        try std.testing.expectApproxEqAbs(-1.21615248827334e-7, geod.lambda.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(3.0042993693628364, geod.eta.x, TEST_TOLERANCE);
    }
}

test "face on" {
    const Dual = ad.DualNumber(f64, 0);
    const FV = FourVector(Dual);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    {
        const x: FV = .{
            .t = .zero,
            .r = .promote(1e4),
            .th = .promote(std.math.degreesToRadians(0)),
            .ph = .zero,
        };

        const geod = NullGeodesic(Dual).fromImpactParameters(
            metric,
            x,
            .zero,
            .promote(-2.0),
        );

        try std.testing.expectApproxEqAbs(0, geod.lambda.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(3.003996, geod.eta.x, TEST_TOLERANCE);

        const radial_roots = potentials.rootsOfAngularPotential(Dual, metric, geod.eta, geod.lambda);
        try std.testing.expectApproxEqAbs(1.0, radial_roots.u_plus.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-3.0160481283207696, radial_roots.u_minus.x, TEST_TOLERANCE);
    }
}

/// Calculates the time and azimuthal velocity from the constants of motion at
/// a particular point in the spacetime.
/// B72, Equations (2.9c) and (2.9d), RHS only, i.e. does not include the Sigma
/// term.
fn velocityTimeAndPhi(
    comptime T: type,
    metric: KerrMetric(T),
    lambda: T,
    r: T,
    theta: T,
) [2]T {
    var dt_dphi: [2]T = undefined;
    const A = T.Algebra;
    const delta = metric.delta(r);
    const sin_theta_squared = A.powi(A.sin(theta), 2);

    const _T = A.sub(A.add(A.powi(r, 2), A.powi(metric.a, 2)), A.mult(metric.a, lambda));
    const phi_term = A.sub(metric.a, A.div(lambda, sin_theta_squared));

    dt_dphi[1] = A.sub(A.div(A.mult(_T, metric.a), delta), phi_term);

    const t_term_1 = A.mult(A.div(_T, delta), A.add(A.powi(r, 2), A.powi(metric.a, 2)));
    const t_term_2 = A.mult(metric.a, A.sub(A.mult(metric.a, sin_theta_squared), lambda));

    dt_dphi[0] = A.sub(t_term_1, t_term_2);

    return dt_dphi;
}

test "velocity time and phi" {
    const Dual = ad.DualNumber(f64, 0);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));

    const dt_dphi = velocityTimeAndPhi(
        Dual,
        metric,
        .promote(0.5),
        .promote(4.0),
        .promote(std.math.degreesToRadians(40)),
    );
    try std.testing.expectApproxEqAbs(31.254996404995804, dt_dphi[0].x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(2.042285341566963, dt_dphi[1].x, TEST_TOLERANCE);
}

/// Compute the four velocity (up to the signs of the r and theta component)
/// from the constants of motion at a particular point in the spacetime.
fn fourVelocity(
    comptime T: type,
    metric: KerrMetric(T),
    eta: T,
    lambda: T,
    r: T,
    theta: T,
) FourVector(T) {
    const A = T.Algebra;
    const sigma = metric.sigma(r, theta);
    const dt_dphi = velocityTimeAndPhi(T, metric, lambda, r, theta);
    // TODO: use the radial roots to speed this computation up a little, and
    // generally cache / share values better.
    const _Vr = potentials.radial(T, metric, eta, lambda, r);
    const _Vth = potentials.angular(T, metric, eta, lambda, theta);

    return .{
        .t = A.div(dt_dphi[0], sigma),
        .r = A.div(A.sqrt(_Vr), sigma),
        .th = A.div(A.sqrt(_Vth), sigma),
        .ph = A.div(dt_dphi[1], sigma),
    };
}

test "four velocities" {
    const Dual = ad.DualNumber(f64, 0);
    const metric: KerrMetric(Dual) = .init(.one, .promote(0.998));
    const v = fourVelocity(
        Dual,
        metric,
        .promote(0.2),
        .promote(0.5),
        .promote(4.0),
        .promote(std.math.degreesToRadians(40)),
    );
    try std.testing.expectApproxEqAbs(1.8845931874964392, v.t.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(0.9873326388047902, v.r.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(0.03951247482924956, v.th.x, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(0.12314437639883175, v.ph.x, TEST_TOLERANCE);
}

test "opposite impact in images" {
    const Dual = DualNumber(f64, 0);
    const V = FourVector(Dual);
    const M = KerrMetric(Dual);
    const G = NullGeodesic(Dual);
    const mid_plane = Dual.promote(std.math.pi / 2.0);

    const metric: M = .init(.one, .promote(0.998));
    const x: V = .{
        .t = .zero,
        .r = .promote(1e4),
        .th = .promote(std.math.degreesToRadians(80)),
        .ph = .zero,
    };

    {
        const geod = G.fromImpactParameters(
            metric,
            x,
            .promote(1),
            .promote(5),
        );
        const result = geod.traceToAngle(metric, mid_plane, .{});
        const total = result.totalAntiderivatives(metric, geod);
        try std.testing.expectApproxEqAbs(
            2.6095718495100977,
            total.coordinateAzimuth(metric, geod).x,
            TEST_TOLERANCE,
        );

        const v = result.velocity(metric, geod);
        try std.testing.expectApproxEqAbs(2.889543633643554, v.t.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.23768040331095439, v.r.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.5702962283823976, v.th.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.09227901411320273, v.ph.x, TEST_TOLERANCE);
    }

    // Now try the opposite impact parameter
    {
        const geod = G.fromImpactParameters(
            metric,
            x,
            .promote(1),
            .promote(-5),
        );
        const result = geod.traceToAngle(metric, mid_plane, .{});
        const total = result.totalAntiderivatives(metric, geod);
        try std.testing.expectApproxEqAbs(
            0.03343569593623379,
            total.coordinateAzimuth(metric, geod).x,
            TEST_TOLERANCE,
        );

        const v = result.velocity(metric, geod);
        try std.testing.expectApproxEqAbs(1.0749787273072524, v.t.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.9859503446777134, v.r.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.006068737264802827, v.th.x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.0011031484214341087, v.ph.x, TEST_TOLERANCE);
    }
}

/// A collection of points along a geodesic as a function of Mino time. This
/// can be used to construct an interpolator to interpolate along geodesic
/// trajectories.
pub fn PathPoints(comptime T: type) type {
    return struct {
        const Self = @This();

        pub const empty: Self = .{};

        pub const Point = struct {
            mino: T.T,
            r: T.T,
            theta: T.T,
        };

        pub const ArrayOfPoints = std.MultiArrayList(Point);

        const SortContext = struct {
            mino_times: []const T.T,
            pub fn lessThan(ctx: SortContext, a_index: usize, b_index: usize) bool {
                return ctx.mino_times[a_index] < ctx.mino_times[b_index];
            }
        };

        knots: ArrayOfPoints = .empty,

        pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
            self.knots.deinit(allocator);
        }

        /// Get an interpolator for this set of path points.
        ///
        /// This will sort the knots inplace before returning the interpolator.
        pub fn interpolator(self: *Self, allocator: std.mem.Allocator) !PathInterpolator(T) {
            const mino_times = self.knots.slice().items(.mino);
            self.knots.sort(SortContext{ .mino_times = mino_times });
            return .initKnots(allocator, self.knots);
        }

        /// Add a new point to the interpolation.
        pub fn addPoint(
            self: *Self,
            allocator: std.mem.Allocator,
            mino: T.T,
            r: T.T,
            theta: T.T,
        ) !void {
            try self.knots.append(allocator, .{
                .mino = mino,
                .r = std.math.log10(r),
                .theta = @cos(theta),
            });
        }
    };
}

/// Interpolate the null-geodesic path from a sparse set of points.
///
/// The interpolator does not support derivatives.
pub fn PathInterpolator(comptime T: type) type {
    return struct {
        const Self = @This();

        pub const Point = PathPoints(T).Point;

        const CubicSpline = dinterp.CubicSplineInterpolator(T.T);

        r_interpolator: CubicSpline,
        theta_interpolator: CubicSpline,

        /// Initialise the interpolator from a set of arrays.
        pub fn init(
            allocator: std.mem.Allocator,
            mino: []const T.T,
            log10_r: []const T.T,
            cos_theta: []const T.T,
        ) !Self {
            const r_interp = try CubicSpline.init(allocator, mino, log10_r);
            errdefer r_interp.deinit(allocator);

            const theta_interp = try CubicSpline.init(allocator, mino, cos_theta);
            errdefer theta_interp.deinit(allocator);

            return .{
                .r_interpolator = r_interp,
                .theta_interpolator = theta_interp,
            };
        }

        /// Initialise the interpolator from a struct-of-arrays of Point.
        pub fn initKnots(
            allocator: std.mem.Allocator,
            knots: PathPoints(T).ArrayOfPoints,
        ) !Self {
            const slices = knots.slice();
            return .init(
                allocator,
                slices.items(.mino),
                slices.items(.r),
                slices.items(.theta),
            );
        }

        pub fn deinit(self: *const Self, allocator: std.mem.Allocator) void {
            self.r_interpolator.deinit(allocator);
            self.theta_interpolator.deinit(allocator);
        }

        /// Interpolate the `r` and `theta` coordinates at a particular Mino
        /// time.
        pub fn interpolate(self: *const Self, mino: T.T) Point {
            return .{
                .mino = mino,
                .r = std.math.pow(
                    T.T,
                    10,
                    self.r_interpolator.interpolate(mino),
                ),
                .theta = std.math.acos(
                    std.math.clamp(self.theta_interpolator.interpolate(mino), -1, 1),
                ),
            };
        }

        /// Force a recalculation of the interpolation weights.
        pub fn updateWeights(
            self: *Self,
            mino: []const T.T,
            log10_r: []const T.T,
            cos_theta: []const T.T,
        ) !void {
            self.r_interpolator.recalculateWeightsFor(
                mino,
                log10_r,
            );
            self.theta_interpolator.recalculateWeightsFor(
                mino,
                cos_theta,
            );

            if (@import("builtin").mode == .Debug) {
                // Some sanity checking:
                for (self.r_interpolator.weights) |weight| {
                    if (!std.math.isFinite(weight)) {
                        return error.BadWeight;
                    }
                }
                for (self.theta_interpolator.weights) |weight| {
                    if (!std.math.isFinite(weight)) {
                        return error.BadWeight;
                    }
                }
            }
        }
    };
}
