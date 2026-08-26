/// Maps different physical values to pixels.
const std = @import("std");

const antiderivatives = @import("antiderivatives.zig");
const potentials = @import("potentials.zig");
const geodesic = @import("geodesic.zig");
const geometry = @import("geometry.zig");
const redshift = @import("redshift.zig");
const orbits = @import("orbits.zig");
const emissivity = @import("emissivity.zig");

const NullGeodesic = geodesic.NullGeodesic;
const TraceResult = geodesic.TraceResult;
const KerrMetric = geometry.KerrMetric;
const FourVector = geometry.FourVector;

const Mapper = @This();

const MapperValue = struct {
    name: [:0]const u8,
    description: []const u8,
    derivatives: bool = true,
    operator: bool = true,
};

fn makeMapperEnumeration() type {
    comptime var enum_fields: []const std.builtin.Type.EnumField = &.{};
    for (Fields, 0..) |f, i| {
        enum_fields = enum_fields ++ .{std.builtin.Type.EnumField{
            .name = f.name,
            .value = i,
        }};
    }
    const _enum = std.builtin.Type.Enum{
        .decls = &.{},
        .fields = enum_fields,
        .is_exhaustive = true,
        .tag_type = u8,
    };
    return @Type(.{ .@"enum" = _enum });
}

pub const Value = makeMapperEnumeration();

value: Value,
operator: ?Operator = null,
derivative: ?Derivative = null,

/// Get the units of this Mapper.
pub fn getUnits(self: Mapper) ?[]const u8 {
    return switch (self.value) {
        .radius => "rg",
        .time => "tg",
        else => null,
    };
}

/// Get a short descriptor of the Mapper.
pub fn getShortDescriptor(self: Mapper) ?[]const u8 {
    return switch (self.value) {
        else => null,
    };
}

/// Convert a string to a Mapper instance.
pub fn fromString(s: []const u8) !Mapper {
    // split the operations
    var itt = std.mem.tokenizeScalar(u8, s, '.');

    const s1 = itt.next().?;
    const v = std.meta.stringToEnum(Value, s1) orelse
        return error.NoSuchMapper;

    var operator: ?Operator = null;
    var deriv: ?Derivative = null;

    while (itt.next()) |op| {
        if (op.len == 2 and op[0] == 'd') {
            if (deriv != null) return error.InvalidOperator;
            deriv = switch (op[1]) {
                'x' => .x,
                'y' => .y,
                else => return error.InvalidDerivative,
            };
        } else {
            if (operator != null) return error.InvalidOperator;
            operator = std.meta.stringToEnum(Operator, op) orelse
                return error.NoSuchOperator;
        }
    }

    if (operator != null) {
        // check we can actually apply an operator
        for (Fields) |field| {
            if (std.mem.eql(u8, @tagName(v), field.name)) {
                if (field.operator == false) {
                    return error.InvalidOperator;
                }
            }
        }
    }

    if (deriv != null) {
        // check we can actually take a derivative of this value
        for (Fields) |field| {
            if (std.mem.eql(u8, @tagName(v), field.name)) {
                if (field.derivatives == false) {
                    return error.InvalidDerivative;
                }
            }
        }
    }

    return .{ .value = v, .operator = operator, .derivative = deriv };
}

// Returns true if this mapper value requires derivative information.
pub fn needsDerivatives(self: Mapper) bool {
    return self.derivative != null;
}

pub fn Context(comptime T: type) type {
    return struct {
        geometry: KerrMetric(T),
        v_source: FourVector(T),
        em: emissivity.EmissivityProfile(T),
        total: ?antiderivatives.TotalAntiderivatives(T) = null,
    };
}

/// Evaluate the mapper with a particular geodesic.
pub fn calculate(
    self: Mapper,
    comptime T: type,
    ctx: Context(T),
    geod: NullGeodesic(T),
    result: TraceResult(T),
) T.T {
    var ctx_copy = ctx;
    return self.calculateAlt(T, &ctx_copy, geod, result);
}

/// Like `calculate` but can modify the context to cache e.g. the total
/// antiderivative computation. This is useful if calling multiple mappers on
/// the same result.
pub fn calculateAlt(
    self: Mapper,
    comptime T: type,
    ctx: *Context(T),
    geod: NullGeodesic(T),
    result: TraceResult(T),
) T.T {
    const A = T.Algebra;
    switch (self.value) {
        .angular_case, .status => {},
        else => {
            switch (result.status) {
                .event_horizon => return 0,
                .infinity => return 0,
                else => {},
            }
        },
    }

    const dual_value: T = b: {
        switch (self.value) {
            .ang_mom => break :b geod.L,
            .carter => break :b geod.Q,
            .redshift => {
                break :b redshift.keplerianRedshiftResult(
                    T,
                    ctx.geometry,
                    geod,
                    result,
                    ctx.v_source,
                );
            },
            .radius => {
                break :b result.r;
            },
            .theta => {
                break :b result.theta;
            },
            .time => {
                ctx.total = ctx.total orelse result.totalAntiderivatives(ctx.geometry, geod);
                break :b ctx.total.?.coordinateTime(ctx.geometry, geod);
            },
            .delta_time => {
                ctx.total = ctx.total orelse result.totalAntiderivatives(ctx.geometry, geod);
                break :b A.sub(ctx.total.?.coordinateTime(ctx.geometry, geod), geod.x_init.r);
            },
            .azimuth => {
                ctx.total = ctx.total orelse result.totalAntiderivatives(ctx.geometry, geod);
                break :b ctx.total.?.coordinateAzimuth(ctx.geometry, geod);
            },
            .mino => {
                break :b result.mino_time;
            },
            .winding => {
                return @floatFromInt(1 + result.winding);
            },
            .flux => {
                const eshift = redshift.keplerianRedshiftResult(
                    T,
                    ctx.geometry,
                    geod,
                    result,
                    ctx.v_source,
                );
                const eshift_cubed = A.powi(eshift, 3);
                break :b A.mult(eshift_cubed, ctx.em.radialValues(result.r).em);
            },
            .emissivity => {
                break :b ctx.em.radialValues(result.r).em;
            },
            .case => {
                return @floatFromInt(1 + @intFromEnum(result.state.radial_case));
            },
            .angular_case => {
                return @floatFromInt(1 + @intFromEnum(result.state.angular_case));
            },
            .status => {
                return @floatFromInt(1 + @intFromEnum(result.status));
            },
            .r_sign => {
                return @floatFromInt(geod.radial_sign);
            },
            .init_theta_sign => {
                return @floatFromInt(geod.theta_sign);
            },
            .theta_sign => {
                return result.final_theta_sign;
            },
            .local_theta, .local_phi => {
                const projected_r = A.mult(result.r, A.sin(result.theta));
                const v_medium = orbits.keplerianPlungingAlt(T, ctx.geometry, projected_r);
                const angles = result.localAngles(ctx.geometry, geod, v_medium);
                const v = switch (self.value) {
                    .local_theta => angles.theta,
                    .local_phi => angles.phi,
                    else => unreachable,
                };
                break :b v;
            },
            .pot_X => {
                break :b result.potential_X;
            },
            inline .I_0, .I_1, .I_2, .I_plus, .I_minus, .am_X => |field| {
                ctx.total = ctx.total orelse result.totalAntiderivatives(ctx.geometry, geod);
                break :b @field(ctx.total.?.radial, @tagName(field));
            },
            inline .G_phi_init, .G_phi_final, .G_t_init, .G_t_final, .G_phi_half, .G_t_half => |field| {
                ctx.total = ctx.total orelse result.totalAntiderivatives(ctx.geometry, geod);
                break :b @field(ctx.total.?.angular, @tagName(field));
            },
            inline .u_plus, .u_minus => |field| {
                break :b @field(result.state.angular_roots, @tagName(field));
            },
            inline .r1, .r2, .r3, .r4 => |field| {
                break :b @field(result.state.radial_roots, @tagName(field)).arg;
            },
            .theta_0, .theta_1 => {
                switch (result.state.angular_case) {
                    .normal => {
                        if (self.value == .theta_0) {
                            break :b A.acos(A.sqrt(result.state.angular_roots.u_plus));
                        } else {
                            break :b A.acos(A.sqrt(result.state.angular_roots.u_plus).neg());
                        }
                    },
                    .vortical => {
                        if (self.value == .theta_0) {
                            break :b A.acos(A.sqrt(result.state.angular_roots.u_minus));
                        } else {
                            break :b A.acos(A.sqrt(result.state.angular_roots.u_plus));
                        }
                    },
                }
            },
            .fr, .ftheta, .fphi => {
                const v_medium = orbits.circularFourVelocity(T, ctx.geometry, result.r);
                const pol_vec = result.polarisationVector(ctx.geometry, geod, v_medium);
                const v = switch (self.value) {
                    .fr => pol_vec.r,
                    .ftheta => pol_vec.th,
                    .fphi => pol_vec.ph,
                    else => unreachable,
                };
                break :b v;
            },
            .pol_angle => {
                const v_medium = orbits.keplerianPlungingAlt(T, ctx.geometry, result.r);
                const pol_ang = result.polarisationAngle(ctx.geometry, geod, v_medium);
                break :b pol_ang;
            },
            .pol_degree => {
                const v_medium = orbits.keplerianPlungingAlt(T, ctx.geometry, result.r);
                const pol_ang = result.polarisationDegree(ctx.geometry, geod, v_medium);
                break :b pol_ang;
            },
            .stokes_x => {
                const v_medium = orbits.keplerianPlungingAlt(T, ctx.geometry, result.r);
                const xy = result.polarisationXYAtInfinity(ctx.geometry, geod, v_medium);
                break :b xy.x;
            },
            .stokes_y => {
                const v_medium = orbits.keplerianPlungingAlt(T, ctx.geometry, result.r);
                const xy = result.polarisationXYAtInfinity(ctx.geometry, geod, v_medium);
                break :b xy.y;
            },
            .theta_0_time, .theta_1_time => {
                const times = result.state.minoTimeToAngularTurns(
                    0,
                    @floatFromInt(geod.theta_sign),
                );
                if (self.value == .theta_0_time) {
                    break :b times.tau_0;
                } else {
                    break :b times.tau_1;
                }
            },
            .r_potential => {
                break :b potentials.radialFromRoots(T, result.state.radial_roots, result.r);
            },
            .elliptic_sc, .elliptic_sn, .elliptic_cn, .elliptic_dn => |elliptic| {
                if (result.state.angular_case == .vortical) {
                    unreachable;
                }
                const integrals = result.rv.integrals;
                const v = switch (elliptic) {
                    .elliptic_sc => integrals.sc,
                    .elliptic_sn => integrals.sn,
                    .elliptic_cn => integrals.cn,
                    .elliptic_dn => integrals.dn,
                    else => unreachable,
                };
                break :b v;
            },
            .vt, .vr, .vtheta, .vphi => {
                const velocity = result.velocity(ctx.geometry, geod);
                const v = switch (self.value) {
                    .vt => velocity.t,
                    .vr => velocity.r,
                    .vtheta => velocity.th,
                    .vphi => velocity.ph,
                    else => unreachable,
                };
                break :b v;
            },
            .G_theta => {
                const angular_values = antiderivatives.angularFromCache(
                    T,
                    result.state.angular_case,
                    result.state.angular_cache,
                    result.theta,
                );
                break :b angular_values.G_theta_final;
            },
            .init_G_theta => {
                break :b result.state.angular_cache.G_theta_init;
            },
            .G_theta_half => {
                break :b result.state.angular_cache.G_theta_half_libration;
            },
            .u_ratio => {
                break :b A.div(
                    result.state.angular_roots.u_plus,
                    result.state.angular_roots.u_minus,
                );
            },
        }
    };

    var value = b: {
        if (self.derivative) |d| {
            if (T.N < 2) {
                unreachable;
            }
            const index: usize = if (d == .x) 0 else 1;
            break :b dual_value.dx[index];
        }
        break :b dual_value.x;
    };

    if (self.operator) |op| {
        value = switch (op) {
            .log => std.math.log10(@abs(value)),
            .cos => @cos(value),
            .sqrt => std.math.sqrt(@abs(value)),
            .sign => std.math.sign(value),
            .inv => 1 / value,
            .mod2pi => @mod(value, std.math.pi * 2.0),
            .abs => @abs(value),
        };
    }

    return value;
}

pub const Operator = enum {
    log,
    cos,
    sqrt,
    sign,
    inv,
    mod2pi,
    abs,
};

pub const Derivative = enum {
    x,
    y,
};

pub const Fields = [_]MapperValue{
    .{
        .name = "ang_mom",
        .description = "The angular momentum constant of motion.",
    },
    .{
        .name = "carter",
        .description = "Carter's constant.",
    },
    .{
        .name = "redshift",
        .description = "The energyshift along the geodesic.",
    },
    .{
        .name = "radius",
        .description = "The radial coordinate at the endpoint of the geodesic.",
    },
    .{
        .name = "theta",
        .description = "The poloidal coordinate at the endpoint of the geodesic.",
    },
    .{
        .name = "time",
        .description = "The time coordinate.",
    },
    .{
        .name = "delta_time",
        .description = "The time coordinate with the observer radius subtracted.",
    },
    .{
        .name = "azimuth",
        .description = "The azimuthal coordinate.",
    },
    .{
        .name = "case",
        .description = "The radial case (I, II, III, IV) determined by the radial roots.",
        .derivatives = false,
        .operator = false,
    },
    .{
        .name = "angular_case",
        .description = "The angular case (normal or vortical).",
        .derivatives = false,
        .operator = false,
    },
    .{
        .name = "status",
        .description = "The status code of the geodesic.",
        .derivatives = false,
        .operator = false,
    },
    .{
        .name = "mino",
        .description = "The Mino time. This is equivalently I_0.",
    },
    .{
        .name = "winding",
        .description = "The winding number of the geodesic.",
        .derivatives = false,
        .operator = false,
    },
    .{
        .name = "flux",
        .description = "The flux, in arbitrary units.",
    },
    .{
        .name = "emissivity",
        .description = "The emissivity at a particular point on the disc, in arbitrary units.",
    },
    .{
        .name = "local_theta",
        .description = "The cosine of the local elevation angle of the photon's four-velocity in the rest frame of the disc.",
    },
    .{
        .name = "local_phi",
        .description = "The cosine of the local azimuthal angle of the photon's four-velocity in the rest frame of the disc.",
    },
    .{
        .name = "r_potential",
        .description = "The value of the radial potential R(r) at the endpoint of the geodesic.",
    },
    .{
        .name = "r_sign",
        .description = "The sign of the radial momentum.",
        .derivatives = false,
        .operator = false,
    },
    .{
        .name = "theta_sign",
        .description = "The sign of the poloidal momentum at the endpoint of the geodesic.",
        .derivatives = false,
        .operator = false,
    },
    .{
        .name = "init_theta_sign",
        .description = "The sign of the poloidal momentum at the startpoint of the geodesic.",
        .derivatives = false,
        .operator = false,
    },
    .{
        .name = "vt",
        .description = "The time component of the geodesic's four-velocity.",
    },
    .{
        .name = "vr",
        .description = "The radial component of the geodesic's four-velocity.",
    },
    .{
        .name = "vtheta",
        .description = "The theta component of the geodesic's four-velocity.",
    },
    .{
        .name = "vphi",
        .description = "The theta component of the geodesic's four-velocity.",
    },
    .{
        .name = "I_0",
        .description = "The radial antiderivative I_0.",
    },
    .{
        .name = "I_1",
        .description = "The radial antiderivative I_1.",
    },
    .{
        .name = "I_2",
        .description = "The radial antiderivative I_2.",
    },
    .{
        .name = "I_plus",
        .description = "The radial antiderivative I_plus.",
    },
    .{
        .name = "I_minus",
        .description = "The radial antiderivative I_minus.",
    },
    .{
        .name = "G_phi_init",
        .description = "The angular antiderivative G_phi at the start.",
    },
    .{
        .name = "G_phi_final",
        .description = "The angular antiderivative G_phi at the endpoint.",
    },
    .{
        .name = "G_phi_half",
        .description = "The angular antiderivative G_phi over a half-libration.",
    },
    .{
        .name = "G_t_init",
        .description = "The angular antiderivative G_t at the start.",
    },
    .{
        .name = "G_t_final",
        .description = "The angular antiderivative G_t at the endpoint.",
    },
    .{
        .name = "G_t_half",
        .description = "The angular antiderivative G_t over a half-libration.",
    },
    .{
        .name = "G_theta",
        .description = "The angular antiderivative G_theta at the endpoint of the geodesic.",
    },
    .{
        .name = "init_G_theta",
        .description = "The angular antiderivative G_theta at the start of the geodesic.",
    },
    .{
        .name = "G_theta_half",
        .description = "The half-libration value of the G_theta antiderivative.",
    },
    .{
        .name = "am_X",
        .description = "The ampltiude of the radial X potential.",
    },
    .{
        .name = "pot_X",
        .description = "The value of the radial X potential.",
    },
    .{
        .name = "u_plus",
        .description = "The positive angular root.",
    },
    .{
        .name = "u_minus",
        .description = "The negative angular root.",
    },
    .{
        .name = "r1",
        .description = "The (magnitude of the) first radial root, as it may be a complex number.",
    },
    .{
        .name = "r2",
        .description = "The (magnitude of the) second radial root.",
    },
    .{
        .name = "r3",
        .description = "The (magnitude of the) third radial root.",
    },
    .{
        .name = "r4",
        .description = "The (magnitude of the) fourth radial root.",
    },
    .{
        .name = "pol_angle",
        .description = "The polarisation angle.",
    },
    .{
        .name = "pol_degree",
        .description = "The polarisation degree.",
    },
    .{
        .name = "stokes_x",
        .description = "The normalised Stokes parameter X.",
    },
    .{
        .name = "stokes_y",
        .description = "The normalised Stokes parameter Y.",
    },
    .{
        .name = "fr",
        .description = "The radial component of the polarisation vector.",
    },
    .{
        .name = "ftheta",
        .description = "The elevation component of the polarisation vector.",
    },
    .{
        .name = "fphi",
        .description = "The azimuthal component of the polarisation vector.",
    },
    .{
        .name = "elliptic_sc",
        .description = "The sc Jacobi elliptic integral evaluated at the endpoint.",
    },
    .{
        .name = "elliptic_sn",
        .description = "The sn Jacobi elliptic integral.",
    },
    .{
        .name = "elliptic_cn",
        .description = "The cn Jacobi elliptic integral.",
    },
    .{
        .name = "elliptic_dn",
        .description = "The dn Jacobi elliptic integral.",
    },
    .{
        .name = "theta_0_time",
        .description = "The Mino time at the first polar turning point.",
    },
    .{
        .name = "theta_1_time",
        .description = "The Mino time at the second polar turning point.",
    },
    .{
        .name = "theta_0",
        .description = "The Mino time at the first polar turning point.",
    },
    .{
        .name = "theta_1",
        .description = "The Mino time at the second polar turning point.",
    },
    .{
        .name = "u_ratio",
        .description = "The ratio of u_plus to u_minus",
    },
};
