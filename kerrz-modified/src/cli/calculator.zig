const std = @import("std");
const clippy = @import("clippy");
const kerrz = @import("kerrz");

const utils = @import("./utils.zig");

const Constants = struct {
    const mass_sun = 1.99e30;
    const big_G = 6.67e-11;
    const speed_of_light = 2.99e8;
};

const LOCAL_TETRAD_NAMES: []const []const u8 = &.{
    "e^μ_(t)",
    "e^μ_(r)",
    "e^μ_(θ)",
    "e^μ_(φ)",
};

const COORDINATE_NAMES: []const []const u8 = &.{
    "t",
    "r",
    "θ",
    "φ",
};

pub const short_description = "A simple calculator.";
pub const description =
    \\This calculator is for computing basic quantities related to black holes. It is
    \\currently a minimal ISCO and event horizon radius calculator, but will be
    \\fleshed out with additional features in the future.
    \\
;

pub const Args = clippy.Arguments(&[_]clippy.ArgumentDescriptor{
    .{
        .arg = "--mass mass",
        .argtype = f64,
        .default = "1",
        .help = "The black hole mass.",
    },
    .{
        .arg = "--spin spin",
        .argtype = f64,
        .default = "0.998",
        .help = "The black hole spin.",
    },
    .{
        .arg = "--radius rg",
        .argtype = f64,
        .help = "The radius in rg at which to evaluate e.g. timescales. It is ignored for expressions that are not expecting this argument.",
    },
    .{
        .arg = "--theta ang",
        .argtype = f64,
        .default = "90",
        .help = "The theta angle (degrees) in Boyer-Lindquist coordinates at which to evaluate the calculation.",
    },
    .{
        .arg = "--alpha alpha",
        .argtype = f64,
        .default = "0.1",
        .help = "The alpha-disc anomalous stress parameter (Shakura & Sunyaev, 1973). The default value is chosen from King et al. (2007).",
    },
    .{
        .arg = "expression",
        .default = "all",
        .help = "What to calculate. This can be " ++ utils.makeList(Expression),
    },
    utils.VelocityProfile.ArgDescriptor,
});

const Expression = enum {
    isco,
    horizon,
    photon_orbit,
    gamma,
    timescales,
    dynamical,
    viscous,
    tetrad,
    metric,
    all,
};

const Dual = kerrz.DualNumber(f64, 1);

/// GM / c
fn factor_GM_over_c(solar_mass: f64) f64 {
    return Constants.big_G * solar_mass * Constants.mass_sun / Constants.speed_of_light;
}

/// GM / c^2
fn gravitational_radius(solar_mass: f64) f64 {
    return factor_GM_over_c(solar_mass) / (Constants.speed_of_light);
}

/// GM / c^3
fn gravitational_time(solar_mass: f64) f64 {
    return gravitational_radius(solar_mass) / Constants.speed_of_light;
}

/// Execute the completion helper
pub fn run(out: *std.Io.Writer, allocator: std.mem.Allocator, itt: *clippy.ArgumentIterator) !void {
    var arg_parser = Args.init(itt, .{});
    const args_with_overflow = try arg_parser.parseWithOverflow(allocator);
    defer args_with_overflow.deinit(allocator);
    const args = args_with_overflow.args;

    const expr: Expression = std.meta.stringToEnum(Expression, args.expression) orelse
        return error.InvalidExpression;

    const metric = kerrz.KerrMetric(Dual).init(.one, Dual.promote(args.spin).diff(0));
    const calc: Calculator = .{
        .metric = metric,
        .solar_mass = args.mass,
        .radius = args.radius orelse metric.isco.x,
        .theta = std.math.degreesToRadians(args.theta),
        .alpha = args.alpha,
        .diff = false,
        .vel_prof = try utils.VelocityProfile.fromArgs(args),
    };

    try out.print(
        "M = {d} M_sol, a = {d}, r = {d:.3} rg, θ = {d:.2}°:\nKey: 1 rg = {e:.5} m, 1 tg = {e:.5} s\n\n",
        .{
            calc.solar_mass,
            calc.metric.a.x,
            calc.radius,
            std.math.radiansToDegrees(calc.theta),
            gravitational_radius(calc.solar_mass),
            gravitational_time(calc.solar_mass),
        },
    );
    try calc.calculate(out, expr);
    try out.writeAll("\n");
    try out.flush();
}

pub const Calculator = struct {
    metric: kerrz.KerrMetric(Dual),
    solar_mass: f64,
    radius: f64,
    theta: f64,
    diff: bool,
    alpha: f64,
    vel_prof: kerrz.orbits.VelocityProfiles,

    fn printValue(out: *std.Io.Writer, name: []const u8, value: f64, units: []const u8) !void {
        if (value > 1e4 or value < 1e-3) {
            try out.print("{s:>20} : {e: <12.6} {s}\n", .{ name, value, units });
        } else {
            try out.print("{s:>20} : {d: <12.6} {s}\n", .{ name, value, units });
        }
    }

    fn dynamicalTime(self: Calculator) f64 {
        const keplerian = kerrz.orbits.Keplerian(Dual).init(
            self.metric,
            .promote(self.radius),
            .promote(self.theta),
        );
        // Need the factor 2 pi to convert from angular frequency.
        return std.math.pi * 2.0 / keplerian.angular_frequency().x;
    }

    fn scaleHeightRatio(self: Calculator) f64 {
        _ = self;
        return 0.01;
    }

    pub fn calculate(self: Calculator, out: *std.Io.Writer, expr: Expression) !void {
        switch (expr) {
            .isco => {
                try printValue(out, "r_isco", self.metric.isco.x, "rg");
                try printValue(out, "", self.metric.isco.x * gravitational_radius(self.solar_mass), "m");
            },
            .photon_orbit => {
                const photon_orbit = self.metric.photonOrbit();
                try printValue(out, "r_photon", photon_orbit.x, "rg");
            },
            .horizon => {
                try printValue(out, "r_horizon", self.metric.horizon_radius.x, "rg");
                try printValue(out, "", self.metric.horizon_radius.x * gravitational_radius(self.solar_mass), "m");
            },
            .dynamical => {
                const t_dyn = self.dynamicalTime();
                try printValue(out, "Dynamical / Orbital", 1 / t_dyn, "tg⁻¹");
                const t_dyn_physical = gravitational_time(self.solar_mass) * t_dyn;
                try printValue(out, "", 1 / t_dyn_physical, "Hz");
            },
            .gamma => {
                const ts = self.metric.tangentSpaceAlt(
                    .promote(self.radius),
                    .promote(std.math.pi / 2.0),
                );
                const gamma = kerrz.orbits.lorentzFactorKeplerian(Dual, self.metric, ts);
                try printValue(out, "Gamma", gamma.x, "");
            },
            .viscous => {
                const t_dyn = self.dynamicalTime();
                const height_ratio = self.scaleHeightRatio();
                const t_visc = t_dyn / (height_ratio * height_ratio * self.alpha);
                try printValue(out, "Viscous", t_visc, "tg⁻¹");
                const t_visc_physical = gravitational_time(self.solar_mass) * t_visc;
                try printValue(out, "", t_visc_physical, "Hz");
            },
            .timescales => {
                try self.calculate(out, .dynamical);
                try self.calculate(out, .viscous);
            },
            .metric => {
                const x: kerrz.FourVector(Dual) = .{
                    .t = .zero,
                    .r = .promote(self.radius),
                    .th = .promote(self.theta),
                    .ph = .zero,
                };
                const ts = self.metric.tangentSpace(x);

                try out.writeAll("  Regular g_μν:\n\n");

                for (0..4) |col| {
                    try out.print("{s: >11}", .{COORDINATE_NAMES[col]});
                }
                try out.writeAll("\n");

                const matrix = ts.metric_components.toMatrix();

                for (0..4) |row| {
                    for (0..4) |col| {
                        const index = row + col * 4;
                        const element = matrix.v[index];
                        try out.print("{d: >11.5}", .{element.x});
                    }
                    try out.writeAll("\n");
                }

                const inverse_matrix = ts.metric_components.inverse().toMatrix();

                try out.writeAll("\n  Inverse g^μν:\n\n");

                for (0..4) |row| {
                    for (0..4) |col| {
                        const index = row + col * 4;
                        const element = inverse_matrix.v[index];
                        try out.print("{d: >11.5}", .{element.x});
                    }
                    try out.writeAll("\n");
                }
            },
            .tetrad => {
                const x: kerrz.FourVector(Dual) = .{
                    .t = .zero,
                    .r = .promote(self.radius),
                    .th = .promote(self.theta),
                    .ph = .zero,
                };
                const ts = self.metric.tangentSpace(x);
                const v = self.vel_prof.fourVector(Dual, self.metric, ts);
                const frame = ts.localFrame(v);

                for (0..4) |col| {
                    try out.print("{s: >12}", .{LOCAL_TETRAD_NAMES[col]});
                }
                try out.writeAll("\n");

                for (0..4) |row| {
                    for (0..4) |col| {
                        const index = row + col * 4;
                        const element = frame.m.v[index];
                        try out.print("{d: >10.5}", .{element.x});
                    }
                    try out.writeAll("\n");
                }
            },
            .all => {
                inline for (@typeInfo(Expression).@"enum".fields) |field| {
                    const _f = std.meta.stringToEnum(Expression, field.name).?;
                    switch (_f) {
                        .all, .timescales, .tetrad, .metric => {},
                        else => try self.calculate(out, _f),
                    }
                }
            },
        }
    }
};
