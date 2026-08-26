const std = @import("std");
const ad = @import("zad");
const options = @import("options");

const TEST_TOLERANCE = options.test_numerical_tolerance;

const threads = @import("threads.zig");
const utils = @import("utils.zig");

pub const geodesic = @import("geodesic.zig");
pub const accretion_discs = @import("accretion-discs.zig");
pub const geometry = @import("geometry.zig");
pub const transfer_tables = @import("transfer-tables.zig");
pub const integration_state = @import("./integration-state.zig");
pub const potentials = @import("potentials.zig");
pub const antiderivatives = @import("antiderivatives.zig");
pub const elliptic_integrals = @import("elliptic-integrals.zig");
pub const redshift = @import("redshift.zig");
pub const orbits = @import("orbits.zig");
pub const emissivity = @import("emissivity.zig");
pub const polarisation = @import("polarisation.zig");
pub const interpolations = @import("interpolations.zig");
pub const spectra = @import("spectra.zig");
pub const iterators = @import("iterators.zig");
pub const tools = @import("tools.zig");
pub const matrix = @import("matrix.zig");
pub const continuum = @import("continuum.zig");

pub const ThreadMap = threads.ThreadMap;
pub const version = options.version;
pub const KerrMetric = geometry.KerrMetric;
pub const FourVector = geometry.FourVector;
pub const DualNumber = ad.DualNumber;
pub const NullGeodesic = geodesic.NullGeodesic;
pub const TraceResult = geodesic.TraceResult;
pub const AccretionDisc = accretion_discs.AccretionDisc;
pub const EmissivityCache = emissivity.EmissivityCache;
pub const EmissivityProfile = emissivity.EmissivityProfile;
pub const EmissivityTrace = emissivity.EmissivityTrace;
pub const Status = geodesic.Status;
pub const Mapper = @import("Mapper.zig");
pub const CoronalModelOptions = emissivity.CoronalModelOptions;
pub const CoronalModel = emissivity.CoronalModel;
pub const Matrix = matrix.Matrix;
pub const CunninghamTransferFunctionTable = transfer_tables.CunninghamTransferFunctionTable;

pub const FitsFile = @import("zfits").FitsFile;
pub const addKerrzFITSInfo = utils.addKerrzFITSInfo;

test "all" {
    _ = elliptic_integrals;
    _ = antiderivatives;
    _ = @import("complex.zig");
    _ = redshift;
    _ = threads;
    _ = geometry;
    _ = geodesic;
    _ = accretion_discs;
    _ = transfer_tables;
    _ = orbits;
    _ = emissivity;
    _ = polarisation;
    _ = integration_state;
    _ = interpolations;
    _ = spectra;
    _ = iterators;
    _ = matrix;
    _ = continuum;
    _ = utils;
}

fn testImpactParameters(
    comptime T: type,
    metric: KerrMetric(T),
    r: T.T,
    incl: T.T,
    alpha: T.T,
    beta: T.T,
) FourVector(T) {
    const x: FourVector(T) = .{
        .t = .zero,
        .r = .promote(r),
        .th = .promote(std.math.degreesToRadians(incl)),
        .ph = .zero,
    };
    const geod = NullGeodesic(T).fromImpactParameters(
        metric,
        x,
        .promote(alpha),
        .promote(beta),
    );
    const result = geod.traceToAngle(
        metric,
        .promote(std.math.pi / 2.0),
        .{},
    );

    const total = result.totalAntiderivatives(metric, geod);
    const t = total.coordinateTime(metric, geod);
    const phi = total.coordinateAzimuth(metric, geod);

    return .{
        .t = t,
        .r = result.r,
        .th = result.theta,
        .ph = phi,
    };
}

test "problem cases" {
    // these test cases represent various problematic geodesics that have been
    // encountered
    const Dual = DualNumber(f64, 0);
    const Vec = FourVector(Dual);
    const M = KerrMetric(Dual);
    const G = NullGeodesic(Dual);

    const metric: M = .init(.one, .promote(0.998));
    {
        // This triggers a case II geodesic but r < r4, which should not be
        // allowed for case II. It should therefore be a case I geodesic.
        const x: Vec = .{
            .t = .zero,
            .r = .promote(2.8),
            .th = .promote(std.math.degreesToRadians(45)),
            .ph = .zero,
        };
        const geod = G.fromStationarySkyAngles(
            metric,
            x,
            .promote(1.5707963267948966 + 1e-5),
            .promote(4.71238898038469 + 1e-5),
        );
        const result = geod.traceToAngle(
            metric,
            .promote(std.math.pi / 2.0),
            .{},
        );
        try std.testing.expectEqual(potentials.RadialCase.case_I, result.state.radial_case);
    }

    {
        // This is from the reltrans grtrace test, which produces subtly
        // different values on MacOS and linux.
        const res1 = testImpactParameters(Dual, metric, 1.8e7, 30, 2.2, 2.2);
        try std.testing.expectApproxEqAbs(18000041.18860833, res1.t.x, TEST_TOLERANCE);

        const res2 = testImpactParameters(Dual, metric, 1.8e7, 30, 1.2, -2.2);
        try std.testing.expectApproxEqAbs(18000033.528562296, res2.t.x, TEST_TOLERANCE);

        const res3 = testImpactParameters(Dual, metric, 1.8e7, 1, 2.0, -2.0);
        try std.testing.expectApproxEqAbs(18000036.02186861, res3.t.x, TEST_TOLERANCE);

        const res4 = testImpactParameters(Dual, metric, 1.8e7, 30, -3.2997356434211702, -0.40670451835608296);
        try std.testing.expectApproxEqAbs(18000033.497825205, res4.t.x, TEST_TOLERANCE);

        const res5 = testImpactParameters(Dual, metric, 1.8e7, 30, -3.283356241887061, -0.49626494359898454);
        try std.testing.expectApproxEqAbs(18000033.41745938, res5.t.x, TEST_TOLERANCE);

        const res6 = testImpactParameters(Dual, metric, 1.8e7, 30, -1.7770090146700821, -0.21902287744375143);
        try std.testing.expectApproxEqAbs(18000112.906015214, res6.t.x, TEST_TOLERANCE);

        const res7 = testImpactParameters(Dual, metric, 1.8e7, 30, 3.0, 4.0);
        try std.testing.expectApproxEqAbs(34.66124926134944, res7.t.x - 1.8e7, TEST_TOLERANCE);
    }

    {
        const x: FourVector(Dual) = .{
            .t = .zero,
            .r = .promote(1.8e7),
            .th = .promote(0.5235987755982987),
            .ph = .zero,
        };

        const geod = NullGeodesic(Dual).fromImpactParameters(
            metric,
            x,
            .promote(3.0),
            .promote(4.0),
        );

        const result = geod.traceToAngle(
            metric,
            .promote(std.math.pi / 2.0),
            .{},
        );

        const total = result.totalAntiderivatives(metric, geod);
        const t = total.coordinateTime(metric, geod);
        try std.testing.expectApproxEqAbs(34.66124925017357, t.x - 1.8e7, TEST_TOLERANCE);
    }
}

fn testEqualPositions(
    comptime T: type,
    v: FourVector(T),
    t: T.T,
    r: T.T,
    theta: T.T,
    phi: T.T,
) !void {
    try std.testing.expectApproxEqAbs(
        t,
        v.t.x,
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        r,
        v.r.x,
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        theta,
        v.th.x,
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        phi,
        v.ph.x,
        TEST_TOLERANCE,
    );
}

test "path construction" {
    // These tests are motivated by the need for the turning points to be
    // correctly considered when tracing or reconstructing an individual path
    // of a geodesic.
    const Dual = DualNumber(f64, 0);
    const V = FourVector(Dual);
    const M = KerrMetric(Dual);
    const G = NullGeodesic(Dual);
    {
        // First for the highly spinning case. The points of comparison are
        // computed with Gradus.jl and check that the sign of the coordinates
        // are correct a small degree either side of the turning point.
        const metric: M = .init(.one, .promote(0.998));
        const x: V = .{
            .t = .zero,
            .r = .promote(1e4),
            .ph = .zero,
            .th = .promote(std.math.degreesToRadians(80)),
        };
        const geod = G.fromImpactParameters(metric, x, .promote(-2), .promote(2));
        const builder = geod.traceBuilder(metric, .{});

        try std.testing.expectApproxEqAbs(
            0.47641735137028607,
            builder.theta_1_time.x,
            TEST_TOLERANCE,
        );

        const p0 = builder.atMinoTime(.promote(0.4)).position(metric, geod);
        try testEqualPositions(
            Dual,
            p0,
            10016.924561675976,
            2.823408345399896,
            0.7663313854975166,
            -1.3418461148251275,
        );

        const p1 = builder.atMinoTime(.promote(1.5)).position(metric, geod);
        try testEqualPositions(
            Dual,
            p1,
            10036.875550704603,
            1.313408544377003,
            2.3557516053001675,
            -8.562960764735635,
        );

        const p2 = builder.atMinoTime(.promote(1.7)).position(metric, geod);
        try testEqualPositions(
            Dual,
            p2,
            10041.756280593727,
            1.2582836402485993,
            2.359853386899176,
            -4.422927778074036,
        );
    }

    {
        // This path has two turning points, the second of which was once not
        // handled correctly.
        const metric: M = .init(.one, .promote(0.998));
        const x: V = .{
            .t = .zero,
            .r = .promote(5.0),
            .ph = .zero,
            .th = .promote(std.math.degreesToRadians(0.1)),
        };
        const geod = G.fromStationarySkyAngles(
            metric,
            x,
            .promote(std.math.degreesToRadians(132)),
            .zero,
        );
        const builder = geod.traceBuilder(metric, .{});

        try std.testing.expectApproxEqAbs(
            -0.0003607627901596744,
            builder.theta_0_time.x,
            TEST_TOLERANCE,
        );
        try std.testing.expectApproxEqAbs(
            0.6560911887770564,
            builder.theta_1_time.x,
            TEST_TOLERANCE,
        );
    }
}

test "vortical motion" {
    const Dual = DualNumber(f64, 0);
    const V = FourVector(Dual);
    const M = KerrMetric(Dual);
    const G = NullGeodesic(Dual);

    const metric: M = .init(.one, .promote(0.998));
    const x: V = .{
        .t = .zero,
        .r = .promote(5),
        .th = .promote(std.math.degreesToRadians(0.1)),
        .ph = .zero,
    };

    // Find a pair of geodesics which are very similar but where one is normal
    // angular motion and the other is vortical. Use a small pertubation of the
    // mino time to ensure the results are approximately consistent.
    const geod_normal = G.fromStationarySkyAngles(
        metric,
        x,
        .promote(std.math.degreesToRadians(171.267)),
        .zero,
    );
    const geod_vortical = G.fromStationarySkyAngles(
        metric,
        x,
        .promote(std.math.degreesToRadians(171.268)),
        .zero,
    );

    const builder_normal = geod_normal.traceBuilder(metric, .{});
    const builder_vortical = geod_vortical.traceBuilder(metric, .{});

    try std.testing.expect(builder_normal.state.angular_case == .normal);
    try std.testing.expect(builder_vortical.state.angular_case == .vortical);

    // Apply a small Mino time to the trace
    const result_delta_normal = builder_normal.atMinoTime(.promote(1e-3));
    const result_delta_vortical = builder_vortical.atMinoTime(.promote(1e-3));

    // Test direct quantities
    try std.testing.expectApproxEqAbs(
        result_delta_normal.r.x,
        result_delta_vortical.r.x,
        TEST_TOLERANCE,
    );
    try std.testing.expectApproxEqAbs(
        result_delta_normal.theta.x,
        result_delta_vortical.theta.x,
        TEST_TOLERANCE,
    );

    // Test derived quantities
    const total_normal = result_delta_normal.totalAntiderivatives(
        metric,
        geod_normal,
    );
    const total_vortical = result_delta_vortical.totalAntiderivatives(
        metric,
        geod_vortical,
    );

    try std.testing.expectApproxEqAbs(
        total_normal.coordinateAzimuth(metric, geod_normal).x,
        total_vortical.coordinateAzimuth(metric, geod_vortical).x,
        TEST_TOLERANCE,
    );

    try std.testing.expectApproxEqAbs(
        total_normal.coordinateTime(metric, geod_normal).x,
        0.9757708658549847 * total_vortical.coordinateTime(metric, geod_vortical).x,
        TEST_TOLERANCE,
    );
}

test "mino time to coordinates" {
    const Dual = DualNumber(f64, 0);
    const V = FourVector(Dual);
    const M = KerrMetric(Dual);
    const G = NullGeodesic(Dual);

    const metric: M = .init(.one, .promote(0.998));
    const x: V = .{
        .t = .zero,
        .r = .promote(1e4),
        .th = .promote(std.math.degreesToRadians(80)),
        .ph = .zero,
    };

    const geod = G.fromImpactParameters(
        metric,
        x,
        .promote(10.0),
        .promote(1e-3),
    );

    const path_builder = geod.traceBuilder(metric, .{});

    {
        const mino = path_builder.state.minoTimeToRadius(
            .promote(50.0),
            @floatFromInt(geod.radial_sign),
        );
        const result = path_builder.atMinoTime(mino);
        try std.testing.expectApproxEqAbs(
            50.0,
            result.r.x,
            TEST_TOLERANCE,
        );
    }

    {
        const mino = path_builder.state.minoTimeToRadius(
            .promote(20.0),
            @floatFromInt(geod.radial_sign),
        );
        const result = path_builder.atMinoTime(mino);
        try std.testing.expectApproxEqAbs(
            20.0,
            result.r.x,
            TEST_TOLERANCE,
        );
    }
}
