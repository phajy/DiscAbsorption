const std = @import("std");
const ad = @import("zad");

/// Interpolate value given a weight.
pub inline fn lerpValue(comptime T: type, w: T, y1: T, y2: T) T {
    if (@typeInfo(T) == .float) {
        return w * y2 + (1 - w) * y1;
    } else {
        const A = T.Algebra;
        return A.add(A.mult(w, y2), A.mult(A.sub(.one, w), y1));
    }
}

/// Interpolation weight for use with `lerpValue`.
pub inline fn lerpWeight(comptime T: type, x: T, x1: T, x2: T) T {
    if (@typeInfo(T) == .float) {
        return (x - x1) / (x2 - x1);
    } else {
        const A = T.Algebra;
        return A.div(A.sub(x, x1), A.sub(x2, x1));
    }
}

/// Linear interpolation.
pub inline fn lerp(comptime T: type, x: T, x1: T, x2: T, y1: T, y2: T) T {
    if (@typeInfo(T) == .float) {
        const Dual0 = ad.DualNumber(T, 0);
        return lerpValue(
            Dual0,
            lerpWeight(Dual0, .promote(x), .promote(x1), .promote(x2)),
            .promote(y1),
            .promote(y2),
        ).x;
    } else {
        return lerpValue(T, lerpWeight(T, x, x1, x2), y1, y2);
    }
}

/// Intepolate the missing values (i.e. those that are zero) with values either
/// side with a linear interpolation.
pub fn interpolateZeroes(comptime T: type, values: []T) void {
    var left: usize = 0;

    // Find the first non-zero value:
    for (0.., values) |i, v| {
        if (v != 0) {
            left = i;
            break;
        }
    } else {
        // Everything is zero: can't interpolate anything.
        return;
    }

    // Save the start for interpolating later.
    const start = left;

    var right: usize = 1;

    outer: while (left < values.len) {
        // Find the next zero value:
        for (left..values.len) |i| {
            if (values[i] == 0) {
                left = i - 1;
                break;
            }
        } else {
            left = values.len;
            // No zeroes left.
            break :outer;
        }

        // Find the next non-zero value:
        for (left + 1..values.len) |i| {
            if (values[i] != 0) {
                right = i;
                break;
            }
        } else {
            // No non-zeroes left: need to extrapolate rest.
            break :outer;
        }

        // Linearly interpolate everything between left and right.
        const lt: T = @floatFromInt(left);
        const rt: T = @floatFromInt(right);
        for (left + 1..right) |i| {
            const w = lerpWeight(T, @floatFromInt(i), lt, rt);
            values[i] = values[left] * (1 - w) + values[right] * w;
        }

        // Advance to the next element.
        left = right + 1;
    }

    // TODO: optional periodicity

    // Extrapolate any leading zeroes:
    if (left < values.len and left > 1 and values[left] != 0 and values[left - 1] != 0) {
        const lt: T = @floatFromInt(left - 1);
        const rt: T = @floatFromInt(left);
        for (left + 1..values.len) |i| {
            const w = lerpWeight(T, @floatFromInt(i), lt, rt);
            values[i] = values[left - 1] * (1 - w) + values[left] * w;
        }
    }

    // Extrapolate leading zeroes using the first two non-zero values:
    if (start > 0 and start + 1 < values.len and values[start + 1] != 0) {
        const lt: T = @floatFromInt(start);
        const rt: T = @floatFromInt(start + 1);
        for (0..start) |i| {
            const w = lerpWeight(T, @floatFromInt(i), lt, rt);
            values[i] = values[start] * (1 - w) + values[start + 1] * w;
        }
    }
}

test "interpolate zeroes" {
    {
        var values = [_]f64{ 1, 2, 3, 4, 0, 6, 0, 8, 9, 0, 11, 12, 13, 14, 15 };
        interpolateZeroes(f64, &values);
        try std.testing.expectEqualSlices(
            f64,
            &values,
            &.{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 },
        );
    }

    {
        var values = [_]f64{ 0, 2, 3, 4, 0, 6, 0, 8, 9, 0, 11, 12, 13, 0, 0 };
        interpolateZeroes(f64, &values);
        try std.testing.expectEqualSlices(
            f64,
            &values,
            &.{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 },
        );
    }

    {
        var values = [_]f64{ 0, 1, 0 };
        interpolateZeroes(f64, &values);
        try std.testing.expectEqualSlices(
            f64,
            &values,
            &.{ 0, 1, 0 },
        );
    }
}

pub fn Indices(comptime T: type) type {
    return struct {
        left: usize,
        right: usize,
        weight: T,
        // Whether 2π (i.e. the modulus) was subtracted from the left point.
        left_mod: bool,
    };
}

/// Used to find the first element of sorted array that has been
/// cyclic-permuted some number of times.
///
/// The index returned is the smallest element.
fn indexOfFirstOrdered(comptime T: type, axis: []const T, x: anytype, comptime keyFunc: fn (T) @TypeOf(x)) usize {
    var low: usize = 0;
    var high: usize = axis.len;
    const start = keyFunc(axis[0]);
    while (low < high) {
        const mid = low + (high - low) / 2;
        if (start < keyFunc(axis[mid])) {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    return low;
}

/// For periodic interpolations e.g. over azimuth.
pub fn periodicInterpolate(
    comptime T: type,
    axis: []const T,
    x: anytype,
    comptime keyFunc: fn (T) @TypeOf(x),
) Indices(@TypeOf(x)) {
    const _FunctionTable = struct {
        fn keyPartition(ctx: @TypeOf(x), t: T) bool {
            return keyFunc(t) < ctx;
        }
    };

    // Although the angles are sorted, the first angle is not
    // necessarily the smallest angle. So we first need to know where
    // the zero point is and then pivot based on that.
    var slice_start: usize = 0;
    var slice_end: usize = axis.len;

    if (keyFunc(axis[axis.len - 1]) < keyFunc(axis[0])) {
        const pivot = indexOfFirstOrdered(T, axis, x, keyFunc);
        if (x < keyFunc(axis[0])) {
            // The index will be after the pivot:
            slice_start = pivot;
        } else {
            // The index will be before the pivot:
            slice_end = pivot;
        }
    }

    const i = slice_start + std.sort.partitionPoint(
        T,
        axis[slice_start..slice_end],
        x,
        _FunctionTable.keyPartition,
    );

    var left_i = i -| 1;
    var right_i = i;

    // Edge cases periodicity:
    if (right_i >= axis.len) {
        right_i = 0;
    }

    if (right_i == 0) {
        left_i = axis.len - 1;
    }

    var phi1 = keyFunc(axis[left_i]);
    const phi2 = keyFunc(axis[right_i]);

    // Handle modulus:
    var left_mod: bool = false;
    if (phi1 > phi2) {
        phi1 -= std.math.pi * 2.0;
        left_mod = true;
    }

    const w = lerpWeight(@TypeOf(x), x, phi1, phi2);
    return .{
        .left = left_i,
        .right = right_i,
        .weight = w,
        .left_mod = left_mod,
    };
}
