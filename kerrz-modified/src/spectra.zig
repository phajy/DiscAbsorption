const std = @import("std");
const ad = @import("zad");
const TEST_TOLERANCE = @import("options").test_numerical_tolerance;

const interpolations = @import("interpolations.zig");

/// Evaluate a blackbody spectrum at a particular energy and given a particular
/// temperature kT.
///
/// The units are those of specific intensity `I_E`.
pub fn blackbody(comptime T: type, energy: T, kT: T) T {
    const A = T.Algebra;

    const denom = A.sub(A.exp(A.div(energy, kT)), .one);
    const numerator = A.mult(.promote(2), A.powi(energy, 3));

    return A.div(numerator, denom);
}

test "blackbody" {
    const Dual = ad.DualNumber(f64, 0);
    const v = blackbody(Dual, .promote(5.0), .promote(3.0));
    try std.testing.expectApproxEqAbs(58.21412951524655, v.x, TEST_TOLERANCE);
}

/// A structure that represents (a table for) the reflection spectrum.
pub const DiscSpectrum = struct {
    energy_grid: []const f64,
    radii: []const f64,
    fluxes: []const []const f64,
    /// A temporary buffer that is used for interpolations.
    flux_buffer: []f64,

    pub fn deinit(self: DiscSpectrum, allocator: std.mem.Allocator) void {
        for (self.fluxes) |f_grid| {
            allocator.free(f_grid);
        }
        allocator.free(self.fluxes);
        allocator.free(self.flux_buffer);
        allocator.free(self.energy_grid);
        allocator.free(self.radii);
    }

    fn sortedPredicate(this: f64, other: f64) bool {
        return this > other;
    }

    /// Interpolate the spectral profile at a particular radius. Returns null
    /// if the radius is out of bounds.
    pub fn atRadius(self: DiscSpectrum, r: f64) ?[]const f64 {
        if (r < self.radii[0] or r > self.radii[self.radii.len - 1]) {
            return null;
        }

        const index = std.sort.partitionPoint(
            f64,
            self.radii[0 .. self.radii.len - 2],
            r,
            sortedPredicate,
        );

        const r1 = self.radii[index];
        const r2 = self.radii[index + 1];

        const flux_1 = self.fluxes[index];
        const flux_2 = self.fluxes[index + 1];

        const w = interpolations.lerpWeight(f64, r, r1, r2);
        for (self.flux_buffer, flux_1, flux_2) |*buf, f1, f2| {
            buf.* = interpolations.lerpValue(f64, w, f1, f2);
        }

        return self.flux_buffer;
    }

    /// Read a disc spectrum from a file. The format of the file should be as
    /// defined below. Most values are comma seperated, using the new-line to
    /// denote the end of a particular radius's flux grid.
    ///
    ///     g0, g1, g2, ...              : The energy grid of of the disc.
    ///     r0: f0, f1, f2, f3, ...      : The flux grid at radius r0.
    ///     r1: f0, f1, f2, f3, ...      : The flux grid at radius r1.
    ///     ...
    ///
    /// Caller owns the memory and must call deinit.
    pub fn parseFromFile(allocator: std.mem.Allocator, path: []const u8) !DiscSpectrum {
        const dir = std.fs.cwd();
        const stat = try dir.statFile(path);
        const contents = try std.fs.cwd().readFileAlloc(allocator, path, stat.size);
        defer allocator.free(contents);

        // The radius of each spectrum.
        var radii = std.ArrayList(f64).empty;
        defer radii.deinit(allocator);

        // The spectrum of each radius.
        var spectra = std.ArrayList([]const f64).empty;
        defer spectra.deinit(allocator);
        errdefer for (spectra.items) |item| {
            allocator.free(item);
        };

        var line_itt = std.mem.tokenizeAny(u8, contents, "\n");

        // Parse the energy grid.
        const g_line = line_itt.next() orelse return error.InvalidDiscSpectrum;
        const e_grid = try readCSV(allocator, g_line);
        errdefer allocator.free(e_grid);

        // Obtain the radius and the flux for that radius.
        while (line_itt.next()) |line| {
            if (line.len == 0) continue;
            const split = std.mem.indexOf(u8, line, ":") orelse
                return error.InvalidDiscSpectrum;
            const radius_data = line[0..split];
            const spectrum_data = line[split + 1 .. line.len];

            try radii.append(allocator, try std.fmt.parseFloat(f64, radius_data));
            const spectrum = try readCSV(allocator, spectrum_data);
            try spectra.append(allocator, spectrum);
        }

        const buffer = try allocator.alloc(f64, spectra.items[0].len);
        errdefer allocator.free(buffer);

        const _radii = try radii.toOwnedSlice(allocator);
        errdefer allocator.free(_radii);

        const _fluxes = try spectra.toOwnedSlice(allocator);
        errdefer allocator.free(_fluxes);

        return .{
            .energy_grid = e_grid,
            .flux_buffer = buffer,
            .radii = _radii,
            .fluxes = _fluxes,
        };
    }
};

/// Parse a given line of `data` as a row in the CSV file. Caller owns the
/// memory.
fn readCSV(allocator: std.mem.Allocator, data: []const u8) ![]f64 {
    var list = std.ArrayList(f64).empty;
    defer list.deinit(allocator);
    var itt = std.mem.tokenizeScalar(u8, data, ',');
    while (itt.next()) |untrimmed_token| {
        const token = std.mem.trim(u8, untrimmed_token, " ");
        const value = try std.fmt.parseFloat(f64, token);
        try list.append(allocator, value);
    }
    return list.toOwnedSlice(allocator);
}

test "spectral interpolations" {
    var dir = std.testing.tmpDir(.{});
    defer dir.cleanup();
    {
        const file = try dir.dir.createFile("temp.ds", .{});
        defer file.close();
        try file.writeAll(
            \\0,1,2,3,4,5,6,7,8,9
            \\1: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
            \\2: 0, 2, 4, 6, 8, 10, 12, 14, 16, 18
            \\3: 0, 3, 6, 9, 12, 15, 18, 21, 24, 27
            \\4: 0, 4, 8, 12, 16, 20, 24, 28, 32, 36
        );
    }
    const tmppath = try dir.dir.realpathAlloc(
        std.testing.allocator,
        "temp.ds",
    );
    defer std.testing.allocator.free(tmppath);

    const ds = try DiscSpectrum.parseFromFile(std.testing.allocator, tmppath);
    defer ds.deinit(std.testing.allocator);

    try std.testing.expectEqualSlices(f64, &.{ 1, 2, 3, 4 }, ds.radii);
    try std.testing.expectEqualSlices(f64, &.{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 }, ds.fluxes[0]);
}

/// Convolves a and b on domains `g` and `x` respectively, by rebinning if
/// necessary. Assumes the `output` is on domain `x`.
///
/// Does not zero the output slice. Caller must do so if needed.
///
/// Approximately follows the method of constructing a Toeplitz matrix on an
/// (irregular) grid for the convolution, but without storing the matrix in
/// memory.
pub fn convolve(
    comptime T: type,
    output: []T,
    g: []const T,
    a: []const T,
    x: []const T,
    b: []const T,
) void {
    std.debug.assert(output.len == b.len);

    // Find window of non-zero.
    const window_start = blk: {
        for (0..a.len) |i| if (a[i] > 0)
            break :blk if (i > 0) i - 1 else i;
        break :blk 0;
    };
    const window_end = blk: {
        for (0..a.len) |i| {
            const j = a.len - i - 1;
            if (a[j] > 0)
                break :blk if (j < a.len - 1) j + 1 else j;
        }
        break :blk 0;
    };

    // Find the extremal.
    const g_min = g[window_start];
    const g_max = g[window_end];

    for (0..output.len) |i| {
        const avg = 2 / (x[i + 1] + x[i]);
        for (0..output.len) |j| {
            // Output bin extremes shifted.
            const low = x[j] * avg;
            const high = x[j + 1] * avg;

            // Skip if outside of the window.
            if (high < g_min or low > g_max) continue;

            // Integrate the window.
            const weight = convolution_weight(T, window_start, window_end, g, a, low, high);

            // Accumulate the dot product.
            output[j] += weight * b[i];
        }
    }
}

fn convolution_weight(
    comptime T: type,
    start: usize,
    end: usize,
    g: []const T,
    a: []const T,
    low: T,
    high: T,
) T {
    var weight: T = 0;
    for (start..end) |w| {
        const bin_low = g[w];
        const bin_high = g[w + 1];
        const bin_width = bin_high - bin_low;
        // Check different cases to calculate overlap amount.
        // 1. Check if bin is not in output bin:
        //    bin_low |---| bin_high     low |---| high
        if (bin_high < low) {
            continue;
        } else if (bin_low > high) {
            // Break the loop early.
            break;
        }
        var overlap: T = 1;
        // 2. Check if bin is bigger than output bin:
        //    bin_low |...| low |.....| high |...| bin_high
        if (bin_low < low and bin_high > high) {
            overlap = (high - low) / bin_width;
        }
        // 3. Check if bin is striding:
        //    bin_low |---| low |.....| bin_high |---| high
        else if (bin_low < low and bin_high < high) {
            overlap = (bin_high - low) / bin_width;
        }
        // or
        //    low |---| bin_low |.....| high |---| bin_high
        else if (bin_low > low and bin_high > high) {
            overlap = (high - bin_low) / bin_width;
        }
        // 4. Bin must be contained, and overlap is 1:
        //    low |---| bin_low |.....| bin_high |---| high
        weight += a[w] * bin_width * overlap;
    }
    return weight;
}
