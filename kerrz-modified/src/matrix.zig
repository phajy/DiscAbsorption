const std = @import("std");
const zad = @import("zad");
const options = @import("options");

const TEST_TOLERANCE = options.test_numerical_tolerance;

/// A matrix implementation designed to live on the stack. The matrix is
/// row-major.
pub fn SMatrix(comptime T: type, comptime nRows: usize, comptime nCols: usize) type {
    return struct {
        const Self = @This();
        v: [nRows * nCols]T,

        fn getRow(self: Self, row_index: usize) [nCols]T {
            var row: [nCols]T = undefined;
            for (0..nCols) |i| {
                const index = row_index * nCols + i;
                row[i] = self.v[index];
            }
            return row;
        }

        fn getCol(self: Self, col_index: usize) [nRows]T {
            var col: [nRows]T = undefined;
            for (0..nRows) |i| {
                const index = i * nCols + col_index;
                col[i] = self.v[index];
            }
            return col;
        }

        /// Multiply to the right (M * v) with a vector v.
        pub fn multR(self: Self, vec: [nCols]T) [nRows]T {
            const A = T.Algebra;
            var out: [nRows]T = undefined;
            for (0..nRows) |i| {
                const row = self.getRow(i);
                out[i] = .zero;
                for (row, vec) |r, v| {
                    out[i] = A.add(A.mult(r, v), out[i]);
                }
            }
            return out;
        }

        /// Multiply to the right after transposing (M.T * v) with a vector v.
        pub fn multR_T(self: Self, vec: [nRows]T) [nCols]T {
            const A = T.Algebra;
            var out: [nCols]T = undefined;
            for (0..nCols) |i| {
                const col = self.getCol(i);
                out[i] = .zero;
                for (col, vec) |r, v| {
                    out[i] = A.add(A.mult(r, v), out[i]);
                }
            }
            return out;
        }
    };
}

test "SMatrix" {
    const Dual = zad.DualNumber(f64, 0);
    const Mat = SMatrix(Dual, 3, 3);
    const mat: Mat = .{ .v = [9]Dual{
        .one,  .one,  .zero,
        .one,  .zero, .one,
        .zero, .zero, .one,
    } };

    const vec = [3]Dual{ .promote(0.5), .promote(0.3), .promote(-1.1) };
    {
        const out = mat.multR(vec);
        try std.testing.expectApproxEqAbs(0.8, out[0].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.6, out[1].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-1.1, out[2].x, TEST_TOLERANCE);
    }
    {
        // Transposed
        const out = mat.multR_T(vec);
        try std.testing.expectApproxEqAbs(0.8, out[0].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(0.5, out[1].x, TEST_TOLERANCE);
        try std.testing.expectApproxEqAbs(-0.8, out[2].x, TEST_TOLERANCE);
    }
}

/// A column-major matrix.
pub fn Matrix(comptime T: type) type {
    return struct {
        const Self = @This();

        values: []T,
        n_rows: usize,
        n_cols: usize,

        pub fn init(allocator: std.mem.Allocator, n_cols: usize, n_rows: usize) !Self {
            const values = try allocator.alloc(T, n_cols * n_rows);
            errdefer allocator.free(values);
            return .{
                .values = values,
                .n_cols = n_cols,
                .n_rows = n_rows,
            };
        }

        pub fn deinit(self: *const Self, allocator: std.mem.Allocator) void {
            allocator.free(self.values);
        }

        /// Compute the index.
        inline fn indexOf(self: *const Self, row: usize, col: usize) usize {
            return self.n_rows * col + row;
        }

        pub fn getColumn(self: *Self, col: usize) []T {
            const start = self.n_rows * col;
            return self.values[start .. start + self.n_rows];
        }

        /// Get a pointer to a particular element.
        pub fn getPtr(self: *Self, row: usize, col: usize) *T {
            return &self.values[self.indexOf(row, col)];
        }

        /// Get an element.
        pub fn get(self: *const Self, row: usize, col: usize) T {
            return self.values[self.indexOf(row, col)];
        }

        /// Sum over all columns. Returns an array with length equal to the
        /// number of columns. Caller owns the memory.
        pub fn sumColumns(self: *const Self, allocator: std.mem.Allocator) ![]T {
            const out = try allocator.alloc(T, self.n_cols);
            errdefer out.deinit(allocator);
            @memset(out, 0);

            for (0..self.n_cols) |col| {
                for (0..self.n_rows) |row| {
                    out[col] += self.values[self.indexOf(row, col)];
                }
            }

            return out;
        }

        /// Sum over all rows. Returns an array with length equal to the number
        /// of rows. Caller owns the memory.
        pub fn sumRows(self: *const Self, allocator: std.mem.Allocator) ![]T {
            const out = try allocator.alloc(T, self.n_cols);
            errdefer out.deinit(allocator);
            @memset(out, 0);

            for (0..self.n_rows) |row| {
                for (0..self.n_cols) |col| {
                    out[row] += self.values[self.indexOf(row, col)];
                }
            }

            return out;
        }
    };
}
