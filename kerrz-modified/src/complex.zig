const std = @import("std");

/// Represents a complex number in complex exponential format, that is, through
/// a magnitude and an argument angle.
pub fn ComplexNumber(comptime T: type) type {
    return struct {
        /// How far a value can be into the imaginary plane before it's
        /// considered no longer a pure real number.
        const rounding_tolerance = 1e-8;

        const Self = @This();
        mag: T,
        arg: T,

        /// Adapt the dual number type.
        pub fn adapt(self: Self, comptime NewT: type) ComplexNumber(NewT) {
            return .{ .mag = .adaptFrom(self.mag), .arg = .adaptFrom(self.arg) };
        }

        pub fn fromReIm(re: T, im: T) Self {
            const A = T.Algebra;
            const mag = A.sqrt(A.add(A.powi(re, 2), A.powi(im, 2)));
            const arg = A.atan2(im, re);
            return .{ .mag = mag, .arg = arg };
        }

        /// Calculate the real component.
        pub fn real(self: Self) T {
            const A = T.Algebra;
            return A.mult(self.mag, A.cos(self.arg));
        }

        /// Calculate the imaginary component.
        pub fn imag(self: Self) T {
            const A = T.Algebra;
            return A.mult(self.mag, A.sin(self.arg));
        }

        /// Add to the real component.
        pub fn addRe(self: Self, re: T) Self {
            const A = T.Algebra;
            return .fromReIm(A.add(re, self.real()), self.imag());
        }

        /// Subtract one complex number from another, returning a - b.
        pub fn sub(a: Self, b: Self) Self {
            const A = T.Algebra;
            return .fromReIm(
                A.sub(a.real(), b.real()),
                A.sub(a.imag(), b.imag()),
            );
        }

        /// Rotate the complex value in the imaginary plane.
        pub fn rotateArg(self: Self, angle: T) Self {
            const A = T.Algebra;
            return .{ .mag = self.mag, .arg = A.add(self.arg, angle) };
        }

        /// Negate the complex value, sending re to -re, and im to -im.
        pub fn neg(self: Self) Self {
            return self.rotateArg(.promote(std.math.pi));
        }

        /// Take the square root of the complex value. Since there is an
        /// ambiguity of 180 degrees in the argument, the `neg` function can be
        /// used to obtain the other root.
        pub fn sqrt(self: Self) Self {
            const A = T.Algebra;
            return .{
                .mag = A.sqrt(self.mag),
                .arg = A.div(self.arg, .promote(2)),
            };
        }

        /// Take the cube root of the complex value. Since there is an
        /// ambiguity of 120 degrees in the argument, the `rotateArg` can be
        /// used to obtain the other roots.
        pub fn cuberoot(self: Self) Self {
            const A = T.Algebra;
            return .{
                .mag = A.cuberoot(self.mag),
                .arg = A.mult(self.arg, .promote(1.0 / 3.0)),
            };
        }

        /// Returns true if the imaginary component of the number is within
        /// some tolerance of zero.
        pub fn isReal(self: Self) bool {
            const error_tolerance = @sqrt(std.math.floatEps(T.T));
            return std.math.approxEqAbs(T.T, self.imag().x, 0.0, error_tolerance);
        }
    };
}
