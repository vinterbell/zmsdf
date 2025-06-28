//! A 2-dimensional euclidean floating-point vector.
const Vector2 = @This();

const std = @import("std");

x: f64,
y: f64,

pub const zero: Vector2 = .{ .x = 0, .y = 0 };

pub inline fn init(x: f64, y: f64) Vector2 {
    return .{ .x = x, .y = y };
}

pub inline fn squaredLength(self: Vector2) f64 {
    return self.x * self.x + self.y * self.y;
}

pub inline fn length(self: Vector2) f64 {
    return @sqrt(self.squaredLength());
}

/// Returns the normalized vector - one that has the same direction but unit length.
/// `allow_zero` false by default
pub inline fn normalize(self: Vector2, allow_zero: bool) Vector2 {
    const len = self.length();
    if (len != 0) {
        return .{ .x = self.x / len, .y = self.y / len };
    }
    return .{ .x = 0, .y = if (allow_zero) 0 else 1 };
}

/// Returns a vector with the same length that is orthogonal to this one.
/// `polarity` true by default
pub inline fn getOrthogonal(self: Vector2, polarity: bool) Vector2 {
    if (polarity) {
        return .{ .x = -self.y, .y = self.x };
    } else {
        return .{ .x = self.y, .y = -self.x };
    }
}

/// Returns a vector with unit length that is orthogonal to this one.
/// `polarity` true by default
/// `allow_zero` false by default
pub inline fn getOrthonormal(self: Vector2, polarity: bool, allow_zero: bool) Vector2 {
    const len = self.length();
    if (len != 0) {
        if (polarity) {
            return .{ .x = -self.y / len, .y = self.x / len };
        } else {
            return .{ .x = self.y / len, .y = -self.x / len };
        }
    }
    const zero_value = if (allow_zero) 0 else 1;
    return .{ .x = 0, .y = if (polarity) zero_value else -zero_value };
}

pub inline fn isNonZero(self: Vector2) bool {
    return self.x != 0 or self.y != 0;
}

pub inline fn isZero(self: Vector2) bool {
    return self.x == 0 and self.y == 0;
}

pub inline fn add(self: Vector2, other: Vector2) Vector2 {
    return .{ .x = self.x + other.x, .y = self.y + other.y };
}

pub inline fn subtract(self: Vector2, other: Vector2) Vector2 {
    return .{ .x = self.x - other.x, .y = self.y - other.y };
}

pub inline fn multiply(self: Vector2, other: Vector2) Vector2 {
    return .{ .x = self.x * other.x, .y = self.y * other.y };
}

pub inline fn divide(self: Vector2, other: Vector2) Vector2 {
    return .{ .x = self.x / other.x, .y = self.y / other.y };
}

pub inline fn multiplyByScalar(self: Vector2, value: f64) Vector2 {
    return .{ .x = self.x * value, .y = self.y * value };
}

pub inline fn divideByScalar(self: Vector2, value: f64) Vector2 {
    return .{ .x = self.x / value, .y = self.y / value };
}

pub inline fn dot(a: Vector2, b: Vector2) f64 {
    return a.x * b.x + a.y * b.y;
}

pub inline fn cross(a: Vector2, b: Vector2) f64 {
    return a.x * b.y - a.y * b.x;
}

pub inline fn eql(a: Vector2, b: Vector2) bool {
    return a.x == b.x and a.y == b.y;
}

pub inline fn lerp(a: Vector2, b: Vector2, t: f64) Vector2 {
    return .{ .x = std.math.lerp(a.x, b.x, t), .y = std.math.lerp(a.y, b.y, t) };
}
