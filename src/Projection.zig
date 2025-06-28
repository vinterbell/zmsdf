//! A transformation from shape coordinates to pixel coordinates.
const Projection = @This();

const std = @import("std");
const zmsdf = @import("root.zig");

scale: zmsdf.Vector2,
translate: zmsdf.Vector2,

pub const default: Projection = .{
    .scale = .init(1, 1),
    .translate = zmsdf.Vector2.zero,
};

pub fn init(scale: zmsdf.Vector2, translate: zmsdf.Vector2) Projection {
    return .{ .scale = scale, .translate = translate };
}

/// Converts the shape coordinate to pixel coordinate.
pub fn project(self: Projection, coord: zmsdf.Vector2) zmsdf.Vector2 {
    return self.scale.multiply(coord.add(self.translate));
}

/// Converts the pixel coordinate to shape coordinate.
pub fn unproject(self: Projection, coord: zmsdf.Vector2) zmsdf.Vector2 {
    return coord.divide(self.scale).subtract(self.translate);
}

/// Converts the vector to pixel coordinate space.
pub fn projectVector(self: Projection, vector: zmsdf.Vector2) zmsdf.Vector2 {
    return self.scale.multiply(vector);
}

/// Converts the vector from pixel coordinate space.
pub fn unprojectVector(self: Projection, vector: zmsdf.Vector2) zmsdf.Vector2 {
    return vector.divide(self.scale);
}

/// Converts the X-coordinate from shape to pixel coordinate space.
pub fn projectX(self: Projection, x: f64) f64 {
    return self.scale.x * (x + self.translate.x);
}

/// Converts the Y-coordinate from shape to pixel coordinate space.
pub fn projectY(self: Projection, y: f64) f64 {
    return self.scale.y * (y + self.translate.y);
}

/// Converts the X-coordinate from pixel to shape coordinate space.
pub fn unprojectX(self: Projection, x: f64) f64 {
    return x / self.scale.x - self.translate.x;
}

/// Converts the Y-coordinate from pixel to shape coordinate space.
pub fn unprojectY(self: Projection, y: f64) f64 {
    return y / self.scale.y - self.translate.y;
}
