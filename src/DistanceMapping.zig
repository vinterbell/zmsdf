//! Linear transformation of signed distance values.
const DistanceMapping = @This();

const zmsdf = @import("root.zig");

scale: f64,
translate: f64,

pub fn inverseRange(range: zmsdf.Range) DistanceMapping {
    const range_width = range.upper - range.lower;
    return .{
        .scale = range_width,
        .translate = range.lower / (if (range_width != 0) range_width else 1),
    };
}

pub fn fromRange(range: zmsdf.Range) DistanceMapping {
    return .{
        .scale = 1.0 / (range.upper - range.lower),
        .translate = -range.lower,
    };
}

pub fn of(self: DistanceMapping, distance: f64) f64 {
    return self.scale * (distance + self.translate);
}

pub fn ofDelta(self: DistanceMapping, delta: f64) f64 {
    return self.scale * delta;
}

pub fn inverse(self: DistanceMapping) DistanceMapping {
    return .{
        .scale = 1.0 / self.scale,
        .translate = -self.scale * self.translate,
    };
}
