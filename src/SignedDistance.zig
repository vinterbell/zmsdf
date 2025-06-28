//! Represents a signed distance and alignment, which together can be compared to uniquely determine the closest edge segment.
const SignedDistance = @This();

const std = @import("std");

const util = @import("util.zig");

distance: f64,
dot: f64,

pub fn init(distance: f64, dot: f64) SignedDistance {
    return .{ .distance = distance, .dot = dot };
}

pub const zero: SignedDistance = .{
    .distance = util.f64_min,
    .dot = 0,
};

pub inline fn lessThan(a: SignedDistance, b: SignedDistance) bool {
    return @abs(a.distance) < @abs(b.distance) or (@abs(a.distance) == @abs(b.distance) and a.dot < b.dot);
}

pub inline fn greaterThan(a: SignedDistance, b: SignedDistance) bool {
    return @abs(a.distance) > @abs(b.distance) or (@abs(a.distance) == @abs(b.distance) and a.dot > b.dot);
}

pub inline fn lessThanOrEqual(a: SignedDistance, b: SignedDistance) bool {
    return @abs(a.distance) < @abs(b.distance) or (@abs(a.distance) == @abs(b.distance) and a.dot <= b.dot);
}

pub inline fn greaterThanOrEqual(a: SignedDistance, b: SignedDistance) bool {
    return @abs(a.distance) > @abs(b.distance) or (@abs(a.distance) == @abs(b.distance) and a.dot >= b.dot);
}
