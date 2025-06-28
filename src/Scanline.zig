//! Represents a horizontal scanline intersecting a shape.
const Scanline = @This();

const std = @import("std");

const util = @import("util.zig");

intersections: []Intersection = &.{},
last_index: i32 = 0,

pub const Intersection = struct {
    /// X coordinate of the intersection.
    x: f64,
    /// Normalized Y direction of the oriented edge at the point of intersection.
    direction: i32,

    pub const invalid: Intersection = .{
        .x = util.f64_min,
        .direction = 0,
    };
};

/// Fill rule dictates how intersection total is interpreted during rasterization.
pub const FillRule = enum {
    nonzero,
    odd, // "even-odd"
    positive,
    negative,

    pub fn interpret(rule: FillRule, intersections: i32) bool {
        return switch (rule) {
            .nonzero => intersections != 0,
            .odd => intersections & 1 != 0,
            .positive => intersections > 0,
            .negative => intersections < 0,
        };
    }
};

pub fn preprocess(self: *Scanline) void {
    self.last_index = 0;
    if (self.intersections.len > 0) {
        std.sort.insertion(
            Intersection,
            self.intersections,
            {},
            compareIntersections,
        );
        var total_direction: i32 = 0;
        for (self.intersections) |*intersection| {
            total_direction += intersection.direction;
            intersection.direction = total_direction;
        }
    }
}

pub fn compareIntersections(_: void, a: Intersection, b: Intersection) bool {
    return util.polarity(a.x - b.x) == .pos;
}

pub fn moveTo(self: *Scanline, x: f64) i32 {
    if (self.intersections.len == 0)
        return -1;
    var index = self.last_index;
    if (x < self.intersections[@intCast(index)].x) {
        while (true) {
            if (index == 0) {
                self.last_index = 0;
                return -1;
            }
            index -= 1;
            if (x >= self.intersections[@intCast(index)].x) {
                break;
            }
        }
    } else {
        while (index < self.intersections.len - 1 and x >= self.intersections[@intCast(index + 1)].x) {
            index += 1;
        }
    }
    self.last_index = index;
    return index;
}

pub fn countIntersections(self: *Scanline, x: f64) i32 {
    return self.moveTo(x) + 1;
}

pub fn sumIntersections(self: *Scanline, x: f64) i32 {
    const index = self.moveTo(x);
    if (index >= 0)
        return self.intersections[@intCast(index)].direction;
    return 0;
}

pub fn filled(self: *Scanline, x: f64, rule: FillRule) bool {
    return rule.interpret(self.sumIntersections(x));
}

pub fn overlap(
    a: *const Scanline,
    b: *const Scanline,
    from_x: f64,
    to_x: f64,
    rule: FillRule,
) f64 {
    var total: f64 = 0.0;
    var a_inside = false;
    var b_inside = false;
    var ai: usize = 0;
    var bi: usize = 0;
    var ax = if (a.intersections.len > 0) a.intersections[ai].x else to_x;
    var bx = if (b.intersections.len > 0) b.intersections[bi].x else to_x;
    while (ax < from_x or bx < from_x) {
        const x_next = @min(ax, bx);
        if (ax == x_next and ai < a.intersections.len) {
            a_inside = rule.interpret(a.intersections[ai].direction);
            ai += 1;
            ax = if (ai < a.intersections.len) a.intersections[ai].x else to_x;
        }
        if (bx == x_next and bi < b.intersections.len) {
            b_inside = rule.interpret(b.intersections[bi].direction);
            bi += 1;
            bx = if (bi < b.intersections.len) b.intersections[bi].x else to_x;
        }
    }
    var x = from_x;
    while (ax < to_x or bx < to_x) {
        const x_next = @min(ax, bx);
        if (a_inside == b_inside) {
            total += x_next - x;
        }
        if (ax == x_next and ai < a.intersections.len) {
            a_inside = rule.interpret(a.intersections[ai].direction);
            ai += 1;
            ax = if (ai < a.intersections.len) a.intersections[ai].x else to_x;
        }
        if (bx == x_next and bi < b.intersections.len) {
            b_inside = rule.interpret(b.intersections[bi].direction);
            bi += 1;
            bx = if (bi < b.intersections.len) b.intersections[bi].x else to_x;
        }
        x = x_next;
    }
    if (a_inside == b_inside) {
        total += to_x - x;
    }
    return total;
}
