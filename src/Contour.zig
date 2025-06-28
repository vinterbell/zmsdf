const Contour = @This();

const std = @import("std");
const zmsdf = @import("root.zig");

const util = @import("util.zig");

// please use this when modifying edges
allocator: std.mem.Allocator,
edges: std.ArrayListUnmanaged(zmsdf.EdgeSegment),

pub fn init(allocator: std.mem.Allocator) Contour {
    return .{
        .allocator = allocator,
        .edges = .empty,
    };
}

pub fn deinit(self: *Contour) void {
    self.edges.deinit(self.allocator);
}

pub fn bound(self: Contour, bounds: *zmsdf.Bounds) void {
    for (self.edges.items) |edge| {
        edge.bound(bounds);
    }
}

pub fn boundMiters(self: Contour, bounds: *zmsdf.Bounds, border: f64, miter_limit: f64, polarity: i32) void {
    if (self.edges.items.len == 0) return;
    var prev_dir = self.edges.getLast().direction(1.0).normalize(true);
    for (self.edges.items) |edge| {
        const dir = edge.direction(0.0).normalize(true).multiplyByScalar(-1);
        if (zmsdf.Vector2.cross(prev_dir, dir) * @as(f64, @floatFromInt(polarity)) >= 0) {
            var miter_len = miter_limit;
            const q = 0.5 * (1 - zmsdf.Vector2.dot(prev_dir, dir));
            if (q > 0) {
                miter_len = @min(1 / @sqrt(q), miter_limit);
            }
            const miter = edge.point(0.0).add((prev_dir.add(dir).normalize(true).multiplyByScalar(border * miter_len)));
            util.pointBounds(miter, bounds);
        }
        prev_dir = edge.direction(1.0).normalize(true);
    }
}

pub fn winding(self: Contour) zmsdf.Polarity {
    if (self.edges.items.len == 0) return .pos;
    var total: f64 = 0.0;
    if (self.edges.items.len == 1) {
        const a = self.edges.items[0].point(0);
        const b = self.edges.items[0].point(1.0 / 3.0);
        const c = self.edges.items[0].point(2.0 / 3.0);
        total += shoelace(a, b);
        total += shoelace(b, c);
        total += shoelace(c, a);
    } else if (self.edges.items.len == 2) {
        const a = self.edges.items[0].point(0);
        const b = self.edges.items[0].point(1.0 / 2.0);
        const c = self.edges.items[1].point(1.0 / 2.0);
        const d = self.edges.items[1].point(1.0);
        total += shoelace(a, b);
        total += shoelace(b, c);
        total += shoelace(c, d);
        total += shoelace(d, a);
    } else {
        var prev = self.edges.getLast().point(0);
        for (self.edges.items) |edge| {
            const cur = edge.point(0);
            total += shoelace(prev, cur);
            prev = cur;
        }
    }

    return .of(total);
}

fn shoelace(a: zmsdf.Vector2, b: zmsdf.Vector2) f64 {
    return (b.x - a.x) * (b.y + a.y);
}

pub fn reverse(self: *Contour) void {
    std.mem.reverse(zmsdf.EdgeSegment, self.edges.items);
    for (self.edges.items) |*edge| {
        edge.* = edge.reverse();
    }
}
