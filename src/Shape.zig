//! Vector shape representation.
//!
//! If you want to add contours, use the array list `contours`. `deinit` makes the assumption that everything uses the same allocator.
const Shape = @This();

const std = @import("std");
const zmsdf = @import("root.zig");

pub const corner_dot_epsilon: f64 = 0.000001;
// moves control points slightly more than necessary to account for floating-point errors
const deconverge_overshoot: f64 = 1.11111111111111111; // 1 + 1/9
const corner_dot_epsilon_minus_one: f64 = corner_dot_epsilon - 1;

const Bounds = zmsdf.Bounds;

allocator: std.mem.Allocator,
/// The list of contours the shape consists of.
contours: std.ArrayListUnmanaged(zmsdf.Contour),
/// Specifies whether the shape uses bottom-to-top (false) or top-to-bottom (true) Y coordinates.
inverse_y_axis: bool,

pub fn init(allocator: std.mem.Allocator) Shape {
    return .{
        .allocator = allocator,
        .contours = .empty,
        .inverse_y_axis = false,
    };
}

pub fn deinit(self: *Shape) void {
    for (self.contours.items) |*contour| {
        contour.deinit();
    }
    self.contours.deinit(self.allocator);
}

pub fn validate(self: Shape) bool {
    for (self.contours.items) |contour| {
        if (contour.edges.items.len == 0) continue;
        var corner = contour.edges.getLast().point(1.0);
        for (contour.edges.items) |edge| {
            if (!edge.point(0.0).eql(corner)) return false;
            corner = edge.point(1.0);
        }
    }
    return true;
}

fn deconvergeEdge(edge: *zmsdf.EdgeSegment, param: i32, vector: zmsdf.Vector2) void {
    switch (edge.kind) {
        .quadratic, .cubic => {
            if (edge.kind == .quadratic) {
                edge.* = edge.convertToCubic();
            }
            const p = &edge.points;
            switch (param) {
                0 => {
                    p[1] = p[1].add(vector.multiplyByScalar(p[1].subtract(p[0]).length()));
                },
                1 => {
                    p[2] = p[2].add(vector.multiplyByScalar(p[2].subtract(p[3]).length()));
                },
                else => {},
            }
        },
        else => {},
    }
}

pub fn normalize(shape: *Shape) !void {
    for (shape.contours.items) |*contour| {
        if (contour.edges.items.len == 1) {
            const edge = contour.edges.items[0];
            const parts = edge.splitInThirds();
            contour.edges.clearRetainingCapacity();
            try contour.edges.appendSlice(contour.allocator, &parts);
        } else {
            var prev_edge = &contour.edges.items[contour.edges.items.len - 1];
            for (contour.edges.items) |*edge| {
                const prev_dir = prev_edge.direction(1.0).normalize(false);
                const cur_dir = edge.direction(0.0).normalize(false);
                if (zmsdf.Vector2.dot(prev_dir, cur_dir) < corner_dot_epsilon_minus_one) {
                    const factor = deconverge_overshoot *
                        @sqrt(1 - (corner_dot_epsilon_minus_one) * (corner_dot_epsilon_minus_one)) / (corner_dot_epsilon_minus_one);
                    var axis = (cur_dir.subtract(prev_dir)).normalize(false).multiplyByScalar(factor);
                    if (zmsdf.convergentCurveOrdering(prev_edge, edge) == .neg) {
                        axis = axis.multiplyByScalar(-1);
                    }
                    deconvergeEdge(prev_edge, 1, axis.getOrthonormal(true, false));
                    deconvergeEdge(edge, 0, axis.getOrthonormal(false, false));
                }
                prev_edge = edge;
            }
        }
    }
}

pub fn bound(self: Shape, bounds: *Bounds) void {
    for (self.contours.items) |contour| {
        contour.bound(bounds);
    }
}

pub fn boundMiters(self: Shape, bounds: *Bounds, border: f64, miter_limit: f64, polarity: i32) void {
    for (self.contours.items) |contour| {
        contour.boundMiters(bounds, border, miter_limit, polarity);
    }
}

pub const large_value: f64 = 1e240;
pub fn getBounds(self: Shape, border: f64, miter_limit: f64, polarity: i32) Bounds {
    var bounds: Bounds = .{
        .left = large_value,
        .bottom = large_value,
        .right = -large_value,
        .top = -large_value,
    };
    self.bound(&bounds);
    if (border > 0) {
        bounds.left -= border;
        bounds.bottom -= border;
        bounds.right += border;
        bounds.top += border;
        if (miter_limit > 0) {
            self.boundMiters(&bounds, border, miter_limit, polarity);
        }
    }
    return bounds;
}

pub fn scanline(self: Shape, allocator: std.mem.Allocator, y: f64) !zmsdf.Scanline {
    var intersections: std.ArrayListUnmanaged(zmsdf.Scanline.Intersection) = .empty;
    errdefer intersections.deinit(allocator);

    for (self.contours.items) |contour| {
        for (contour.edges.items) |edge| {
            const ints = edge.scanlineIntersections(y);
            try intersections.appendSlice(allocator, ints.constSlice());
        }
    }
    return .{
        .intersections = intersections.items,
    };
}

pub fn edgeCount(self: Shape) usize {
    var total: usize = 0;
    for (self.contours.items) |contour| {
        total += contour.edges.items.len;
    }
    return total;
}

const ContourAndIntersection = struct {
    intersection: zmsdf.Scanline.Intersection,
    contour_index: usize,
};

pub fn compareContourIntersections(_: void, a: ContourAndIntersection, b: ContourAndIntersection) bool {
    return zmsdf.Scanline.compareIntersections({}, a.intersection, b.intersection);
}

/// any memory allocated here is temporary and will be freed after the function returns
/// maybe a good idea to use an arena or something similar
pub fn orientContours(self: *Shape, temp: std.mem.Allocator) !void {
    // An irrational number to minimize chance of intersecting a corner or other point of interest
    const ratio = 0.5 * (@sqrt(5.0) - 1);

    var orientations: std.ArrayListUnmanaged(i32) = try .initCapacity(temp, self.contours.items.len);
    defer orientations.deinit(temp);
    orientations.items.len = self.contours.items.len;
    @memset(orientations.items, 0);

    var intersections: std.ArrayListUnmanaged(ContourAndIntersection) = .empty;
    defer intersections.deinit(temp);

    for (self.contours.items, orientations.items) |*contour, orientation| {
        if (!(orientation == 0 and contour.edges.items.len > 0)) {
            continue;
        }
        // Find a Y that crosses the contour
        const y0 = contour.edges.items[0].point(0).y;
        var y1 = y0;
        for (contour.edges.items) |edge| {
            y1 = @max(y1, edge.point(1).y);
            if (y0 != y1) {
                break; // found a non-horizontal edge
            }
        }
        for (contour.edges.items) |edge| {
            y1 = @max(y1, edge.point(ratio).y); // in case all endpoints are in a horizontal line
            if (y0 != y1) {
                break; // found a non-horizontal edge
            }
        }
        const y = std.math.lerp(y0, y1, ratio);
        for (self.contours.items, 0..) |other_contour, other_contour_index| {
            for (other_contour.edges.items) |edge| {
                const n = edge.scanlineIntersections(y);
                for (n.constSlice()) |int| {
                    try intersections.append(temp, .{
                        .intersection = int,
                        .contour_index = other_contour_index,
                    });
                }
            }
        }
        if (intersections.items.len == 0) {
            continue; // no intersections found
        }
        std.sort.insertion(
            ContourAndIntersection,
            intersections.items,
            {},
            compareContourIntersections,
        );
        // Disqualify multiple intersections
        for (1..intersections.items.len) |j| {
            if (intersections.items[j].intersection.x == intersections.items[j - 1].intersection.x) {
                intersections.items[j].intersection.direction = 0;
                intersections.items[j - 1].intersection.direction = 0;
            }
        }
        // Inspect scanline and deduce orientations of intersected contours
        for (intersections.items, 0..) |*intersection, j| {
            if (intersection.intersection.direction != 0) {
                const changed = 2 * ((@as(i32, @intCast(j)) & 1) ^ @intFromBool(intersection.intersection.direction > 0)) - 1;
                orientations.items[intersection.contour_index] += changed;
            }
        }
        intersections.clearRetainingCapacity();
    }

    // Reverse contours that have the opposite orientation
    for (self.contours.items, orientations.items) |*contour, orientation| {
        if (orientation < 0) {
            contour.reverse();
        }
    }
}
