const std = @import("std");
const zmsdf = @import("root.zig");

/// Finds the distance between a point and a Shape. ContourCombiner dictates the distance metric and its data type.
pub fn ShapeDistanceFinder(
    comptime ContourCombiner: type,
) type {
    return struct {
        allocator: std.mem.Allocator,
        shape: *const zmsdf.Shape,
        combiner: ContourCombiner,
        shape_edge_cache: []ContourCombiner.EdgeSelectorType.EdgeCache,

        pub const Self = @This();
        pub const DistanceType = ContourCombiner.DistanceType;

        /// Passed shape object must persist until the distance finder is destroyed!
        pub fn init(allocator: std.mem.Allocator, shape: *const zmsdf.Shape) !Self {
            const shape_edge_cache = try allocator.alloc(ContourCombiner.EdgeSelectorType.EdgeCache, shape.edgeCount());
            return .{
                .allocator = allocator,
                .shape = shape,
                .combiner = try .init(allocator, shape),
                .shape_edge_cache = shape_edge_cache,
            };
        }

        pub fn deinit(self: *Self) void {
            self.allocator.free(self.shape_edge_cache);
            self.combiner.deinit();
        }

        /// Finds the distance from origin. Not thread-safe! Is fastest when subsequent queries are close together.
        pub fn distance(self: *Self, origin: zmsdf.Vector2) DistanceType {
            self.combiner.reset(origin);
            var edge_cache = self.shape_edge_cache;

            for (self.shape.contours.items, 0..) |contour, i| {
                if (contour.edges.items.len > 0) {
                    const edge_selector = self.combiner.edgeSelector(i);

                    var prev_edge: *const zmsdf.EdgeSegment = if (contour.edges.items.len >= 2)
                        &contour.edges.items[contour.edges.items.len - 2]
                    else
                        &contour.edges.items[0];
                    var cur_edge: *const zmsdf.EdgeSegment = &contour.edges.items[contour.edges.items.len - 1];
                    for (0..contour.edges.items.len) |edge_index| {
                        const next_edge = &contour.edges.items[edge_index];
                        edge_selector.addEdge(&edge_cache[0], prev_edge, cur_edge, next_edge);
                        edge_cache = edge_cache[1..];
                        prev_edge = cur_edge;
                        cur_edge = next_edge;
                    }
                }
            }

            return self.combiner.distance();
        }

        // / Finds the distance between shape and origin. Does not allocate result cache used to optimize performance of multiple queries.
        // pub fn oneShotDistance(shape: *const zmsdf.Shape, origin: zmsdf.Vector2) DistanceType {}
    };
}
