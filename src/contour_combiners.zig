const std = @import("std");
const zmsdf = @import("root.zig");

const util = @import("util.zig");

pub fn SimpleContourCombiner(comptime EdgeSelector: type) type {
    return struct {
        const Self = @This();

        edge_selector: EdgeSelector,

        pub const EdgeSelectorType = EdgeSelector;
        pub const DistanceType = EdgeSelector.DistanceType;

        pub fn init(allocator: std.mem.Allocator, shape: *const zmsdf.Shape) !Self {
            _ = allocator; // Not used in this implementation
            _ = shape;

            return .{
                .edge_selector = .init,
            };
        }

        pub fn deinit(_: *Self) void {
            // No resources to deinitialize in this simple combiner
        }

        pub fn reset(self: *Self, p: zmsdf.Vector2) void {
            self.edge_selector.reset(p);
        }

        pub fn edgeSelector(self: *Self, i: usize) *EdgeSelector {
            _ = i; // Not used in this implementation
            return &self.edge_selector;
        }

        pub fn distance(self: *const Self) EdgeSelector.DistanceType {
            return self.edge_selector.distance();
        }
    };
}

pub const SimpleTrueDistanceContourCombiner = SimpleContourCombiner(zmsdf.TrueDistanceSelector);
pub const SimplePerpendicularDistanceContourCombiner = SimpleContourCombiner(zmsdf.PerpendicularDistanceSelector);
pub const SimpleMultiDistanceContourCombiner = SimpleContourCombiner(zmsdf.MultiDistanceSelector);
pub const SimpleMultiAndTrueDistanceContourCombiner = SimpleContourCombiner(zmsdf.MultiAndTrueDistanceSelector);

pub fn OverlappingContourCombiner(comptime EdgeSelector: type) type {
    return struct {
        const Self = @This();

        p: zmsdf.Vector2,
        allocator: std.mem.Allocator,
        windings: std.ArrayListUnmanaged(zmsdf.Polarity),
        edge_selectors: std.ArrayListUnmanaged(EdgeSelector),

        pub const EdgeSelectorType = EdgeSelector;
        pub const DistanceType = EdgeSelector.DistanceType;

        pub fn init(allocator: std.mem.Allocator, shape: *const zmsdf.Shape) !Self {
            var windings: std.ArrayListUnmanaged(zmsdf.Polarity) = try .initCapacity(allocator, shape.contours.items.len);
            errdefer windings.deinit(allocator);
            var edge_selectors: std.ArrayListUnmanaged(EdgeSelector) = try .initCapacity(allocator, shape.contours.items.len);
            errdefer edge_selectors.deinit(allocator);

            for (shape.contours.items) |contour| {
                windings.appendAssumeCapacity(contour.winding());
                edge_selectors.appendAssumeCapacity(.init);
            }

            return .{
                .p = .zero,
                .allocator = allocator,
                .windings = windings,
                .edge_selectors = edge_selectors,
            };
        }

        pub fn deinit(self: *Self) void {
            self.windings.deinit(self.allocator);
            self.edge_selectors.deinit(self.allocator);
        }

        pub fn reset(self: *Self, p: zmsdf.Vector2) void {
            self.p = p;
            for (self.edge_selectors.items) |*selector| {
                selector.reset(p);
            }
        }

        pub fn edgeSelector(self: *Self, i: usize) *EdgeSelector {
            return &self.edge_selectors.items[i];
        }

        pub fn distance(self: *const Self) DistanceType {
            const contour_count = self.edge_selectors.items.len;
            var shape_edge_selector: EdgeSelector = .init;
            var inner_edge_selector: EdgeSelector = .init;
            var outer_edge_selector: EdgeSelector = .init;

            shape_edge_selector.reset(self.p);
            inner_edge_selector.reset(self.p);
            outer_edge_selector.reset(self.p);

            for (0..contour_count) |i| {
                const edge_distance = self.edge_selectors.items[i].distance();
                shape_edge_selector.merge(self.edge_selectors.items[i]);
                if (self.windings.items[i] == .pos and resolveDistance(edge_distance) >= 0) {
                    inner_edge_selector.merge(self.edge_selectors.items[i]);
                }
                if (self.windings.items[i] == .neg and resolveDistance(edge_distance) <= 0) {
                    outer_edge_selector.merge(self.edge_selectors.items[i]);
                }
            }

            const shape_distance = shape_edge_selector.distance();
            const inner_distance = inner_edge_selector.distance();
            const outer_distance = outer_edge_selector.distance();
            const inner_scalar_distance = resolveDistance(inner_distance);
            const outer_scalar_distance = resolveDistance(outer_distance);
            var d: DistanceType = undefined;
            initDistance(&d);

            var winding: zmsdf.Polarity = .zero;
            if (inner_scalar_distance >= 0 and @abs(inner_scalar_distance) <= @abs(outer_scalar_distance)) {
                d = inner_distance;
                winding = .pos;
                for (0..contour_count) |i| {
                    if (self.windings.items[i] == .pos) {
                        const contour_distance = self.edge_selectors.items[i].distance();
                        if (@abs(resolveDistance(contour_distance)) < @abs(outer_scalar_distance) and
                            resolveDistance(contour_distance) > resolveDistance(d))
                        {
                            d = contour_distance;
                        }
                    }
                }
            } else if (outer_scalar_distance <= 0 and @abs(outer_scalar_distance) < @abs(inner_scalar_distance)) {
                d = outer_distance;
                winding = .neg;
                for (0..contour_count) |i| {
                    if (self.windings.items[i] == .neg) {
                        const contour_distance = self.edge_selectors.items[i].distance();
                        if (@abs(resolveDistance(contour_distance)) < @abs(inner_scalar_distance) and
                            resolveDistance(contour_distance) < resolveDistance(d))
                        {
                            d = contour_distance;
                        }
                    }
                }
            } else {
                return shape_distance;
            }

            for (0..contour_count) |i| {
                if (self.windings.items[i] != winding) {
                    const contour_distance = self.edge_selectors.items[i].distance();
                    if (resolveDistance(contour_distance) * resolveDistance(d) >= 0 and
                        @abs(resolveDistance(contour_distance)) < @abs(resolveDistance(d)))
                    {
                        d = contour_distance;
                    }
                }
            }
            if (resolveDistance(d) == resolveDistance(shape_distance)) {
                d = shape_distance;
            }
            return d;
        }
    };
}

pub const OverlappingTrueDistanceContourCombiner = OverlappingContourCombiner(zmsdf.TrueDistanceSelector);
pub const OverlappingPerpendicularDistanceContourCombiner = OverlappingContourCombiner(zmsdf.PerpendicularDistanceSelector);
pub const OverlappingMultiDistanceContourCombiner = OverlappingContourCombiner(zmsdf.MultiDistanceSelector);
pub const OverlappingMultiAndTrueDistanceContourCombiner = OverlappingContourCombiner(zmsdf.MultiAndTrueDistanceSelector);

fn initDistanceFloat(distance: *f64) void {
    distance.* = util.f64_min;
}

fn initDistanceMulti(distance: *zmsdf.MultiDistance) void {
    initDistanceFloat(&distance.r);
    initDistanceFloat(&distance.g);
    initDistanceFloat(&distance.b);
}

fn initDistanceMultiAndTrue(distance: *zmsdf.MultiAndTrueDistance) void {
    initDistanceFloat(&distance.r);
    initDistanceFloat(&distance.g);
    initDistanceFloat(&distance.b);
    initDistanceFloat(&distance.a);
}

fn initDistance(distance: anytype) void {
    return switch (@TypeOf(distance)) {
        *zmsdf.MultiDistance => initDistanceMulti(distance),
        *zmsdf.MultiAndTrueDistance => initDistanceMultiAndTrue(distance),
        *f64 => initDistanceFloat(distance),
        else => @compileError("Unsupported distance type: " ++ @typeName(@TypeOf(distance))),
    };
}

fn resolveDistanceFloat(distance: f64) f64 {
    return distance;
}

fn resolveDistanceMulti(distance: zmsdf.MultiDistance) f64 {
    return util.median(f64, distance.r, distance.g, distance.b);
}

fn resolveDistanceMultiAndTrue(distance: zmsdf.MultiAndTrueDistance) f64 {
    return distance.a;
}

fn resolveDistance(distance: anytype) f64 {
    return switch (@TypeOf(distance)) {
        zmsdf.MultiDistance => resolveDistanceMulti(distance),
        zmsdf.MultiAndTrueDistance => resolveDistanceMultiAndTrue(distance),
        f64 => resolveDistanceFloat(distance),
        else => @compileError("Unsupported distance type: " ++ @typeName(@TypeOf(distance))),
    };
}
