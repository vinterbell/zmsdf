const std = @import("std");
const zmsdf = @import("root.zig");

const util = @import("util.zig");

pub const MultiDistance = struct {
    r: f64,
    g: f64,
    b: f64,
};

pub const MultiAndTrueDistance = struct {
    r: f64,
    g: f64,
    b: f64,
    a: f64,

    pub inline fn toMultiDistance(self: MultiAndTrueDistance) MultiDistance {
        return .{ .r = self.r, .g = self.g, .b = self.b };
    }
};

pub const distance_delta_factor = 1.001;

pub const TrueDistanceSelector = struct {
    point: zmsdf.Vector2,
    min_distance: zmsdf.SignedDistance,

    pub const DistanceType = f64;

    pub const init: TrueDistanceSelector = .{
        .point = .zero,
        .min_distance = .zero,
    };

    pub const EdgeCache = struct {
        point: zmsdf.Vector2,
        abs_distance: f64,

        pub fn init() EdgeCache {
            return .{ .point = .zero, .abs_distance = 0 };
        }
    };

    pub fn reset(self: *TrueDistanceSelector, p: zmsdf.Vector2) void {
        const delta = distance_delta_factor * (p.subtract(self.point)).length();
        self.min_distance.distance += util.nonZeroSign(f64, self.min_distance.distance) * delta;
        self.point = p;
    }

    pub fn addEdge(
        self: *TrueDistanceSelector,
        cache: *EdgeCache,
        prev_edge: *const zmsdf.EdgeSegment,
        edge: *zmsdf.EdgeSegment,
        next_edge: *const zmsdf.EdgeSegment,
    ) void {
        _ = prev_edge; // Unused
        _ = next_edge; // Unused
        const delta = distance_delta_factor * (self.point.subtract(cache.point)).length();
        if (cache.abs_distance - delta <= @abs(self.min_distance.distance)) {
            const alpha_and_distance = edge.signedDistance(self.point);
            if (alpha_and_distance.distance.lessThan(self.min_distance)) {
                self.min_distance = alpha_and_distance.distance;
            }
            cache.point = self.point;
            cache.abs_distance = @abs(alpha_and_distance.distance.distance);
        }
    }

    pub fn merge(self: *TrueDistanceSelector, other: TrueDistanceSelector) void {
        if (other.min_distance.lessThan(self.min_distance)) {
            self.min_distance = other.min_distance;
        }
    }

    pub fn distance(self: TrueDistanceSelector) DistanceType {
        return self.min_distance.distance;
    }
};

pub const PerpendicularDistanceSelectorBase = struct {
    min_true_distance: zmsdf.SignedDistance,
    min_negative_perpendicular_distance: f64,
    min_positive_perpendicular_distance: f64,
    near_edge: ?zmsdf.EdgeSegment,
    near_edge_param: f64,

    pub const EdgeCache = struct {
        point: zmsdf.Vector2,
        abs_distance: f64,
        a_domain_distance: f64,
        b_domain_distance: f64,
        a_perpendicular_distance: f64,
        b_perpendicular_distance: f64,

        pub const zero: EdgeCache = .{
            .point = .zero,
            .abs_distance = 0,
            .a_domain_distance = 0,
            .b_domain_distance = 0,
            .a_perpendicular_distance = 0,
            .b_perpendicular_distance = 0,
        };
    };

    pub const init: PerpendicularDistanceSelectorBase = .{
        .min_negative_perpendicular_distance = util.f64_min,
        .min_positive_perpendicular_distance = util.f64_max,
        .min_true_distance = .zero,
        .near_edge = null,
        .near_edge_param = 0,
    };

    /// if there is a perpendicular distance to the edge, it will returned
    pub fn getPerpendicularDistance(
        distance: *f64,
        ep: zmsdf.Vector2,
        edge_dir: zmsdf.Vector2,
    ) bool {
        const ts = zmsdf.Vector2.dot(ep, edge_dir);
        if (ts > 0) {
            const perpendicular_distance = zmsdf.Vector2.cross(ep, edge_dir);
            if (@abs(perpendicular_distance) < @abs(distance.*)) {
                distance.* = perpendicular_distance;
                return true;
            }
        }
        return false;
    }

    pub fn reset(self: *PerpendicularDistanceSelectorBase, delta: f64) void {
        self.min_true_distance.distance += util.nonZeroSign(f64, self.min_true_distance.distance) * delta;
        self.min_negative_perpendicular_distance = -@abs(self.min_true_distance.distance);
        self.min_positive_perpendicular_distance = @abs(self.min_true_distance.distance);
        self.near_edge = null;
        self.near_edge_param = 0;
    }

    pub fn isEdgeRelevant(
        self: *PerpendicularDistanceSelectorBase,
        cache: EdgeCache,
        edge: *const zmsdf.EdgeSegment,
        p: zmsdf.Vector2,
    ) bool {
        _ = edge; // Unused
        const delta = (p.subtract(cache.point).multiplyByScalar(distance_delta_factor)).length();
        return (cache.abs_distance - delta <= @abs(self.min_true_distance.distance) or
            @abs(cache.a_domain_distance) < delta or
            @abs(cache.b_domain_distance) < delta or
            (cache.a_domain_distance > 0 and (if (cache.a_perpendicular_distance < 0)
                cache.a_perpendicular_distance + delta >= self.min_negative_perpendicular_distance
            else
                cache.a_perpendicular_distance - delta <= self.min_positive_perpendicular_distance)) or
            (cache.b_domain_distance > 0 and (if (cache.b_perpendicular_distance < 0)
                cache.b_perpendicular_distance + delta >= self.min_negative_perpendicular_distance
            else
                cache.b_perpendicular_distance - delta <= self.min_positive_perpendicular_distance)));
    }

    pub fn addEdgeTrueDistance(
        self: *PerpendicularDistanceSelectorBase,
        edge: *const zmsdf.EdgeSegment,
        distance: zmsdf.SignedDistance,
        param: f64,
    ) void {
        if (distance.lessThan(self.min_true_distance)) {
            self.min_true_distance = distance;
            self.near_edge = edge.*;
            self.near_edge_param = param;
        }
    }

    pub fn addEdgePerpendicularDistance(
        self: *PerpendicularDistanceSelectorBase,
        distance: f64,
    ) void {
        if (distance <= 0 and distance > self.min_negative_perpendicular_distance) {
            self.min_negative_perpendicular_distance = distance;
        }
        if (distance >= 0 and distance < self.min_positive_perpendicular_distance) {
            self.min_positive_perpendicular_distance = distance;
        }
    }

    pub fn merge(self: *PerpendicularDistanceSelectorBase, other: PerpendicularDistanceSelectorBase) void {
        if (other.min_true_distance.lessThan(self.min_true_distance)) {
            self.min_true_distance = other.min_true_distance;
            self.near_edge = other.near_edge;
            self.near_edge_param = other.near_edge_param;
        }
        if (other.min_negative_perpendicular_distance > self.min_negative_perpendicular_distance) {
            self.min_negative_perpendicular_distance = other.min_negative_perpendicular_distance;
        }
        if (other.min_positive_perpendicular_distance < self.min_positive_perpendicular_distance) {
            self.min_positive_perpendicular_distance = other.min_positive_perpendicular_distance;
        }
    }

    pub fn computeDistance(self: PerpendicularDistanceSelectorBase, p: zmsdf.Vector2) f64 {
        var min_distance = if (self.min_true_distance.distance < 0)
            self.min_negative_perpendicular_distance
        else
            self.min_positive_perpendicular_distance;

        if (self.near_edge) |edge| {
            const distance = edge.distanceToPerpendicularDistance(p, self.near_edge_param) orelse
                self.min_true_distance;
            if (@abs(distance.distance) < @abs(min_distance)) {
                min_distance = distance.distance;
            }
        }
        return min_distance;
    }

    pub fn trueDistance(self: PerpendicularDistanceSelectorBase) zmsdf.SignedDistance {
        return self.min_true_distance;
    }
};

pub const PerpendicularDistanceSelector = struct {
    base: PerpendicularDistanceSelectorBase,
    point: zmsdf.Vector2,

    pub const DistanceType = f64;

    pub const init: PerpendicularDistanceSelector = .{
        .base = .init,
        .point = .zero,
    };

    pub fn reset(self: *PerpendicularDistanceSelector, p: zmsdf.Vector2) void {
        const delta = (p.subtract(self.point).multiplyByScalar(distance_delta_factor)).length();
        self.base.reset(delta);
        self.point = p;
    }

    pub fn addEdge(
        self: *PerpendicularDistanceSelector,
        cache: *PerpendicularDistanceSelectorBase.EdgeCache,
        prev_edge: *const zmsdf.EdgeSegment,
        edge: *const zmsdf.EdgeSegment,
        next_edge: *const zmsdf.EdgeSegment,
    ) void {
        if (self.base.isEdgeRelevant(cache.*, edge, self.point)) {
            const d = edge.signedDistance(self.point);
            self.base.addEdgeTrueDistance(edge, d.distance, d.alpha);
            cache.point = self.point;
            cache.abs_distance = @abs(d.distance.distance);

            const ap = self.point.subtract(edge.point(0));
            const bp = self.point.subtract(edge.point(1));
            const a_dir = edge.direction(0).normalize(true);
            const b_dir = edge.direction(1).normalize(true);
            const prev_dir = prev_edge.direction(1).normalize(true);
            const next_dir = next_edge.direction(0).normalize(true);
            const add = zmsdf.Vector2.dot(ap, (prev_dir.add(a_dir)).normalize(true));
            const bdd = -zmsdf.Vector2.dot(bp, (b_dir.add(next_dir)).normalize(true));

            if (add > 0) {
                var pd = d.distance.distance;
                if (PerpendicularDistanceSelectorBase.getPerpendicularDistance(&pd, ap, a_dir.multiplyByScalar(-1))) {
                    pd = -pd;
                    self.base.addEdgePerpendicularDistance(pd);
                }
                cache.a_perpendicular_distance = pd;
            }
            if (bdd > 0) {
                var pd = d.distance.distance;
                if (PerpendicularDistanceSelectorBase.getPerpendicularDistance(&pd, bp, b_dir)) {
                    self.base.addEdgePerpendicularDistance(pd);
                }
                cache.b_perpendicular_distance = pd;
            }
            cache.a_domain_distance = add;
            cache.b_domain_distance = bdd;
        }
    }

    pub fn distance(self: PerpendicularDistanceSelector) DistanceType {
        return self.base.computeDistance(self.point);
    }

    pub fn trueDistance(self: PerpendicularDistanceSelector) zmsdf.SignedDistance {
        return self.base.trueDistance();
    }

    pub fn merge(self: *PerpendicularDistanceSelector, other: PerpendicularDistanceSelector) void {
        self.base.merge(other.base);
    }
};

pub const MultiDistanceSelector = struct {
    r: PerpendicularDistanceSelectorBase,
    g: PerpendicularDistanceSelectorBase,
    b: PerpendicularDistanceSelectorBase,
    point: zmsdf.Vector2,

    pub const DistanceType = MultiDistance;
    pub const EdgeCache = PerpendicularDistanceSelectorBase.EdgeCache;

    pub const init: MultiDistanceSelector = .{
        .r = .init,
        .g = .init,
        .b = .init,
        .point = .zero,
    };

    pub fn reset(self: *MultiDistanceSelector, p: zmsdf.Vector2) void {
        const delta = (p.subtract(self.point).multiplyByScalar(distance_delta_factor)).length();
        self.r.reset(delta);
        self.g.reset(delta);
        self.b.reset(delta);
        self.point = p;
    }

    pub fn addEdge(
        self: *MultiDistanceSelector,
        cache: *EdgeCache,
        prev_edge: *const zmsdf.EdgeSegment,
        edge: *const zmsdf.EdgeSegment,
        next_edge: *const zmsdf.EdgeSegment,
    ) void {
        if ((edge.color.red_channel and self.r.isEdgeRelevant(cache.*, edge, self.point)) or
            (edge.color.green_channel and self.g.isEdgeRelevant(cache.*, edge, self.point)) or
            (edge.color.blue_channel and self.b.isEdgeRelevant(cache.*, edge, self.point)))
        {
            const d = edge.signedDistance(self.point);
            if (edge.color.red_channel)
                self.r.addEdgeTrueDistance(edge, d.distance, d.alpha);
            if (edge.color.green_channel)
                self.g.addEdgeTrueDistance(edge, d.distance, d.alpha);
            if (edge.color.blue_channel)
                self.b.addEdgeTrueDistance(edge, d.distance, d.alpha);
            cache.point = self.point;
            cache.abs_distance = @abs(d.distance.distance);

            const ap = self.point.subtract(edge.point(0));
            const bp = self.point.subtract(edge.point(1));
            const a_dir = edge.direction(0).normalize(true);
            const b_dir = edge.direction(1).normalize(true);
            const prev_dir = prev_edge.direction(1).normalize(true);
            const next_dir = next_edge.direction(0).normalize(true);
            const add = zmsdf.Vector2.dot(ap, (prev_dir.add(a_dir)).normalize(true));
            const bdd = -zmsdf.Vector2.dot(bp, (b_dir.add(next_dir)).normalize(true));

            if (add > 0) {
                var pd = d.distance.distance;
                if (PerpendicularDistanceSelectorBase.getPerpendicularDistance(&pd, ap, a_dir.multiplyByScalar(-1))) {
                    pd = -pd;
                    if (edge.color.red_channel)
                        self.r.addEdgePerpendicularDistance(pd);
                    if (edge.color.green_channel)
                        self.g.addEdgePerpendicularDistance(pd);
                    if (edge.color.blue_channel)
                        self.b.addEdgePerpendicularDistance(pd);
                }
                cache.a_perpendicular_distance = pd;
            }
            if (bdd > 0) {
                var pd = d.distance.distance;
                if (PerpendicularDistanceSelectorBase.getPerpendicularDistance(&pd, bp, b_dir)) {
                    if (edge.color.red_channel)
                        self.r.addEdgePerpendicularDistance(pd);
                    if (edge.color.green_channel)
                        self.g.addEdgePerpendicularDistance(pd);
                    if (edge.color.blue_channel)
                        self.b.addEdgePerpendicularDistance(pd);
                }
                cache.b_perpendicular_distance = pd;
            }
            cache.a_domain_distance = add;
            cache.b_domain_distance = bdd;
        }
    }

    pub fn merge(self: *MultiDistanceSelector, other: MultiDistanceSelector) void {
        self.r.merge(other.r);
        self.g.merge(other.g);
        self.b.merge(other.b);
    }

    pub fn distance(self: MultiDistanceSelector) DistanceType {
        return .{
            .r = self.r.computeDistance(self.point),
            .g = self.g.computeDistance(self.point),
            .b = self.b.computeDistance(self.point),
        };
    }

    pub fn trueDistance(self: MultiDistanceSelector) zmsdf.SignedDistance {
        var d = self.r.trueDistance();
        if (self.g.trueDistance().lessThan(d))
            d = self.g.trueDistance();
        if (self.b.trueDistance().lessThan(d))
            d = self.b.trueDistance();
        return d;
    }
};

pub const MultiAndTrueDistanceSelector = struct {
    base: MultiDistanceSelector,

    pub const DistanceType = MultiAndTrueDistance;
    pub const EdgeCache = MultiDistanceSelector.EdgeCache;

    pub const init: MultiAndTrueDistanceSelector = .{
        .base = .init,
    };

    pub fn distance(self: MultiAndTrueDistanceSelector) DistanceType {
        const multi_distance = self.base.distance();
        return .{
            .r = multi_distance.r,
            .g = multi_distance.g,
            .b = multi_distance.b,
            .a = self.base.trueDistance().distance,
        };
    }

    pub fn reset(self: *MultiAndTrueDistanceSelector, p: zmsdf.Vector2) void {
        self.base.reset(p);
    }
    pub fn addEdge(
        self: *MultiAndTrueDistanceSelector,
        cache: *MultiDistanceSelector.EdgeCache,
        prev_edge: *const zmsdf.EdgeSegment,
        edge: *const zmsdf.EdgeSegment,
        next_edge: *const zmsdf.EdgeSegment,
    ) void {
        self.base.addEdge(cache, prev_edge, edge, next_edge);
    }
    pub fn merge(self: *MultiAndTrueDistanceSelector, other: MultiAndTrueDistanceSelector) void {
        self.base.merge(other.base);
    }
    pub fn trueDistance(self: MultiAndTrueDistanceSelector) zmsdf.SignedDistance {
        return self.base.trueDistance();
    }
};
