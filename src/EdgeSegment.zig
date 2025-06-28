const std = @import("std");
const zmsdf = @import("root.zig");

const util = @import("util.zig");
const equation_solver = @import("equation_solver.zig");

const EdgeSegment = @This();

kind: Kind,
color: zmsdf.EdgeColor,
points: [4]zmsdf.Vector2,

pub const cubic_search_starts = 4;
pub const cubic_search_steps = 4;

pub const Kind = enum(u8) {
    linear,
    quadratic,
    cubic,
};

const dot = zmsdf.Vector2.dot;
const mixVec2 = zmsdf.Vector2.lerp;

pub fn create2(
    p0: zmsdf.Vector2,
    p1: zmsdf.Vector2,
    color: zmsdf.EdgeColor,
) EdgeSegment {
    return .{
        .kind = .linear,
        .color = color,
        .points = .{
            p0,
            p1,
            .zero,
            .zero,
        },
    };
}

pub fn create3(
    p0: zmsdf.Vector2,
    p1: zmsdf.Vector2,
    p2: zmsdf.Vector2,
    color: zmsdf.EdgeColor,
) EdgeSegment {
    const cross = zmsdf.Vector2.cross(p1.subtract(p0), p2.subtract(p1));
    if (cross == 0.0) return create2(p0, p2, color);
    return .{
        .kind = .quadratic,
        .color = color,
        .points = .{
            p0,
            p1,
            p2,
            .zero,
        },
    };
}

pub fn create4(
    p0: zmsdf.Vector2,
    p1: zmsdf.Vector2,
    p2: zmsdf.Vector2,
    p3: zmsdf.Vector2,
    color: zmsdf.EdgeColor,
) EdgeSegment {
    var p12 = p2.subtract(p1);
    if (zmsdf.Vector2.cross(p1.subtract(p0), p12) == 0 and
        zmsdf.Vector2.cross(p12, p3.subtract(p2)) == 0)
    {
        return create2(p0, p3, color);
    }
    p12 = p1.multiplyByScalar(1.5).subtract(p0.multiplyByScalar(0.5));
    const p23 = p2.multiplyByScalar(1.5).subtract(p3.multiplyByScalar(0.5));
    if (p12.eql(p23)) {
        return create3(p0, p12, p3, color);
    }
    return .{
        .kind = .cubic,
        .color = color,
        .points = .{
            p0,
            p1,
            p2,
            p3,
        },
    };
}

pub fn distanceToPerpendicularDistance(self: EdgeSegment, origin: zmsdf.Vector2, param: f64) ?zmsdf.SignedDistance {
    if (param < 0) {
        const dir = self.direction(0).normalize(false);
        const aq = origin.subtract(self.points[0]);
        const ts = dot(aq, dir);
        if (ts < 0) {
            const perpendicular_distance = zmsdf.Vector2.cross(aq, dir);
            return .{
                .distance = perpendicular_distance,
                .dot = 0,
            };
        }
    } else if (param > 1) {
        const dir = self.direction(1).normalize(false);
        const bq = origin.subtract(self.points[1]);
        const ts = dot(bq, dir);
        if (ts > 0) {
            const perpendicular_distance = zmsdf.Vector2.cross(bq, dir);
            return .{
                .distance = perpendicular_distance,
                .dot = 0,
            };
        }
    }

    return null;
}

pub fn point(self: EdgeSegment, param: f64) zmsdf.Vector2 {
    switch (self.kind) {
        .linear => {
            return self.points[0].lerp(self.points[1], param);
        },
        .quadratic => {
            const p01 = self.points[0].lerp(self.points[1], param);
            const p12 = self.points[1].lerp(self.points[2], param);
            return p01.lerp(p12, param);
        },
        .cubic => {
            const p01 = self.points[0].lerp(self.points[1], param);
            const p12 = self.points[1].lerp(self.points[2], param);
            const p23 = self.points[2].lerp(self.points[3], param);
            const p012 = p01.lerp(p12, param);
            const p123 = p12.lerp(p23, param);
            return p012.lerp(p123, param);
        },
    }
}

pub fn direction(self: EdgeSegment, param: f64) zmsdf.Vector2 {
    switch (self.kind) {
        .linear => {
            return self.points[1].subtract(self.points[0]);
        },
        .quadratic => {
            const tangent = self.points[1].subtract(self.points[0]).lerp(self.points[2].subtract(self.points[1]), param);
            if (tangent.isZero()) {
                return self.points[2].subtract(self.points[0]);
            }
            return tangent;
        },
        .cubic => {
            const tangent = self.points[1].subtract(self.points[0]).lerp(self.points[2].subtract(self.points[1]).lerp(self.points[3].subtract(self.points[2]), param), param);
            if (tangent.isZero()) {
                if (param == 0) return self.points[2].subtract(self.points[0]);
                if (param == 1) return self.points[3].subtract(self.points[1]);
            }
            return tangent;
        },
    }
}

pub fn directionChange(self: EdgeSegment, param: f64) zmsdf.Vector2 {
    switch (self.kind) {
        .linear => {
            return .zero;
        },
        .quadratic => {
            return self.points[2].subtract(self.points[1]).subtract(self.points[1].subtract(self.points[0]));
        },
        .cubic => {
            const p10 = self.points[1].subtract(self.points[0]);
            const p21 = self.points[2].subtract(self.points[1]);
            const p32 = self.points[3].subtract(self.points[2]);
            return p21.subtract(p10).lerp(p32.subtract(p21), param);
        },
    }
}

pub fn length(self: EdgeSegment) ?f64 {
    switch (self.kind) {
        .linear => {
            return self.points[1].subtract(self.points[0]).length();
        },
        .quadratic => {
            const ab = self.points[1].subtract(self.points[0]);
            const br = self.points[2].subtract(self.points[1]).subtract(ab);
            const abab = dot(ab, ab);
            const abbr = dot(ab, br);
            const brbr = dot(br, br);
            const abLen = @sqrt(abab);
            const brLen = @sqrt(brbr);
            const crs = zmsdf.Vector2.cross(ab, br);
            const h = @sqrt(abab + abbr + abbr + brbr);
            return (brLen * ((abbr + brbr) * h - abbr * abLen) +
                crs * crs * @log((brLen * h + abbr + brbr) / (brLen * abLen + abbr))) / (brbr * brLen);
        },
        .cubic => return null,
    }
}

pub fn signedDistance(self: EdgeSegment, origin: zmsdf.Vector2) struct {
    alpha: f64,
    distance: zmsdf.SignedDistance,
} {
    switch (self.kind) {
        .linear => {
            const aq = origin.subtract(self.points[0]);
            const ab = self.points[1].subtract(self.points[0]);
            const alpha = dot(aq, ab) / dot(ab, ab);
            const eq = self.points[if (alpha > 0.5) 1 else 0].subtract(origin);
            const endpoint_distance = eq.length();
            if (alpha > 0 and alpha < 1) {
                const ortho_distance = dot(ab.getOrthonormal(false, false), aq);
                if (@abs(ortho_distance) < endpoint_distance) {
                    return .{
                        .alpha = alpha,
                        .distance = .init(ortho_distance, 0),
                    };
                }
            }
            return .{
                .alpha = alpha,
                .distance = .init(
                    util.nonZeroSign(f64, zmsdf.Vector2.cross(aq, ab)) * endpoint_distance,
                    @abs(dot(ab.normalize(false), eq.normalize(false))),
                ),
            };
        },
        .quadratic => {
            const qa = self.points[0].subtract(origin);
            const ab = self.points[1].subtract(self.points[0]);
            const br = self.points[2].subtract(self.points[1]).subtract(ab);
            const a = dot(br, br);
            const b = 3 * dot(ab, br);
            const c = 2 * dot(ab, ab) + dot(qa, br);
            const d = dot(qa, ab);

            const solutions = equation_solver.solveCubic(a, b, c, d);
            var ep_dir = self.direction(0.0);
            var min_distance = util.nonZeroSign(f64, zmsdf.Vector2.cross(ep_dir, qa)) * qa.length(); // distance from A

            var alpha = -dot(qa, ep_dir) / dot(ep_dir, ep_dir);
            {
                ep_dir = self.direction(1);
                const distance = (self.points[2].subtract(origin)).length(); // distance from B
                if (distance < @abs(min_distance)) {
                    min_distance = util.nonZeroSign(
                        f64,
                        zmsdf.Vector2.cross(ep_dir, self.points[2].subtract(origin)),
                    ) * distance;
                    alpha = dot(origin.subtract(self.points[1]), ep_dir) /
                        dot(ep_dir, ep_dir);
                }
            }

            for (solutions.constSlice()) |solution| {
                if (solution > 0 and solution < 1) {
                    const qe = qa.add(ab.multiplyByScalar(2 * solution)).add(br.multiplyByScalar(solution * solution));
                    const distance = qe.length();
                    if (distance <= @abs(min_distance)) {
                        min_distance = util.nonZeroSign(
                            f64,
                            zmsdf.Vector2.cross(ab.add(br.multiplyByScalar(solution)), qe),
                        ) * distance;
                        alpha = solution;
                    }
                }
            }

            if (alpha >= 0 and alpha <= 1) {
                return .{
                    .alpha = alpha,
                    .distance = .init(min_distance, 0),
                };
            } else if (alpha < 0.5) {
                return .{
                    .alpha = alpha,
                    .distance = .init(min_distance, @abs(dot(self.direction(0).normalize(false), qa.normalize(false)))),
                };
            } else {
                return .{
                    .alpha = alpha,
                    .distance = .init(min_distance, @abs(
                        dot(self.direction(1).normalize(false), (self.points[2].subtract(origin)).normalize(false)),
                    )),
                };
            }
        },
        .cubic => {
            const qa = self.points[0].subtract(origin);
            const ab = self.points[1].subtract(self.points[0]);
            const br = self.points[2].subtract(self.points[1]).subtract(ab);
            const as = self.points[3].subtract(self.points[2]).subtract(self.points[2].subtract(self.points[1])).subtract(br);

            var ep_dir = self.direction(0);
            var min_distance = util.nonZeroSign(f64, zmsdf.Vector2.cross(ep_dir, qa)) * qa.length(); // distance from A
            var alpha = -dot(qa, ep_dir) / dot(ep_dir, ep_dir);
            {
                ep_dir = self.direction(1);
                const distance = (self.points[3].subtract(origin)).length(); // distance from B
                if (distance < @abs(min_distance)) {
                    min_distance = util.nonZeroSign(
                        f64,
                        zmsdf.Vector2.cross(ep_dir, self.points[3].subtract(origin)),
                    ) * distance;
                    alpha = dot(ep_dir.subtract(self.points[3].subtract(origin)), ep_dir) / dot(ep_dir, ep_dir);
                }
            }

            for (0..cubic_search_starts) |i| {
                var t = @as(f64, @floatFromInt(i)) / @as(f64, cubic_search_starts);
                var qe = qa
                    .add(ab.multiplyByScalar(3 * t))
                    .add(br.multiplyByScalar(3 * t * t))
                    .add(as.multiplyByScalar(t * t * t));

                for (0..cubic_search_steps) |_| {
                    const d1 = ab.multiplyByScalar(3)
                        .add(br.multiplyByScalar(6 * t))
                        .add(as.multiplyByScalar(3 * t * t));
                    const d2 = br.multiplyByScalar(6)
                        .add(as.multiplyByScalar(6 * t));
                    t -= dot(qe, d1) / (dot(d1, d1) + dot(qe, d2));
                    if (t <= 0 or t >= 1) {
                        break;
                    }
                    qe = qa
                        .add(ab.multiplyByScalar(3 * t))
                        .add(br.multiplyByScalar(3 * t * t))
                        .add(as.multiplyByScalar(t * t * t));
                    const distance = qe.length();
                    if (distance < @abs(min_distance)) {
                        min_distance = util.nonZeroSign(
                            f64,
                            zmsdf.Vector2.cross(d1, qe),
                        ) * distance;
                        alpha = t;
                    }
                }
            }

            if (alpha >= 0 and alpha <= 1) {
                return .{
                    .alpha = alpha,
                    .distance = .init(min_distance, 0),
                };
            } else if (alpha < 0.5) {
                return .{
                    .alpha = alpha,
                    .distance = .init(min_distance, @abs(dot(self.direction(0).normalize(false), qa.normalize(false)))),
                };
            } else {
                return .{
                    .alpha = alpha,
                    .distance = .init(min_distance, @abs(dot(
                        self.direction(1).normalize(false),
                        (self.points[3].subtract(origin)).normalize(false),
                    ))),
                };
            }
        },
    }
}

pub fn scanlineIntersections(self: EdgeSegment, y: f64) std.BoundedArray(zmsdf.Scanline.Intersection, 3) {
    switch (self.kind) {
        .linear => {
            if ((y >= self.points[0].y and y < self.points[1].y) or (y >= self.points[1].y and y < self.points[0].y)) {
                const param = (y - self.points[0].y) / (self.points[1].y - self.points[0].y);
                const x0 = std.math.lerp(self.points[0].x, self.points[1].x, param);
                // dy[0] = if (self.points[1].y > self.points[0].y) 1 else -1;
                const dy_sign: zmsdf.Polarity = .of(self.points[1].y - self.points[0].y);
                return .{
                    .buffer = .{
                        .{ .x = x0, .direction = dy_sign.asInt() },
                        .invalid,
                        .invalid,
                    },
                    .len = 1,
                };
            }
            return .{};
        },
        .quadratic => {
            var x: [3]f64 = @splat(0);
            var dy_sign: [3]i32 = @splat(0);

            var total: usize = 0;
            var next_dy_sign: zmsdf.Polarity = .of(y - self.points[0].y);
            x[total] = self.points[0].x;
            if (self.points[0].y == y) {
                if (self.points[0].y < self.points[1].y or
                    (self.points[0].y == self.points[1].y and self.points[0].y < self.points[2].y))
                {
                    dy_sign[total] = 1;
                    total += 1;
                } else {
                    next_dy_sign = .pos;
                }
            }
            {
                const ab = self.points[1].subtract(self.points[0]);
                const br = self.points[2].subtract(self.points[1]).subtract(ab);
                var solutions = equation_solver.solveQuadratic(
                    br.y,
                    2 * ab.y,
                    self.points[0].y - y,
                );
                if (solutions.len >= 2 and solutions.buffer[0] > solutions.buffer[1]) {
                    std.mem.swap(f64, &solutions.buffer[0], &solutions.buffer[1]);
                }
                for (solutions.constSlice()) |solution| {
                    if (solution >= 0 and solution <= 1) {
                        x[total] = self.points[0].x + 2 * solution * ab.x + solution * solution * br.x;
                        if (next_dy_sign.asFloat() * (ab.y + solution * br.y) >= 0) {
                            dy_sign[total] = next_dy_sign.asInt();
                            total += 1;
                            next_dy_sign = next_dy_sign.invert();
                        }
                    }
                }
            }
            if (self.points[2].y == y) {
                if (next_dy_sign.isPositive() and total > 0) {
                    total -= 1;
                    next_dy_sign = .neg;
                }
                if ((self.points[2].y < self.points[1].y or
                    (self.points[2].y == self.points[1].y and
                        self.points[2].y < self.points[0].y)) and total < 2)
                {
                    x[total] = self.points[2].x;
                    if (next_dy_sign.isNegative()) {
                        dy_sign[total] = -1;
                        total += 1;
                        next_dy_sign = .pos;
                    }
                }
            }
            if (next_dy_sign != zmsdf.Polarity.of(y - self.points[2].y)) {
                if (total > 0) {
                    total -= 1;
                } else {
                    if (@abs(self.points[2].y - y) < @abs(self.points[0].y - y)) {
                        x[total] = self.points[2].x;
                    }
                    dy_sign[total] = next_dy_sign.asInt();
                    total += 1;
                }
            }

            return .{
                .buffer = .{
                    .{ .x = x[0], .direction = dy_sign[0] },
                    .{ .x = x[1], .direction = dy_sign[1] },
                    .{ .x = x[2], .direction = dy_sign[2] },
                },
                .len = total,
            };
        },
        .cubic => {
            var x: [3]f64 = @splat(0);
            var dy_sign: [3]i32 = @splat(0);

            var total: usize = 0;
            var next_dy_sign: zmsdf.Polarity = .of(y - self.points[0].y);
            x[total] = self.points[0].x;
            if (self.points[0].y == y) {
                if (self.points[0].y < self.points[1].y or
                    (self.points[0].y == self.points[1].y and
                        (self.points[0].y < self.points[2].y or
                            (self.points[0].y == self.points[2].y and
                                self.points[0].y < self.points[3].y))))
                {
                    dy_sign[total] = 1;
                    total += 1;
                } else {
                    next_dy_sign = .pos;
                }
            }
            {
                const ab = self.points[1].subtract(self.points[0]);
                const br = self.points[2].subtract(self.points[1]).subtract(ab);
                const as = self.points[3].subtract(self.points[2]).subtract(self.points[2].subtract(self.points[1])).subtract(br);
                var solutions = equation_solver.solveCubic(
                    as.y,
                    3 * br.y,
                    3 * ab.y,
                    self.points[0].y - y,
                );
                if (solutions.len >= 2) {
                    if (solutions.buffer[0] > solutions.buffer[1]) {
                        std.mem.swap(f64, &solutions.buffer[0], &solutions.buffer[1]);
                    }
                    if (solutions.len >= 3 and solutions.buffer[1] > solutions.buffer[2]) {
                        std.mem.swap(f64, &solutions.buffer[1], &solutions.buffer[2]);
                        if (solutions.buffer[0] > solutions.buffer[1]) {
                            std.mem.swap(f64, &solutions.buffer[0], &solutions.buffer[1]);
                        }
                    }
                }
                for (solutions.constSlice()) |solution| {
                    if (solution >= 0 and solution <= 1) {
                        x[total] = self.points[0].x + 3 * solution * ab.x + 3 * solution * solution * br.x + solution * solution * solution * as.x;
                        if (next_dy_sign.asFloat() * (ab.y + 2 * solution * br.y + solution * solution * as.y) >= 0) {
                            dy_sign[total] = next_dy_sign.asInt();
                            total += 1;
                            next_dy_sign = next_dy_sign.invert();
                        }
                    }
                }
            }
            if (self.points[3].y == y) {
                if (next_dy_sign.isPositive() and total > 0) {
                    total -= 1;
                    next_dy_sign = .neg;
                }
                if ((self.points[3].y < self.points[2].y or
                    (self.points[3].y == self.points[2].y and
                        (self.points[3].y < self.points[1].y or
                            (self.points[3].y == self.points[1].y and
                                self.points[3].y < self.points[0].y)))) and total < 2)
                {
                    x[total] = self.points[3].x;
                    if (next_dy_sign.isNegative()) {
                        dy_sign[total] = -1;
                        total += 1;
                        next_dy_sign = .pos;
                    }
                }
            }
            if (next_dy_sign != zmsdf.Polarity.of(y - self.points[3].y)) {
                if (total > 0) {
                    total -= 1;
                } else {
                    if (@abs(self.points[3].y - y) < @abs(self.points[0].y - y)) {
                        x[total] = self.points[3].x;
                    }
                    dy_sign[total] = next_dy_sign.asInt();
                    total += 1;
                }
            }
            return .{
                .buffer = .{
                    .{ .x = x[0], .direction = dy_sign[0] },
                    .{ .x = x[1], .direction = dy_sign[1] },
                    .{ .x = x[2], .direction = dy_sign[2] },
                },
                .len = total,
            };
        },
    }
}

const pointBounds = util.pointBounds;

pub fn bound(self: EdgeSegment, bounds: *zmsdf.Bounds) void {
    switch (self.kind) {
        .linear => {
            pointBounds(self.points[0], bounds);
            pointBounds(self.points[1], bounds);
        },
        .quadratic => {
            pointBounds(self.points[0], bounds);
            pointBounds(self.points[2], bounds);
            const bot = self.points[1].subtract(self.points[0]).subtract(self.points[2].subtract(self.points[1]));
            if (bot.x != 0) {
                const param = (self.points[1].x - self.points[0].x) / bot.x;
                if (param > 0 and param < 1) {
                    pointBounds(self.point(param), bounds);
                }
            }
            if (bot.y != 0) {
                const param = (self.points[1].y - self.points[0].y) / bot.y;
                if (param > 0 and param < 1) {
                    pointBounds(self.point(param), bounds);
                }
            }
        },
        .cubic => {
            pointBounds(self.points[0], bounds);
            pointBounds(self.points[3], bounds);
            const a0 = self.points[1].subtract(self.points[0]);
            const a1 = self.points[2].subtract(self.points[1]).subtract(a0).multiplyByScalar(2);
            const a2 = self.points[3]
                .subtract(self.points[2].multiplyByScalar(3))
                .add(self.points[1].multiplyByScalar(3))
                .subtract(self.points[0]);
            const solutions = equation_solver.solveQuadratic(
                a2.x,
                a1.x,
                a0.x,
            );
            for (solutions.constSlice()) |param| {
                if (param > 0 and param < 1) {
                    pointBounds(self.point(param), bounds);
                }
            }
            const solutions_y = equation_solver.solveQuadratic(
                a2.y,
                a1.y,
                a0.y,
            );
            for (solutions_y.constSlice()) |param| {
                if (param > 0 and param < 1) {
                    pointBounds(self.point(param), bounds);
                }
            }
        },
    }
}

pub fn reverse(self: EdgeSegment) EdgeSegment {
    switch (self.kind) {
        .linear => {
            return create2(self.points[1], self.points[0], self.color);
        },
        .quadratic => {
            return create3(self.points[2], self.points[1], self.points[0], self.color);
        },
        .cubic => {
            return create4(self.points[3], self.points[2], self.points[1], self.points[0], self.color);
        },
    }
}

pub fn moveStartPoint(self: EdgeSegment, to: zmsdf.Vector2) EdgeSegment {
    switch (self.kind) {
        .linear => {
            return create2(to, self.points[1], self.color);
        },
        .quadratic => {
            const orig_s_dir = self.points[0].subtract(self.points[1]);
            var p1 = self.points[1];
            p1 = p1.add((self.points[2].subtract(self.points[1])).multiplyByScalar(zmsdf.Vector2.cross(self.points[0].subtract(self.points[1]), to.subtract(self.points[0])) /
                zmsdf.Vector2.cross(self.points[0].subtract(self.points[1]), self.points[2].subtract(self.points[1]))));
            if (dot(orig_s_dir, to.subtract(p1)) < 0) {
                p1 = self.points[1];
            }
            return create3(to, p1, self.points[2], self.color);
        },
        .cubic => {
            var p1 = self.points[1];
            p1 = p1.add(to.subtract(self.points[0]));
            return create4(to, p1, self.points[2], self.points[3], self.color);
        },
    }
}

pub fn moveEndPoint(self: EdgeSegment, to: zmsdf.Vector2) EdgeSegment {
    switch (self.kind) {
        .linear => {
            return create2(self.points[0], to, self.color);
        },
        .quadratic => {
            const orig_e_dir = self.points[2].subtract(self.points[1]);
            var p1 = self.points[1];
            p1 = p1.add((to.subtract(self.points[2])).multiplyByScalar(zmsdf.Vector2.cross(orig_e_dir, to.subtract(self.points[1])) /
                zmsdf.Vector2.cross(orig_e_dir, self.points[0].subtract(self.points[1]))));
            if (zmsdf.Vector2.dot(orig_e_dir, self.points[2].subtract(p1)) < 0) {
                p1 = self.points[1];
            }
            return create3(self.points[0], p1, to, self.color);
        },
        .cubic => {
            var p2 = self.points[2];
            p2 = p2.add(to.subtract(self.points[3]));
            return create4(self.points[0], self.points[1], p2, to, self.color);
        },
    }
}

pub fn splitInThirds(self: EdgeSegment) [3]EdgeSegment {
    switch (self.kind) {
        .linear => {
            return .{
                create2(self.points[0], self.point(1.0 / 3.0), self.color),
                create2(self.point(1.0 / 3.0), self.point(2.0 / 3.0), self.color),
                create2(self.point(2.0 / 3.0), self.points[1], self.color),
            };
        },
        .quadratic => {
            return .{
                create3(
                    self.points[0],
                    self.points[0].lerp(self.points[1], 1.0 / 3.0),
                    self.point(1.0 / 3.0),
                    self.color,
                ),
                create3(
                    self.point(1.0 / 3.0),
                    mixVec2(
                        mixVec2(
                            self.points[0],
                            self.points[1],
                            5.0 / 9.0,
                        ),
                        mixVec2(
                            self.points[1],
                            self.points[2],
                            4.0 / 9.0,
                        ),
                        0.5,
                    ),
                    self.point(2.0 / 3.0),
                    self.color,
                ),
                create2(
                    self.point(2.0 / 3.0),
                    mixVec2(
                        self.points[1],
                        self.points[2],
                        2.0 / 3.0,
                    ),
                    self.color,
                ),
            };
        },
        .cubic => {
            const part0: EdgeSegment = .create4(
                self.points[0],
                if (self.points[0].eql(self.points[1]))
                    self.points[0]
                else
                    mixVec2(self.points[0], self.points[1], 1.0 / 3.0),
                mixVec2(
                    mixVec2(self.points[0], self.points[1], 1.0 / 3.0),
                    mixVec2(self.points[1], self.points[2], 1.0 / 3.0),
                    1.0 / 3.0,
                ),
                self.point(1.0 / 3.0),
                self.color,
            );
            const part1: EdgeSegment = .create4(
                self.point(1.0 / 3.0),
                mixVec2(
                    mixVec2(
                        mixVec2(self.points[0], self.points[1], 1.0 / 3.0),
                        mixVec2(self.points[1], self.points[2], 1.0 / 3.0),
                        1.0 / 3.0,
                    ),
                    mixVec2(
                        mixVec2(self.points[1], self.points[2], 1.0 / 3.0),
                        mixVec2(self.points[2], self.points[3], 1.0 / 3.0),
                        1.0 / 3.0,
                    ),
                    2.0 / 3.0,
                ),
                mixVec2(
                    mixVec2(
                        mixVec2(self.points[0], self.points[1], 2.0 / 3.0),
                        mixVec2(self.points[1], self.points[2], 2.0 / 3.0),
                        2.0 / 3.0,
                    ),
                    mixVec2(
                        mixVec2(self.points[1], self.points[2], 2.0 / 3.0),
                        mixVec2(self.points[2], self.points[3], 2.0 / 3.0),
                        2.0 / 3.0,
                    ),
                    1.0 / 3.0,
                ),
                self.point(2.0 / 3.0),
                self.color,
            );
            const part2: EdgeSegment = .create4(
                self.point(2.0 / 3.0),
                mixVec2(
                    mixVec2(self.points[1], self.points[2], 2.0 / 3.0),
                    mixVec2(self.points[2], self.points[3], 2.0 / 3.0),
                    2.0 / 3.0,
                ),
                if (self.points[2].eql(self.points[3]))
                    self.points[3]
                else
                    mixVec2(self.points[2], self.points[3], 2.0 / 3.0),
                self.points[3],
                self.color,
            );
            return .{
                part0,
                part1,
                part2,
            };
        },
    }
}

pub fn convertToCubic(self: EdgeSegment) EdgeSegment {
    switch (self.kind) {
        .linear => {
            // idk
            @panic("Cannot convert linear edge segment to cubic");
            // return undefined;
            // return create4(
            //     self.points[0],
            //     mixVec2(self.points[0], self.points[1], 1.0 / 3.0),
            //     mixVec2(self.points[1], self.points[0], 1.0 / 3.0),
            //     self.points[1],
            //     self.color,
            // );
        },
        .quadratic => {
            return create4(
                self.points[0],
                mixVec2(
                    self.points[0],
                    self.points[1],
                    2.0 / 3.0,
                ),
                mixVec2(
                    self.points[1],
                    self.points[2],
                    1.0 / 3.0,
                ),
                self.points[2],
                self.color,
            );
        },
        .cubic => {
            return self; // Already cubic
        },
    }
}
