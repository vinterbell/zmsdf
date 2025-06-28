const std = @import("std");
const zmsdf = @import("root.zig");

pub fn interpolate(
    comptime T: type,
    comptime Channels: u8,
    output: []T,
    bitmap: zmsdf.BitmapConstRef(T, Channels),
    pos: zmsdf.Vector2,
) void {
    const p = pos.subtract(.init(0.5, 0.5));
    const l: i32 = @intFromFloat(@floor(p.x));
    const b: i32 = @intFromFloat(@floor(p.y));
    const r: i32 = l + 1;
    const t: i32 = b + 1;
    const lr: f64 = p.x - @as(f32, @floatFromInt(l));
    const bt: f64 = p.y - @as(f32, @floatFromInt(b));
    const clamped_l = std.math.clamp(l, 0, @as(i32, bitmap.width) - 1);
    const clamped_r = std.math.clamp(r, 0, @as(i32, bitmap.width) - 1);
    const clamped_b = std.math.clamp(b, 0, @as(i32, bitmap.height) - 1);
    const clamped_t = std.math.clamp(t, 0, @as(i32, bitmap.height) - 1);

    const lb_bitmap = bitmap.get(clamped_l, clamped_b);
    const rb_bitmap = bitmap.get(clamped_r, clamped_b);
    const lt_bitmap = bitmap.get(clamped_l, clamped_t);
    const rt_bitmap = bitmap.get(clamped_r, clamped_t);
    for (output, 0..) |*out, i| {
        const bottom = std.math.lerp(
            lb_bitmap[i],
            rb_bitmap[i],
            lr,
        );
        const top = std.math.lerp(
            lt_bitmap[i],
            rt_bitmap[i],
            lr,
        );
        out.* = std.math.lerp(
            bottom,
            top,
            bt,
        );
    }
}

pub fn nonZeroSign(
    comptime T: type,
    value: T,
) T {
    const signed = std.math.sign(value);
    if (signed == 0) {
        return 1;
    } else {
        return signed;
    }
}

pub fn pointBounds(p: zmsdf.Vector2, bounds: *zmsdf.Bounds) void {
    if (p.x < bounds.left) bounds.left = p.x;
    if (p.y < bounds.bottom) bounds.bottom = p.y;
    if (p.x > bounds.right) bounds.right = p.x;
    if (p.y > bounds.top) bounds.top = p.y;
}

pub fn median(comptime T: type, a: T, b: T, c: T) T {
    return @max(@min(a, b), @min(@max(a, b), c));
}

pub fn polarity(x: anytype) zmsdf.Polarity {
    const result = std.math.sign(x);
    if (result > 0) {
        return .pos;
    } else if (result < 0) {
        return .neg;
    } else {
        return .zero;
    }
}

pub const f64_max = std.math.floatMax(f64);
pub const f64_min = -std.math.floatMax(f64);
