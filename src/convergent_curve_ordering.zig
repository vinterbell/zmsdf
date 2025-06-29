const std = @import("std");
const zmsdf = @import("root.zig");

/// For curves a, b converging at P = a->point(1) = b->point(0) with the same (opposite) direction,
/// determines the relative ordering in which they exit P (i.e. whether a is to the left or right of
/// b at the smallest positive radius around P)
pub fn convergentCurveOrdering(a: *const zmsdf.EdgeSegment, b: *const zmsdf.EdgeSegment) zmsdf.Polarity {
    var control_points: [12]zmsdf.Vector2 = undefined;
    const corner: []zmsdf.Vector2 = control_points[4..];
    const a_cp_tmp: []zmsdf.Vector2 = control_points[8..];
    const a_order: u32 = a.kind.order();
    const b_order: u32 = b.kind.order();

    if (!(a_order >= 1 and a_order <= 3 and b_order >= 1 and b_order <= 3)) {
        // Not implemented - only linear, quadratic, and cubic curves supported
        return .zero;
    }

    @memcpy(a_cp_tmp, &a.points);
    @memcpy(corner, &b.points);

    if (!a_cp_tmp[a_order].eql(corner[0])) {
        return .zero;
    }

    _ = simplifyDegenerateCurve(a_cp_tmp[0..a_order]);
    _ = simplifyDegenerateCurve(corner[0..b_order]);
    for (0..a_order) |i| {
        corner[i - a_order] = a_cp_tmp[i];
    }
    return convergentCurveOrderingCorners(&control_points, 4);
}

const cross = zmsdf.Vector2.cross;
const polarity = zmsdf.util.polarity;

// returns the new order
fn simplifyDegenerateCurve(control_points: []zmsdf.Vector2) ?u32 {
    if (control_points.len == 3 and (control_points[1].eql(control_points[0]) or control_points[1].eql(control_points[2])) and
        (control_points[2].eql(control_points[0]) or control_points[2].eql(control_points[3])))
    {
        control_points[1] = control_points[2];
        return 1;
    }
    if (control_points.len == 2 and (control_points[1].eql(control_points[0]) or control_points[1].eql(control_points[2]))) {
        control_points[1] = control_points[2];
        return 0;
    }
    if (control_points.len == 1 and control_points[0].eql(control_points[1])) {
        return 0;
    }
    return null;
}

fn convergentCurveOrderingCorners(corners: []const zmsdf.Vector2, point_index: usize) zmsdf.Polarity {
    const corner = corners[point_index];
    const control_points_before = point_index;
    const control_points_after = corners.len - point_index - 1;
    if (!(control_points_before > 0 and control_points_after > 0))
        return .zero;

    var a1: zmsdf.Vector2 = .zero;
    var a2: zmsdf.Vector2 = .zero;
    var a3: zmsdf.Vector2 = .zero;
    var b1: zmsdf.Vector2 = .zero;
    var b2: zmsdf.Vector2 = .zero;
    var b3: zmsdf.Vector2 = .zero;
    a1 = corners[point_index - 1].subtract(corner);
    b1 = corners[point_index + 1].subtract(corner);
    if (control_points_before >= 2)
        a2 = corners[point_index - 2].subtract(corners[point_index - 1]).subtract(a1);
    if (control_points_after >= 2)
        b2 = corners[point_index + 2].subtract(corners[point_index + 1]).subtract(b1);
    if (control_points_before >= 3) {
        a3 = corners[point_index - 3].subtract(corners[point_index - 2])
            .subtract(corners[point_index - 2].subtract(corners[point_index - 1])).subtract(a2);
        a2 = a2.multiplyByScalar(3);
    }
    if (control_points_after >= 3) {
        b3 = corners[point_index + 3].subtract(corners[point_index + 2])
            .subtract(corners[point_index + 2].subtract(corners[point_index + 1])).subtract(b2);
        b2 = b2.multiplyByScalar(3);
    }
    a1 = a1.multiplyByScalar(@floatFromInt(control_points_before));
    b1 = b1.multiplyByScalar(@floatFromInt(control_points_after));
    if (a1.isNonZero() and b1.isNonZero()) {
        const as = a1.length();
        const bs = b1.length();
        const d1 = as * cross(a1, b2) + bs * cross(a2, b1);
        if (d1 != 0) {
            return polarity(d1);
        }
        const d2 = as * as * cross(a1, b3) + as * bs * cross(a2, b2) + bs * bs * cross(a3, b1);
        if (d2 != 0) {
            return polarity(d2);
        }
        const d3 = as * cross(a2, b3) + bs * cross(a3, b2);
        if (d3 != 0) {
            return polarity(d3);
        }
        return polarity(cross(a3, b3));
    }
    // Degenerate curve after corner (control point after corner equals corner)
    var s: zmsdf.Polarity = .pos;
    if (a1.isNonZero()) {
        // Swap aN <-> bN and handle in if (b1)
        b1 = a1;
        a1 = b2;
        b2 = a2;
        a2 = a1;
        a1 = b3;
        b3 = a3;
        a3 = a1;
        s = .neg; // make sure to also flip output
    }
    // Degenerate curve before corner (control point before corner equals corner)
    if (b1.isNonZero()) {
        // Two-and-a-half-th derivative
        const d1 = cross(a3, b1);
        if (d1 != 0) {
            return s.multiply(polarity(d1));
        }
        // Third derivative
        const d2 = cross(a2, b2);
        if (d2 != 0) {
            return s.multiply(polarity(d2));
        }
        // Three-and-a-half-th derivative
        const d3 = cross(a3, b2);
        if (d3 != 0) {
            return s.multiply(polarity(d3));
        }
        // Fourth derivative
        const d4 = cross(a2, b3);
        if (d4 != 0) {
            return s.multiply(polarity(d4));
        }
        // Four-and-a-half-th derivative
        const d5 = cross(a3, b3);
        if (d5 != 0) {
            return s.multiply(polarity(d5));
        }
    }
    // Degenerate curves on both sides of the corner (control point before and after corner equals corner)
    {
        // Two-and-a-half-th derivative
        const d1 = @sqrt(a2.length()) * cross(a2, b3) + @sqrt(b2.length()) * cross(a3, b2);
        if (d1 != 0) {
            return polarity(d1);
        }
        // Third derivative
        return polarity(cross(a3, b3));
    }
}
