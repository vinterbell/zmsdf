const std = @import("std");
const zmsdf = @import("root.zig");

const TrueType = zmsdf.TrueType;

/// The scaling applied to font glyph coordinates when loading a glyph
pub const FontCoordinateScaling = enum {
    /// The coordinates are kept as the integer values native to the font file
    none,
    /// The coordinates will be normalized to the em size, i.e. 1 = 1 em
    em_normalized,
};

pub fn getFontCoordinateScaling(
    ttf: *const TrueType,
    scaling: FontCoordinateScaling,
) f64 {
    return switch (scaling) {
        .none => 1.0,
        .em_normalized => if (ttf.units_per_em != 0) 1.0 / @as(f64, @floatFromInt(ttf.units_per_em)) else 1.0,
    };
}

fn msVector(
    x: i16,
    y: i16,
    scale: f64,
) zmsdf.Vector2 {
    return .{
        .x = @as(f64, @floatFromInt(x)) * scale,
        .y = @as(f64, @floatFromInt(y)) * scale,
    };
}

fn readTrueTypeOutline(
    shape: *zmsdf.Shape,
    outline: []const TrueType.Vertex,
    scale: f64,
) !void {
    var position: zmsdf.Vector2 = .zero;

    var contour_start_ends: [256]struct { usize, usize } = @splat(undefined);
    var contour_count: usize = 0;

    for (outline) |vertex| {
        if (vertex.type == .vmove) {
            contour_count += 1;
        }
    }

    var ci: usize = 0;
    for (0..outline.len + 1) |i| {
        if (i >= outline.len) {
            contour_start_ends[ci].@"1" = i;
        } else if (outline[i].type == .vmove) {
            if (i > 0) {
                contour_start_ends[ci].@"1" = i;
                ci += 1;
            }
            contour_start_ends[ci].@"0" = i;
        }
    }

    for (0..contour_count) |contour_index| {
        var contour: zmsdf.Contour = .init(shape.allocator);
        const start, const end = contour_start_ends[contour_index];
        for (start..end) |vertex_index| {
            const vertex = outline[vertex_index];
            const next_point = msVector(vertex.x, vertex.y, scale);
            const control_point_1 = msVector(vertex.cx, vertex.cy, scale);
            const control_point_2 = msVector(vertex.cx1, vertex.cy1, scale);

            switch (vertex.type) {
                .vmove => {
                    position = next_point;
                },
                .vline => {
                    try contour.edges.append(
                        contour.allocator,
                        .create2(position, next_point, .white),
                    );
                    position = next_point;
                },
                .vcurve => {
                    const pos_same_control = position.eql(control_point_1);
                    const control_same_next = control_point_1.eql(next_point);
                    const new_control = if (pos_same_control or control_same_next)
                        position.add(next_point).multiplyByScalar(0.5)
                    else
                        control_point_1;
                    try contour.edges.append(
                        contour.allocator,
                        .create3(position, new_control, next_point, .white),
                    );
                    position = next_point;
                },
                .vcubic => {
                    try contour.edges.append(
                        contour.allocator,
                        .create4(position, control_point_1, control_point_2, next_point, .white),
                    );
                    position = next_point;
                },
                else => @panic("Unsupported vertex type"),
            }
        }
        if (contour.edges.items.len > 0) {
            try shape.contours.append(shape.allocator, contour);
        } else {
            contour.deinit();
        }
    }
}

pub fn loadGlyph(
    ttf: *const TrueType,
    temp: std.mem.Allocator,
    shape: *zmsdf.Shape,
    glyph_index: TrueType.GlyphIndex,
    scaling: FontCoordinateScaling,
) !void {
    shape.inverse_y_axis = false;
    var arena: std.heap.ArenaAllocator = .init(temp);
    defer arena.deinit();

    const allocator = arena.allocator();

    const scale = getFontCoordinateScaling(ttf, scaling);
    const outline = try ttf.glyphShape(allocator, glyph_index);
    defer allocator.free(outline);

    // Read the outline and populate the shape
    try readTrueTypeOutline(shape, outline, scale);
}
