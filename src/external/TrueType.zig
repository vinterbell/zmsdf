const builtin = @import("builtin");
const native_endian = builtin.cpu.arch.endian();

const std = @import("std");
const readInt = std.mem.readInt;
const Allocator = std.mem.Allocator;
const assert = std.debug.assert;
const ArrayList = std.ArrayListUnmanaged;

const TrueType = @This();
const debug_todo = true;

table_offsets: [@typeInfo(TableId).@"enum".fields.len]u32,
ttf_bytes: []const u8,
index_map: u32,
units_per_em: u16,
index_to_loc_format: u16,
glyphs_len: u32,

pub const GlyphIndex = enum(u32) { _ };

pub const TableId = enum {
    cmap,
    loca,
    head,
    glyf,
    hhea,
    hmtx,
    kern,
    GPOS,
    maxp,

    fn asInt(id: TableId) u32 {
        const array4: [4]u8 = @tagName(id).*;
        return @bitCast(array4);
    }
};

const PlatformId = enum(u16) {
    unicode = 0,
    mac = 1,
    iso = 2,
    microsoft = 3,
};

const MicrosoftEncodingId = enum(u16) {
    symbol = 0,
    unicode_bmp = 1,
    shiftjis = 2,
    unicode_full = 10,
};

pub fn load(bytes: []const u8) !TrueType {
    // Find tables.
    var table_offsets = [1]u32{0} ** @typeInfo(TableId).@"enum".fields.len;
    const tables_len = readInt(u16, bytes[4..][0..2], .big);
    for (0..tables_len) |i| {
        const loc = 12 + 16 * i;
        const id: TableId = switch (readInt(u32, bytes[loc..][0..4], native_endian)) {
            TableId.cmap.asInt() => .cmap,
            TableId.loca.asInt() => .loca,
            TableId.head.asInt() => .head,
            TableId.glyf.asInt() => .glyf,
            TableId.hhea.asInt() => .hhea,
            TableId.hmtx.asInt() => .hmtx,
            TableId.kern.asInt() => .kern,
            TableId.GPOS.asInt() => .GPOS,
            TableId.maxp.asInt() => .maxp,
            else => continue,
        };
        table_offsets[@intFromEnum(id)] = readInt(u32, bytes[loc + 8 ..][0..4], .big);
    }

    if (table_offsets[@intFromEnum(TableId.cmap)] == 0) return error.MissingRequiredTable;
    if (table_offsets[@intFromEnum(TableId.loca)] == 0) return error.MissingRequiredTable;
    if (table_offsets[@intFromEnum(TableId.head)] == 0) return error.MissingRequiredTable;
    if (table_offsets[@intFromEnum(TableId.glyf)] == 0) return error.MissingRequiredTable;
    if (table_offsets[@intFromEnum(TableId.hhea)] == 0) return error.MissingRequiredTable;
    if (table_offsets[@intFromEnum(TableId.hmtx)] == 0) return error.MissingRequiredTable;

    const maxp = table_offsets[@intFromEnum(TableId.maxp)];
    const glyphs_len = if (maxp == 0) 0xffff else readInt(u16, bytes[maxp + 4 ..][0..2], .big);

    const cmap = table_offsets[@intFromEnum(TableId.cmap)];
    const cmap_tables_len = readInt(u16, bytes[cmap + 2 ..][0..2], .big);
    const index_map = im: {
        var i = cmap_tables_len;
        while (true) {
            i -= 1;
            if (i == 0) return error.IndexMapMissing;
            const encoding_record = cmap + 4 + 8 * i;
            const platform_id = readInt(u16, bytes[encoding_record..][0..2], .big);
            switch (platform_id) {
                @intFromEnum(PlatformId.microsoft) => switch (readInt(u16, bytes[encoding_record + 2 ..][0..2], .big)) {
                    @intFromEnum(MicrosoftEncodingId.unicode_bmp),
                    @intFromEnum(MicrosoftEncodingId.unicode_full),
                    => {
                        break :im cmap + readInt(u32, bytes[encoding_record + 4 ..][0..4], .big);
                    },
                    else => continue,
                },
                @intFromEnum(PlatformId.unicode) => {
                    break :im cmap + readInt(u32, bytes[encoding_record + 4 ..][0..4], .big);
                },
                else => continue,
            }
        }
    };

    const head = table_offsets[@intFromEnum(TableId.head)];
    const units_per_em = readInt(u16, bytes[head + 18 ..][0..2], .big);
    const index_to_loc_format = readInt(u16, bytes[head + 50 ..][0..2], .big);

    return .{
        .table_offsets = table_offsets,
        .ttf_bytes = bytes,
        .index_map = index_map,
        .units_per_em = units_per_em,
        .index_to_loc_format = index_to_loc_format,
        .glyphs_len = glyphs_len,
    };
}

pub fn codepointGlyphIndex(tt: *const TrueType, codepoint: u21) ?GlyphIndex {
    const bytes = tt.ttf_bytes;
    const index_map = tt.index_map;
    const format = readInt(u16, bytes[index_map..][0..2], .big);
    switch (format) {
        0 => {
            const n = readInt(u16, bytes[index_map + 2 ..][0..2], .big);
            if (codepoint < n - 6)
                return @enumFromInt(bytes[index_map + 6 + codepoint]);

            return null;
        },
        2 => {
            if (debug_todo) @panic("TODO implement high-byte mapping for japanese/chinese/korean");
            return null;
        },
        4 => {
            const seg_count = readInt(u16, bytes[index_map + 6 ..][0..2], .big) >> 1;
            var search_range = readInt(u16, bytes[index_map + 8 ..][0..2], .big) >> 1;
            var entry_selector = readInt(u16, bytes[index_map + 10 ..][0..2], .big);
            const range_shift = readInt(u16, bytes[index_map + 12 ..][0..2], .big) >> 1;

            // Do a binary search of the segments.
            const end_count = index_map + 14;
            var search = end_count;

            if (codepoint > 0xffff)
                return null;

            // They lie from end_count .. end_count + seg_count but search_range
            // is the nearest power of two.
            if (codepoint >= readInt(u16, bytes[search + range_shift * 2 ..][0..2], .big))
                search += range_shift * 2;

            // Now decrement to bias correctly to find smallest.
            search -= 2;
            while (entry_selector > 0) {
                search_range >>= 1;
                const end = readInt(u16, bytes[search + search_range * 2 ..][0..2], .big);
                if (codepoint > end)
                    search += search_range * 2;
                entry_selector -= 1;
            }
            search += 2;

            const item: u16 = @intCast((search - end_count) >> 1);

            const start = readInt(u16, bytes[index_map + 14 + seg_count * 2 + 2 + 2 * item ..][0..2], .big);
            const last = readInt(u16, bytes[end_count + 2 * item ..][0..2], .big);
            if (codepoint < start or codepoint > last)
                return null;

            const offset = readInt(u16, bytes[index_map + 14 + seg_count * 6 + 2 + 2 * item ..][0..2], .big);
            if (offset == 0)
                return @enumFromInt(@as(i32, codepoint) + readInt(i16, bytes[index_map + 14 + seg_count * 4 + 2 + 2 * item ..][0..2], .big));

            return @enumFromInt(readInt(u16, bytes[offset + (codepoint - start) * 2 + index_map + 14 + seg_count * 6 + 2 + 2 * item ..][0..2], .big));
        },
        6 => {
            const first = readInt(u16, bytes[index_map + 6 ..][0..2], .big);
            const count = readInt(u16, bytes[index_map + 8 ..][0..2], .big);
            if (codepoint >= first and codepoint < first + count)
                return @enumFromInt(readInt(u16, bytes[index_map + 10 + (codepoint - first) * 2 ..][0..2], .big));

            return null;
        },
        12, 13 => {
            const ngroups = readInt(u32, bytes[index_map + 12 ..][0..4], .big);
            var low: u32 = 0;
            var high: u32 = ngroups;
            // Binary search the right group.
            while (low < high) {
                const mid = low + ((high - low) >> 1); // rounds down, so low <= mid < high
                const off = index_map + 16 + mid * 12;
                const start_char = readInt(u32, bytes[off..][0..4], .big);
                const end_char = readInt(u32, bytes[off + 4 ..][0..4], .big);
                if (codepoint < start_char) {
                    high = mid;
                } else if (codepoint > end_char) {
                    low = mid + 1;
                } else {
                    const start_glyph = readInt(u32, bytes[off + 8 ..][0..4], .big);
                    return @enumFromInt(start_glyph + if (format == 12) codepoint - start_char else 0);
                }
            }
            return null;
        },
        else => {
            if (debug_todo) @panic("TODO implement glyphIndex for more formats");
            return null;
        },
    }
}

pub const GlyphBitmap = struct {
    width: u16,
    height: u16,
    /// Offset in pixel space from the glyph origin to the left of the bitmap.
    off_x: i16,
    /// Offset in pixel space from the glyph origin to the top of the bitmap.
    off_y: i16,

    pub const empty: GlyphBitmap = .{
        .width = 0,
        .height = 0,
        .off_x = 0,
        .off_y = 0,
    };
};

pub const GlyphBitmapError = error{ OutOfMemory, GlyphNotFound };

/// Caller owns returned memory.
pub fn glyphBitmap(
    tt: *const TrueType,
    gpa: Allocator,
    /// Appended to the list.
    /// Stored left-to-right, top-to-bottom. 8 bits per pixel. 0 is
    /// transparent, 255 is opaque.
    pixels: *std.ArrayListUnmanaged(u8),
    glyph: GlyphIndex,
    scale_x: f32,
    scale_y: f32,
) GlyphBitmapError!GlyphBitmap {
    return glyphBitmapSubpixel(tt, gpa, pixels, glyph, scale_x, scale_y, 0, 0);
}

const Bitmap = struct {
    w: u32,
    h: u32,
    stride: u32,
    pixels: []u8,
};

/// Caller owns returned memory.
pub fn glyphBitmapSubpixel(
    tt: *const TrueType,
    gpa: Allocator,
    /// Appended to the list.
    /// Stored left-to-right, top-to-bottom. 8 bits per pixel. 0 is
    /// transparent, 255 is opaque.
    pixels: *std.ArrayListUnmanaged(u8),
    glyph: GlyphIndex,
    scale_x: f32,
    scale_y: f32,
    shift_x: f32,
    shift_y: f32,
) GlyphBitmapError!GlyphBitmap {
    const vertices = try glyphShape(tt, gpa, glyph);
    defer gpa.free(vertices);

    assert(scale_x != 0);
    assert(scale_y != 0);

    const box = glyphBitmapBoxSubpixel(tt, glyph, scale_x, scale_y, shift_x, shift_y);

    const w: u32 = @intCast(box.x1 - box.x0);
    const h: u32 = @intCast(box.y1 - box.y0);

    if (w == 0 or h == 0) return .empty;

    var gbm: Bitmap = .{
        .w = w,
        .h = h,
        .stride = w,
        .pixels = try pixels.addManyAsSlice(gpa, w * h),
    };
    errdefer pixels.shrinkRetainingCapacity(pixels.items.len - gbm.pixels.len);

    try rasterize(gpa, &gbm, 0.35, vertices, scale_x, scale_y, shift_x, shift_y, box.x0, box.y0, true);

    return .{
        .width = @intCast(gbm.w),
        .height = @intCast(gbm.h),
        .off_x = @intCast(box.x0),
        .off_y = @intCast(box.y0),
    };
}

pub fn scaleForPixelHeight(tt: *const TrueType, height: f32) f32 {
    const vm = tt.verticalMetrics();
    const fheight: f32 = @floatFromInt(vm.ascent - vm.descent);
    return height / fheight;
}

pub const VerticalMetrics = struct {
    /// The coordinate above the baseline the font extends.
    ascent: i16,
    /// The coordinate below the baseline the font extends (typically negative).
    descent: i16,
    /// The spacing between one row's descent and the next row's ascent.
    line_gap: i16,
};

/// A typical expression for advancing the vertical position is
/// `ascent - descent + lineGap`. These are expressed in unscaled coordinates,
/// which are typically then multiplied by the scale factor for a given font size.
pub fn verticalMetrics(tt: *const TrueType) VerticalMetrics {
    const bytes = tt.ttf_bytes;
    const hhea = tt.table_offsets[@intFromEnum(TableId.hhea)];
    return .{
        .ascent = readInt(i16, bytes[hhea + 4 ..][0..2], .big),
        .descent = readInt(i16, bytes[hhea + 6 ..][0..2], .big),
        .line_gap = readInt(i16, bytes[hhea + 8 ..][0..2], .big),
    };
}

pub const HMetrics = struct {
    /// The offset from the current horizontal position to the next horizontal
    /// position in unscaled coordinates.
    advance_width: i16,
    /// The offset from the current horizontal position to the left edge of the
    /// character in unscaled coordinates.
    left_side_bearing: i16,
};

pub fn glyphHMetrics(tt: *const TrueType, glyph: GlyphIndex) HMetrics {
    const glyph_index = @intFromEnum(glyph);
    const bytes = tt.ttf_bytes;
    const hhea = tt.table_offsets[@intFromEnum(TableId.hhea)];
    const hmtx = tt.table_offsets[@intFromEnum(TableId.hmtx)];
    const n_long_h_metrics = readInt(u16, bytes[hhea + 34 ..][0..2], .big);
    if (glyph_index < n_long_h_metrics) return .{
        .advance_width = readInt(i16, bytes[hmtx + 4 * glyph_index ..][0..2], .big),
        .left_side_bearing = readInt(i16, bytes[hmtx + 4 * glyph_index + 2 ..][0..2], .big),
    };
    return .{
        .advance_width = readInt(i16, bytes[hmtx + 4 * (n_long_h_metrics - 1) ..][0..2], .big),
        .left_side_bearing = readInt(i16, bytes[hmtx + 4 * n_long_h_metrics + 2 * (glyph_index - n_long_h_metrics) ..][0..2], .big),
    };
}

/// An additional amount to advance the horizontal coordinate between the two
/// provided glyphs.
pub fn glyphKernAdvance(tt: *const TrueType, a: GlyphIndex, b: GlyphIndex) i16 {
    const gpos = tt.table_offsets[@intFromEnum(TableId.GPOS)];
    if (gpos > 0) return glyphKernAdvanceGpos(tt, a, b);
    const kern = tt.table_offsets[@intFromEnum(TableId.kern)];
    if (kern > 0) return glyphKernAdvanceKern(tt, a, b);
    return 0;
}

fn glyphKernAdvanceGpos(tt: *const TrueType, a: GlyphIndex, b: GlyphIndex) i16 {
    const bytes = tt.ttf_bytes;
    const gpos = tt.table_offsets[@intFromEnum(TableId.GPOS)];
    assert(gpos > 0);

    if (readInt(u16, bytes[gpos + 0 ..][0..2], .big) != 1) return 0; // Major version 1
    if (readInt(u16, bytes[gpos + 2 ..][0..2], .big) != 0) return 0; // Minor version 0

    const lookup_list_offset: u16 = readInt(u16, bytes[gpos + 8 ..][0..2], .big);
    const lookup_list = gpos + lookup_list_offset;
    const lookup_count: u16 = readInt(u16, bytes[lookup_list..][0..2], .big);

    for (0..lookup_count) |i| {
        const lookup_offset = readInt(u16, bytes[lookup_list + 2 + 2 * i ..][0..2], .big);
        const lookup_table = lookup_list + lookup_offset;

        const lookup_type = readInt(u16, bytes[lookup_table..][0..2], .big);
        const sub_table_count = readInt(u16, bytes[lookup_table + 4 ..][0..2], .big);
        const sub_table_offsets = lookup_table + 6;
        if (lookup_type != 2) // Pair Adjustment Positioning Subtable
            continue;

        for (0..sub_table_count) |sti| {
            const subtable_offset = readInt(u16, bytes[sub_table_offsets + 2 * sti ..][0..2], .big);
            const table = lookup_table + subtable_offset;
            const pos_format = readInt(u16, bytes[table..][0..2], .big);
            const coverage_offset = readInt(u16, bytes[table + 2 ..][0..2], .big);
            const coverage_index = coverageIndex(bytes, table + coverage_offset, a) orelse continue;

            switch (pos_format) {
                1 => {
                    const value_format_1 = readInt(u16, bytes[table + 4 ..][0..2], .big);
                    const value_format_2 = readInt(u16, bytes[table + 6 ..][0..2], .big);
                    if (value_format_1 == 4 and value_format_2 == 0) {
                        const value_record_pair_size_in_bytes = 2;
                        const pair_set_count = readInt(u16, bytes[table + 8 ..][0..2], .big);
                        const pair_pos_offset = readInt(u16, bytes[table + 10 + 2 * coverage_index ..][0..2], .big);
                        const pair_value_table = table + pair_pos_offset;
                        const pair_value_count = readInt(u16, bytes[pair_value_table..][0..2], .big);
                        const pair_value_array = pair_value_table + 2;

                        if (coverage_index >= pair_set_count) return 0;

                        const needle = @intFromEnum(b);
                        var r: u32 = pair_value_count - 1;
                        var l: u32 = 0;

                        // Binary search.
                        while (l <= r) {
                            const m = (l + r) >> 1;
                            const pair_value = pair_value_array + (2 + value_record_pair_size_in_bytes) * m;
                            const second_glyph = readInt(u16, bytes[pair_value..][0..2], .big);
                            const straw = second_glyph;
                            if (needle < straw) {
                                if (m == 0) break;
                                r = m - 1;
                            } else if (needle > straw) {
                                l = m + 1;
                            } else {
                                return readInt(i16, bytes[pair_value + 2 ..][0..2], .big);
                            }
                        }
                    } else {
                        if (debug_todo) @panic("TODO implement more glyphKernAdvanceGpos");
                        return 0;
                    }
                },
                2 => {
                    const value_format_1 = readInt(u16, bytes[table + 4 ..][0..2], .big);
                    const value_format_2 = readInt(u16, bytes[table + 6 ..][0..2], .big);
                    if (value_format_1 == 4 and value_format_2 == 0) {
                        const class_def10_offset = readInt(u16, bytes[table + 8 ..][0..2], .big);
                        const class_def20_offset = readInt(u16, bytes[table + 10 ..][0..2], .big);
                        const glyph1class = glyphClass(bytes, table + class_def10_offset, a);
                        const glyph2class = glyphClass(bytes, table + class_def20_offset, b);

                        const class1_count = readInt(u16, bytes[table + 12 ..][0..2], .big);
                        const class2_count = readInt(u16, bytes[table + 14 ..][0..2], .big);

                        if (glyph1class >= class1_count) return 0; // malformed
                        if (glyph2class >= class2_count) return 0; // malformed

                        const class1_records = table + 16;
                        const class2_records = class1_records + 2 * (glyph1class * class2_count);
                        return readInt(i16, bytes[class2_records + 2 * glyph2class ..][0..2], .big);
                    } else {
                        if (debug_todo) @panic("TODO implement more glyphKernAdvanceGpos");
                        return 0;
                    }
                },
                else => {
                    if (debug_todo) @panic("TODO implement more glyphKernAdvanceGpos");
                    return 0;
                },
            }
        }
    }

    return 0;
}

fn glyphKernAdvanceKern(tt: *const TrueType, a: GlyphIndex, b: GlyphIndex) i16 {
    const bytes = tt.ttf_bytes;
    const kern = tt.table_offsets[@intFromEnum(TableId.kern)];
    assert(kern > 0);
    // we only look at the first table. it must be 'horizontal' and format 0.
    if (readInt(u16, bytes[kern + 2 ..][0..2], .big) < 1) // number of tables, need at least 1
        return 0;
    if (readInt(u16, bytes[kern + 8 ..][0..2], .big) != 1) // horizontal flag must be set in format
        return 0;

    var l: u32 = 0;
    var r: u32 = readInt(u16, bytes[kern + 10 ..][0..2], .big) - 1;
    const needle: u32 = @intFromEnum(a) << 16 | @intFromEnum(b);
    while (l <= r) {
        const m: u32 = (l + r) >> 1;
        const straw: u32 = readInt(u32, bytes[kern + 18 + (m * 6) ..][0..4], .big); // note: unaligned read
        if (needle < straw) {
            r = m - 1;
        } else if (needle > straw) {
            l = m + 1;
        } else {
            return readInt(i16, bytes[kern + 22 + (m * 6) ..][0..2], .big);
        }
    }
    return 0;
}

pub const Vertex = struct {
    x: i16,
    y: i16,
    cx: i16,
    cy: i16,
    cx1: i16,
    cy1: i16,
    type: Type,

    pub const Type = enum(u8) {
        vmove = 1,
        vline = 2,
        vcurve = 3,
        vcubic = 4,
        _,
    };

    fn set(v: *Vertex, ty: Type, x: i32, y: i32, cx: i32, cy: i32) void {
        v.type = ty;
        v.x = @intCast(x);
        v.y = @intCast(y);
        v.cx = @intCast(cx);
        v.cy = @intCast(cy);
    }
};

pub fn glyphShape(tt: *const TrueType, gpa: Allocator, glyph: GlyphIndex) GlyphBitmapError![]Vertex {
    const bytes = tt.ttf_bytes;
    const g = try glyfOffset(tt, glyph);
    var vertices: ArrayList(Vertex) = .empty;
    defer vertices.deinit(gpa);
    const n_contours_signed = readInt(i16, bytes[g..][0..2], .big);

    if (n_contours_signed > 0) {
        const n_contours: u16 = @intCast(n_contours_signed);
        const contours_end_pts: u32 = g + 10;
        const ins: i32 = readInt(u16, bytes[g + 10 + n_contours * 2 ..][0..2], .big);
        var points: u32 = @intCast(g + 10 + @as(i64, n_contours) * 2 + 2 + ins);

        const n: u32 = 1 + readInt(u16, bytes[contours_end_pts + n_contours * 2 - 2 ..][0..2], .big);

        // A loose bound on how many vertices we might need.
        const m: u32 = n + 2 * n_contours;
        try vertices.resize(gpa, m);

        var next_move: i32 = 0;
        var flagcount: u8 = 0;

        // in first pass, we load uninterpreted data into the allocated array
        // above, shifted to the end of the array so we won't overwrite it when
        // we create our final data starting from the front

        // Starting offset for uninterpreted data, regardless of how m ends up being calculated.
        const off: u32 = m - n;

        // first load flags
        {
            var flags: u8 = 0;
            for (0..n) |i| {
                if (flagcount == 0) {
                    flags = bytes[points];
                    points += 1;
                    if ((flags & 8) != 0) {
                        flagcount = bytes[points];
                        points += 1;
                    }
                } else {
                    flagcount -= 1;
                }
                vertices.items[off + i].type = @enumFromInt(flags);
            }
        }

        // now load x coordinates
        var x: i32 = 0;
        for (0..n) |i| {
            const flags = @intFromEnum(vertices.items[off + i].type);
            if ((flags & 2) != 0) {
                const dx: i16 = bytes[points];
                points += 1;
                x += if ((flags & 16) != 0) dx else -dx;
            } else {
                if ((flags & 16) == 0) {
                    x += readInt(i16, bytes[points..][0..2], .big);
                    points += 2;
                }
            }
            vertices.items[off + i].x = @intCast(x);
        }

        // now load y coordinates
        var y: i32 = 0;
        for (0..n) |i| {
            const flags = @intFromEnum(vertices.items[off + i].type);
            if ((flags & 4) != 0) {
                const dy: i16 = bytes[points];
                points += 1;
                y += if ((flags & 32) != 0) dy else -dy;
            } else {
                if ((flags & 32) == 0) {
                    y += readInt(i16, bytes[points..][0..2], .big);
                    points += 2;
                }
            }
            vertices.items[off + i].y = @intCast(y);
        }

        // now convert them to our format
        var num_vertices: u32 = 0;
        var sx: i32 = 0;
        var sy: i32 = 0;
        var cx: i32 = 0;
        var cy: i32 = 0;
        var scx: i32 = 0;
        var scy: i32 = 0;
        var i: u32 = 0;
        var j: u32 = 0;
        var start_off: bool = false;
        var was_off: bool = false;
        while (i < n) : (i += 1) {
            const flags = @intFromEnum(vertices.items[off + i].type);
            x = @intCast(vertices.items[off + i].x);
            y = @intCast(vertices.items[off + i].y);

            if (next_move == i) {
                if (i != 0)
                    num_vertices = closeShape(vertices.items, num_vertices, was_off, start_off, sx, sy, scx, scy, cx, cy);

                // now start the new one
                start_off = (flags & 1) == 0;
                if (start_off) {
                    // if we start off with an off-curve point, then when we need to find a point on the curve
                    // where we can start, and we need to save some state for when we wraparound.
                    scx = x;
                    scy = y;
                    if ((@intFromEnum(vertices.items[off + i + 1].type) & 1) == 0) {
                        // next point is also a curve point, so interpolate an on-point curve
                        sx = (x + vertices.items[off + i + 1].x) >> 1;
                        sy = (y + vertices.items[off + i + 1].y) >> 1;
                    } else {
                        // otherwise just use the next point as our start point
                        sx = vertices.items[off + i + 1].x;
                        sy = vertices.items[off + i + 1].y;
                        i += 1; // we're using point i+1 as the starting point, so skip it
                    }
                } else {
                    sx = x;
                    sy = y;
                }
                vertices.items[num_vertices].set(.vmove, sx, sy, 0, 0);
                num_vertices += 1;
                was_off = false;
                next_move = 1 + readInt(u16, bytes[contours_end_pts + j * 2 ..][0..2], .big);
                j += 1;
            } else {
                if ((flags & 1) == 0) { // if it's a curve
                    if (was_off) {
                        // two off-curve control points in a row means interpolate an on-curve midpoint
                        vertices.items[num_vertices].set(.vcurve, (cx + x) >> 1, (cy + y) >> 1, cx, cy);
                        num_vertices += 1;
                    }
                    cx = x;
                    cy = y;
                    was_off = true;
                } else {
                    if (was_off)
                        vertices.items[num_vertices].set(.vcurve, x, y, cx, cy)
                    else
                        vertices.items[num_vertices].set(.vline, x, y, 0, 0);
                    num_vertices += 1;
                    was_off = false;
                }
            }
        }
        num_vertices = closeShape(vertices.items, num_vertices, was_off, start_off, sx, sy, scx, scy, cx, cy);
        vertices.shrinkRetainingCapacity(num_vertices);
    } else if (n_contours_signed < 0) {
        // Compound shapes.
        var more = true;
        var comp = g + 10;
        while (more) {
            var mtx: [6]f32 = .{ 1, 0, 0, 1, 0, 0 };

            const flags = readCursor(u16, bytes, &comp);
            const gidx: GlyphIndex = @enumFromInt(readCursor(u16, bytes, &comp));

            if ((flags & 2) != 0) { // XY values
                if ((flags & 1) != 0) { // shorts
                    mtx[4] = @floatFromInt(readCursor(i16, bytes, &comp));
                    mtx[5] = @floatFromInt(readCursor(i16, bytes, &comp));
                } else {
                    mtx[4] = @floatFromInt(readCursor(i8, bytes, &comp));
                    mtx[5] = @floatFromInt(readCursor(i8, bytes, &comp));
                }
            } else {
                if (debug_todo) @panic("TODO handle matching point");
            }
            if ((flags & (1 << 3)) != 0) { // WE_HAVE_A_SCALE
                mtx[0] = @as(f32, @floatFromInt(readCursor(i16, bytes, &comp))) / 16384.0;
                mtx[1] = 0;
                mtx[2] = 0;
                mtx[3] = mtx[0];
            } else if ((flags & (1 << 6)) != 0) { // WE_HAVE_AN_X_AND_YSCALE
                mtx[0] = @as(f32, @floatFromInt(readCursor(i16, bytes, &comp))) / 16384.0;
                mtx[1] = 0;
                mtx[2] = 0;
                mtx[3] = @as(f32, @floatFromInt(readCursor(i16, bytes, &comp))) / 16384.0;
            } else if ((flags & (1 << 7)) != 0) { // WE_HAVE_A_TWO_BY_TWO
                mtx[0] = @as(f32, @floatFromInt(readCursor(i16, bytes, &comp))) / 16384.0;
                mtx[1] = @as(f32, @floatFromInt(readCursor(i16, bytes, &comp))) / 16384.0;
                mtx[2] = @as(f32, @floatFromInt(readCursor(i16, bytes, &comp))) / 16384.0;
                mtx[3] = @as(f32, @floatFromInt(readCursor(i16, bytes, &comp))) / 16384.0;
            }

            // Find transformation scales.
            const m: f32 = @sqrt(mtx[0] * mtx[0] + mtx[1] * mtx[1]);
            const n: f32 = @sqrt(mtx[2] * mtx[2] + mtx[3] * mtx[3]);

            // Get indexed glyph.
            const comp_verts = try glyphShape(tt, gpa, gidx);
            defer gpa.free(comp_verts);
            if (comp_verts.len > 0) {
                // Transform vertices.
                for (comp_verts) |*v| {
                    {
                        const x: f32 = @floatFromInt(v.x);
                        const y: f32 = @floatFromInt(v.y);
                        v.x = @intFromFloat(m * (mtx[0] * x + mtx[2] * y + mtx[4]));
                        v.y = @intFromFloat(n * (mtx[1] * x + mtx[3] * y + mtx[5]));
                    }
                    {
                        const x: f32 = @floatFromInt(v.cx);
                        const y: f32 = @floatFromInt(v.cy);
                        v.cx = @intFromFloat(m * (mtx[0] * x + mtx[2] * y + mtx[4]));
                        v.cy = @intFromFloat(n * (mtx[1] * x + mtx[3] * y + mtx[5]));
                    }
                }
                try vertices.appendSlice(gpa, comp_verts);
            }
            more = (flags & (1 << 5)) != 0;
        }
    }
    return vertices.toOwnedSlice(gpa);
}

fn glyfOffset(tt: *const TrueType, glyph: GlyphIndex) error{GlyphNotFound}!u32 {
    const bytes = tt.ttf_bytes;
    const glyph_index: u32 = @intFromEnum(glyph);

    assert(glyph_index < tt.glyphs_len);
    assert(tt.index_to_loc_format < 2);

    const glyf = tt.table_offsets[@intFromEnum(TableId.glyf)];
    const loca = tt.table_offsets[@intFromEnum(TableId.loca)];
    const g1, const g2 = if (tt.index_to_loc_format == 0) .{
        glyf + readInt(u16, bytes[loca + glyph_index * 2 ..][0..2], .big) * 2,
        glyf + readInt(u16, bytes[loca + glyph_index * 2 + 2 ..][0..2], .big) * 2,
    } else .{
        glyf + readInt(u32, bytes[loca + glyph_index * 4 ..][0..4], .big),
        glyf + readInt(u32, bytes[loca + glyph_index * 4 + 4 ..][0..4], .big),
    };
    if (g1 == g2) return error.GlyphNotFound;
    return g1;
}

const BitmapBox = struct {
    x0: i32,
    y0: i32,
    x1: i32,
    y1: i32,
};

fn glyphBitmapBoxSubpixel(
    tt: *const TrueType,
    glyph: GlyphIndex,
    scale_x: f32,
    scale_y: f32,
    shift_x: f32,
    shift_y: f32,
) BitmapBox {
    const box = glyphBox(tt, glyph) catch |err| switch (err) {
        error.GlyphNotFound => return .{ .x0 = 0, .y0 = 0, .x1 = 0, .y1 = 0 }, // e.g. space character
    };
    return .{
        // move to integral bboxes (treating pixels as little squares, what pixels get touched)?
        .x0 = @intFromFloat(@floor(@as(f32, @floatFromInt(box.x0)) * scale_x + shift_x)),
        .y0 = @intFromFloat(@floor(@as(f32, @floatFromInt(-box.y1)) * scale_y + shift_y)),
        .x1 = @intFromFloat(@ceil(@as(f32, @floatFromInt(box.x1)) * scale_x + shift_x)),
        .y1 = @intFromFloat(@ceil(@as(f32, @floatFromInt(-box.y0)) * scale_y + shift_y)),
    };
}

fn glyphBox(tt: *const TrueType, glyph: GlyphIndex) error{GlyphNotFound}!BitmapBox {
    const bytes = tt.ttf_bytes;
    const g = try glyfOffset(tt, glyph);
    return .{
        .x0 = readInt(i16, bytes[g + 2 ..][0..2], .big),
        .y0 = readInt(i16, bytes[g + 4 ..][0..2], .big),
        .x1 = readInt(i16, bytes[g + 6 ..][0..2], .big),
        .y1 = readInt(i16, bytes[g + 8 ..][0..2], .big),
    };
}

fn rasterize(
    gpa: Allocator,
    result: *Bitmap,
    flatness_in_pixels: f32,
    vertices: []Vertex,
    scale_x: f32,
    scale_y: f32,
    shift_x: f32,
    shift_y: f32,
    off_x: i32,
    off_y: i32,
    invert: bool,
) Allocator.Error!void {
    const scale = @min(scale_x, scale_y);
    var windings = try flattenCurves(gpa, vertices, flatness_in_pixels / scale);
    defer windings.deinit(gpa);
    try rasterizeInner(gpa, result, windings.points, windings.contour_lengths, scale_x, scale_y, shift_x, shift_y, off_x, off_y, invert);
}

const Edge = struct {
    x0: f32,
    y0: f32,
    x1: f32,
    y1: f32,
    invert: bool,

    const Sort = struct {
        fn lessThan(ctx: Sort, a: Edge, b: Edge) bool {
            _ = ctx;
            return a.y0 < b.y0;
        }
    };
};

fn rasterizeInner(
    gpa: Allocator,
    result: *Bitmap,
    pts: []Point,
    wcount: []u32,
    scale_x: f32,
    scale_y: f32,
    shift_x: f32,
    shift_y: f32,
    off_x: i32,
    off_y: i32,
    invert: bool,
) Allocator.Error!void {
    const y_scale_inv: f32 = if (invert) -scale_y else scale_y;

    // now we have to blow out the windings into explicit edge lists
    const edge_alloc_n = n: {
        var n: u32 = 1; // Add an extra one as a sentinel.
        for (wcount) |elem| n += elem;
        break :n n;
    };

    const e = try gpa.alloc(Edge, edge_alloc_n);
    defer gpa.free(e);

    var n: u32 = 0;
    var m: u32 = 0;
    for (wcount) |wcount_elem| {
        const p: []Point = pts[m..];
        m += wcount_elem;
        var j: u32 = wcount_elem - 1;
        var k: u32 = 0;
        while (k < wcount_elem) : ({
            j = k;
            k += 1;
        }) {
            var a = k;
            var b = j;
            // skip the edge if horizontal
            if (p[j].y == p[k].y)
                continue;
            // add edge from j to k to the list
            e[n].invert = false;
            if (if (invert) p[j].y > p[k].y else p[j].y < p[k].y) {
                e[n].invert = true;
                a = j;
                b = k;
            }
            e[n].x0 = p[a].x * scale_x + shift_x;
            e[n].y0 = (p[a].y * y_scale_inv + shift_y);
            e[n].x1 = p[b].x * scale_x + shift_x;
            e[n].y1 = (p[b].y * y_scale_inv + shift_y);
            n += 1;
        }
    }
    // now sort the edges by their highest point (should snap to integer, and then by x)
    std.mem.sortUnstable(Edge, e[0..n], Edge.Sort{}, Edge.Sort.lessThan);

    // now, traverse the scanlines and find the intersections on each scanline, use xor winding rule
    try rasterizeSortedEdges(gpa, result, e[0 .. n + 1], off_x, off_y);
}

const Point = struct {
    x: f32,
    y: f32,
};

const FlattenedCurves = struct {
    points: []Point,
    contour_lengths: []u32,

    const empty: FlattenedCurves = .{
        .points = &.{},
        .contour_lengths = &.{},
    };

    fn deinit(fc: *FlattenedCurves, gpa: Allocator) void {
        gpa.free(fc.points);
        gpa.free(fc.contour_lengths);
        fc.* = undefined;
    }
};

fn flattenCurves(
    gpa: Allocator,
    vertices: []const Vertex,
    objspace_flatness: f32,
) error{OutOfMemory}!FlattenedCurves {
    var points: ArrayList(Point) = .empty;
    defer points.deinit(gpa);
    var contour_lengths: ArrayList(u32) = .empty;
    defer contour_lengths.deinit(gpa);

    const objspace_flatness_squared = objspace_flatness * objspace_flatness;

    var start: u32 = 0;
    var x: f32 = 0;
    var y: f32 = 0;
    for (vertices) |v| {
        sw: switch (v.type) {
            .vmove => {
                if (points.items.len > 0) {
                    try contour_lengths.append(gpa, @intCast(points.items.len - start));
                    start = @intCast(points.items.len);
                }

                continue :sw .vline;
            },
            .vline => {
                x = @floatFromInt(v.x);
                y = @floatFromInt(v.y);
                try points.append(gpa, .{ .x = x, .y = y });
            },
            .vcurve => {
                try tesselateCurve(
                    gpa,
                    &points,
                    x,
                    y,
                    @floatFromInt(v.cx),
                    @floatFromInt(v.cy),
                    @floatFromInt(v.x),
                    @floatFromInt(v.y),
                    objspace_flatness_squared,
                    0,
                );
                x = @floatFromInt(v.x);
                y = @floatFromInt(v.y);
            },
            .vcubic => {
                try tesselateCubic(
                    gpa,
                    &points,
                    x,
                    y,
                    @floatFromInt(v.cx),
                    @floatFromInt(v.cy),
                    @floatFromInt(v.cx1),
                    @floatFromInt(v.cy1),
                    @floatFromInt(v.x),
                    @floatFromInt(v.y),
                    objspace_flatness_squared,
                    0,
                );
                x = @floatFromInt(v.x);
                y = @floatFromInt(v.y);
            },
            _ => continue,
        }
    }
    try contour_lengths.append(gpa, @intCast(points.items.len - start));

    return .{
        .points = try points.toOwnedSlice(gpa),
        .contour_lengths = try contour_lengths.toOwnedSlice(gpa),
    };
}

/// tessellate until threshold p is happy... @TODO warped to compensate for non-linear stretching
fn tesselateCurve(
    gpa: Allocator,
    points: *ArrayList(Point),
    x0: f32,
    y0: f32,
    x1: f32,
    y1: f32,
    x2: f32,
    y2: f32,
    objspace_flatness_squared: f32,
    n: u32,
) Allocator.Error!void {
    // midpoint
    const mx: f32 = (x0 + 2 * x1 + x2) / 4;
    const my: f32 = (y0 + 2 * y1 + y2) / 4;
    // versus directly drawn line
    const dx: f32 = (x0 + x2) / 2 - mx;
    const dy: f32 = (y0 + y2) / 2 - my;
    if (n > 16) // 65536 segments on one curve better be enough!
        return;
    if (dx * dx + dy * dy > objspace_flatness_squared) { // half-pixel error allowed... need to be smaller if AA
        try tesselateCurve(gpa, points, x0, y0, (x0 + x1) / 2.0, (y0 + y1) / 2.0, mx, my, objspace_flatness_squared, n + 1);
        try tesselateCurve(gpa, points, mx, my, (x1 + x2) / 2.0, (y1 + y2) / 2.0, x2, y2, objspace_flatness_squared, n + 1);
    } else {
        try points.append(gpa, .{ .x = x2, .y = y2 });
    }
}

fn tesselateCubic(
    gpa: Allocator,
    points: *ArrayList(Point),
    x0: f32,
    y0: f32,
    x1: f32,
    y1: f32,
    x2: f32,
    y2: f32,
    x3: f32,
    y3: f32,
    objspace_flatness_squared: f32,
    n: u32,
) Allocator.Error!void {
    // According to Dougall Johnson, this "flatness" calculation is just
    // made-up nonsense that seems to work well enough.
    const dx0 = x1 - x0;
    const dy0 = y1 - y0;
    const dx1 = x2 - x1;
    const dy1 = y2 - y1;
    const dx2 = x3 - x2;
    const dy2 = y3 - y2;
    const dx = x3 - x0;
    const dy = y3 - y0;
    const longlen = @sqrt(dx0 * dx0 + dy0 * dy0) + @sqrt(dx1 * dx1 + dy1 * dy1) + @sqrt(dx2 * dx2 + dy2 * dy2);
    const shortlen = @sqrt(dx * dx + dy * dy);
    const flatness_squared = longlen * longlen - shortlen * shortlen;

    if (n > 16) // 65536 segments on one curve better be enough!
        return;

    if (flatness_squared > objspace_flatness_squared) {
        const x01 = (x0 + x1) / 2;
        const y01 = (y0 + y1) / 2;
        const x12 = (x1 + x2) / 2;
        const y12 = (y1 + y2) / 2;
        const x23 = (x2 + x3) / 2;
        const y23 = (y2 + y3) / 2;

        const xa = (x01 + x12) / 2;
        const ya = (y01 + y12) / 2;
        const xb = (x12 + x23) / 2;
        const yb = (y12 + y23) / 2;

        const mx = (xa + xb) / 2;
        const my = (ya + yb) / 2;

        try tesselateCubic(gpa, points, x0, y0, x01, y01, xa, ya, mx, my, objspace_flatness_squared, n + 1);
        try tesselateCubic(gpa, points, mx, my, xb, yb, x23, y23, x3, y3, objspace_flatness_squared, n + 1);
    } else {
        try points.append(gpa, .{ .x = x3, .y = y3 });
    }
}

fn sizedTrapezoidArea(height: f32, top_width: f32, bottom_width: f32) f32 {
    assert(top_width >= 0);
    assert(bottom_width >= 0);
    return (top_width + bottom_width) / 2.0 * height;
}

fn positionTrapezoidArea(height: f32, tx0: f32, tx1: f32, bx0: f32, bx1: f32) f32 {
    return sizedTrapezoidArea(height, tx1 - tx0, bx1 - bx0);
}

fn sizedTriangleArea(height: f32, width: f32) f32 {
    return height * width / 2;
}

const ActiveEdge = struct {
    next: ?*ActiveEdge,
    fx: f32,
    fdx: f32,
    fdy: f32,
    direction: f32,
    sy: f32,
    ey: f32,
};

/// Directly anti-alias rasterize edges without supersampling.
fn rasterizeSortedEdges(
    gpa: Allocator,
    result: *Bitmap,
    edges: []Edge,
    off_x: i32,
    off_y: i32,
) Allocator.Error!void {
    var arena_allocator = std.heap.ArenaAllocator.init(gpa);
    defer arena_allocator.deinit();
    const arena = arena_allocator.allocator();

    var active: ?*ActiveEdge = null;
    var scanline_buffer: [350]f32 = undefined;

    const needed_scanline_len = result.w * 2 + 1;
    assert(scanline_buffer.len >= needed_scanline_len);

    const scanline = scanline_buffer[0..result.w];
    const scanline2 = scanline_buffer[result.w..][0 .. result.w + 1];

    var y: i32 = off_y;
    edges[edges.len - 1].y0 = @floatFromInt((off_y + @as(i32, @intCast(result.h))) + 1);

    var j: u32 = 0;
    var e: u32 = 0;
    while (j < result.h) {
        // find center of pixel for this scanline
        const scan_y_top: f32 = @floatFromInt(y);
        const scan_y_bottom: f32 = @floatFromInt(y + 1);
        var step: *?*ActiveEdge = &active;

        @memset(scanline, 0);
        @memset(scanline2, 0);

        // update all active edges;
        // remove all active edges that terminate before the top of this scanline
        while (step.*) |z| {
            if (z.ey <= scan_y_top) {
                step.* = z.next; // delete from list
                assert(z.direction != 0);
                z.direction = 0;
                arena.destroy(z);
            } else {
                step = &z.next; // advance through list
            }
        }

        // insert all edges that start before the bottom of this scanline
        while (edges[e].y0 <= scan_y_bottom) {
            if (edges[e].y0 != edges[e].y1) {
                const z: *ActiveEdge = try newActive(arena, edges[e], off_x, scan_y_top);
                if (j == 0 and off_y != 0) {
                    z.ey = @max(z.ey, scan_y_top);
                }
                // If we get really unlucky a tiny bit of an edge can be
                // out of bounds.
                assert(z.ey >= scan_y_top);

                // Insert at front.
                z.next = active;
                active = z;
            }
            e += 1;
        }

        if (active) |a| fillActiveEdges(scanline, scanline2, result.w, a, scan_y_top);

        {
            var sum: f32 = 0;
            for (scanline, scanline2[0..result.w], result.pixels[j * result.stride ..][0..result.w]) |s, s2, *p| {
                sum += s2;
                p.* = @intFromFloat(@min(@abs(s + sum) * 255 + 0.5, 255));
            }
        }
        // advance all the edges
        step = &active;
        while (step.*) |z| {
            z.fx += z.fdx; // advance to position for current scanline
            step = &z.next; // advance through list
        }

        y += 1;
        j += 1;
    }
}

fn closeShape(
    vertices: []Vertex,
    vertices_len_start: u32,
    was_off: bool,
    start_off: bool,
    sx: i32,
    sy: i32,
    scx: i32,
    scy: i32,
    cx: i32,
    cy: i32,
) u32 {
    var vertices_len = vertices_len_start;
    if (start_off) {
        if (was_off) {
            vertices[vertices_len].set(.vcurve, (cx + scx) >> 1, (cy + scy) >> 1, cx, cy);
            vertices_len += 1;
        }
        vertices[vertices_len].set(.vcurve, sx, sy, scx, scy);
        vertices_len += 1;
    } else {
        if (was_off) {
            vertices[vertices_len].set(.vcurve, sx, sy, cx, cy);
            vertices_len += 1;
        } else {
            vertices[vertices_len].set(.vline, sx, sy, 0, 0);
            vertices_len += 1;
        }
    }
    return vertices_len;
}

fn readCursor(comptime I: type, bytes: []const u8, cursor: *u32) I {
    const start = cursor.*;
    const result = readInt(I, bytes[start..][0..@sizeOf(I)], .big);
    cursor.* = start + @sizeOf(I);
    return result;
}

fn newActive(arena: Allocator, e: Edge, off_x: i32, start_point: f32) Allocator.Error!*ActiveEdge {
    const z = try arena.create(ActiveEdge);
    const dxdy: f32 = (e.x1 - e.x0) / (e.y1 - e.y0);
    z.* = .{
        .fdx = dxdy,
        .fdy = if (dxdy != 0.0) (1.0 / dxdy) else 0.0,
        .fx = (e.x0 + dxdy * (start_point - e.y0)) - @as(f32, @floatFromInt(off_x)),
        .direction = if (e.invert) 1.0 else -1.0,
        .sy = e.y0,
        .ey = e.y1,
        .next = null,
    };
    return z;
}

fn fillActiveEdges(scanline: []f32, scanline_fill: []f32, len: u32, start_edge: *ActiveEdge, y_top: f32) void {
    const y_bottom: f32 = y_top + 1;
    var opt_e: ?*ActiveEdge = start_edge;
    while (opt_e) |e| : (opt_e = e.next) {
        // brute force every pixel

        // compute intersection points with top & bottom
        assert(e.ey >= y_top);

        if (e.fdx == 0) {
            const x0 = e.fx;
            if (x0 < @as(f32, @floatFromInt(len))) {
                if (x0 >= 0) {
                    handleClippedEdge(scanline, @intFromFloat(x0), e, x0, y_top, x0, y_bottom);
                    handleClippedEdge(scanline_fill, @intFromFloat(x0 + 1), e, x0, y_top, x0, y_bottom);
                } else {
                    handleClippedEdge(scanline_fill, 0, e, x0, y_top, x0, y_bottom);
                }
            }
        } else {
            var x0: f32 = e.fx;
            var dx: f32 = e.fdx;
            var xb: f32 = x0 + dx;
            var dy: f32 = e.fdy;
            assert(e.sy <= y_bottom);
            assert(e.ey >= y_top);

            // Compute endpoints of line segment clipped to this scanline (if the
            // line segment starts on this scanline. x0 is the intersection of the
            // line with y_top, but that may be off the line segment.
            var x_top: f32, var sy0: f32 = if (e.sy > y_top) .{
                x0 + dx * (e.sy - y_top),
                e.sy,
            } else .{
                x0,
                y_top,
            };

            var x_bottom: f32, var sy1: f32 = if (e.ey < y_bottom) .{
                x0 + dx * (e.ey - y_top),
                e.ey,
            } else .{
                xb,
                y_bottom,
            };

            if (x_top >= 0 and x_bottom >= 0 and
                x_top < @as(f32, @floatFromInt(len)) and x_bottom < @as(f32, @floatFromInt(len)))
            {
                // from here on, we don't have to range check x values

                if (@trunc(x_top) == @trunc(x_bottom)) {
                    // simple case, only spans one pixel
                    const x: u32 = @intFromFloat(x_top);
                    const height: f32 = (sy1 - sy0) * e.direction;
                    assert(x < len);
                    scanline[x] += positionTrapezoidArea(height, x_top, @floatFromInt(x + 1), x_bottom, @floatFromInt(x + 1));
                    scanline_fill[x + 1] += height; // everything right of this pixel is filled
                } else {
                    // covers 2+ pixels
                    if (x_top > x_bottom) {
                        // flip scanline vertically; signed area is the same
                        sy0 = y_bottom - (sy0 - y_top);
                        sy1 = y_bottom - (sy1 - y_top);
                        std.mem.swap(f32, &sy0, &sy1);
                        std.mem.swap(f32, &x_bottom, &x_top);
                        dx = -dx;
                        dy = -dy;
                        std.mem.swap(f32, &x0, &xb);
                    }
                    assert(dy >= 0);
                    assert(dx >= 0);

                    const x1: u32 = @intFromFloat(x_top);
                    const x2: u32 = @intFromFloat(x_bottom);
                    const x1p1f: f32 = @floatFromInt(x1 + 1);
                    const x2f: f32 = @floatFromInt(x2);
                    // compute intersection with y axis at x1+1
                    var y_crossing: f32 = y_top + dy * (x1p1f - x0);

                    // compute intersection with y axis at x2
                    var y_final: f32 = y_top + dy * (x2f - x0);

                    //           x1    x_top                            x2    x_bottom
                    //     y_top  +------|-----+------------+------------+--------|---+------------+
                    //            |            |            |            |            |            |
                    //            |            |            |            |            |            |
                    //       sy0  |      Txxxxx|............|............|............|............|
                    // y_crossing |            *xxxxx.......|............|............|............|
                    //            |            |     xxxxx..|............|............|............|
                    //            |            |     /-   xx*xxxx........|............|............|
                    //            |            | dy <       |    xxxxxx..|............|............|
                    //   y_final  |            |     \-     |          xx*xxx.........|............|
                    //       sy1  |            |            |            |   xxxxxB...|............|
                    //            |            |            |            |            |            |
                    //            |            |            |            |            |            |
                    //  y_bottom  +------------+------------+------------+------------+------------+
                    //
                    // goal is to measure the area covered by '.' in each pixel

                    // if x2 is right at the right edge of x1, y_crossing can blow up, github #1057
                    // @TODO: maybe test against sy1 rather than y_bottom?
                    if (y_crossing > y_bottom)
                        y_crossing = y_bottom;

                    const sign: f32 = e.direction;

                    // area of the rectangle covered from sy0..y_crossing
                    var area: f32 = sign * (y_crossing - sy0);

                    // area of the triangle (x_top,sy0), (x1+1,sy0), (x1+1,y_crossing)
                    scanline[x1] += sizedTriangleArea(area, x1p1f - x_top);

                    // check if final y_crossing is blown up; no test case for this
                    if (y_final > y_bottom) {
                        y_final = y_bottom;
                        dy = (y_final - y_crossing) / (x2f - x1p1f); // if denom=0, y_final = y_crossing, so y_final <= y_bottom
                    }

                    // in second pixel, area covered by line segment found in first pixel
                    // is always a rectangle 1 wide * the height of that line segment; this
                    // is exactly what the variable 'area' stores. it also gets a contribution
                    // from the line segment within it. the THIRD pixel will get the first
                    // pixel's rectangle contribution, the second pixel's rectangle contribution,
                    // and its own contribution. the 'own contribution' is the same in every pixel except
                    // the leftmost and rightmost, a trapezoid that slides down in each pixel.
                    // the second pixel's contribution to the third pixel will be the
                    // rectangle 1 wide times the height change in the second pixel, which is dy.

                    const step: f32 = sign * dy * 1; // dy is dy/dx, change in y for every 1 change in x,
                    // which multiplied by 1-pixel-width is how much pixel area changes for each step in x
                    // so the area advances by 'step' every time

                    for (scanline[x1 + 1 .. x2]) |*s| {
                        s.* += area + step / 2; // area of trapezoid is 1*step/2
                        area += step;
                    }
                    assert(@abs(area) <= 1.01); // accumulated error from area += step unless we round step down
                    assert(sy1 > y_final - 0.01);

                    // area covered in the last pixel is the rectangle from all the pixels to the left,
                    // plus the trapezoid filled by the line segment in this pixel all the way to the right edge
                    scanline[x2] += area + sign * positionTrapezoidArea(sy1 - y_final, x2f, x2f + 1.0, x_bottom, x2f + 1.0);

                    // the rest of the line is filled based on the total height of the line segment in this pixel
                    scanline_fill[x2 + 1] += sign * (sy1 - sy0);
                }
            } else {
                // if edge goes outside of box we're drawing, we require
                // clipping logic. since this does not match the intended use
                // of this library, we use a different, very slow brute
                // force implementation
                // note though that this does happen some of the time because
                // x_top and x_bottom can be extrapolated at the top & bottom of
                // the shape and actually lie outside the bounding box
                for (0..len) |x_usize| {
                    const x: u32 = @intCast(x_usize);
                    // cases:
                    //
                    // there can be up to two intersections with the pixel. any intersection
                    // with left or right edges can be handled by splitting into two (or three)
                    // regions. intersections with top & bottom do not necessitate case-wise logic.
                    //
                    // the old way of doing this found the intersections with the left & right edges,
                    // then used some simple logic to produce up to three segments in sorted order
                    // from top-to-bottom. however, this had a problem: if an x edge was epsilon
                    // across the x border, then the corresponding y position might not be distinct
                    // from the other y segment, and it might ignored as an empty segment. to avoid
                    // that, we need to explicitly produce segments based on x positions.

                    // rename variables to clearly-defined pairs
                    const y0: f32 = y_top;
                    const x1: f32 = @floatFromInt(x);
                    const x2: f32 = @floatFromInt(x + 1);
                    const x3: f32 = xb;
                    const y3: f32 = y_bottom;

                    // x = e.x + e.dx * (y-y_top)
                    // (y-y_top) = (x - e.x) / e.dx
                    // y = (x - e.x) / e.dx + y_top
                    const y1: f32 = (x1 - x0) / dx + y_top;
                    const y2: f32 = (x1 + 1 - x0) / dx + y_top;

                    if (x0 < x1 and x3 > x2) { // three segments descending down-right
                        handleClippedEdge(scanline, x, e, x0, y0, x1, y1);
                        handleClippedEdge(scanline, x, e, x1, y1, x2, y2);
                        handleClippedEdge(scanline, x, e, x2, y2, x3, y3);
                    } else if (x3 < x1 and x0 > x2) { // three segments descending down-left
                        handleClippedEdge(scanline, x, e, x0, y0, x2, y2);
                        handleClippedEdge(scanline, x, e, x2, y2, x1, y1);
                        handleClippedEdge(scanline, x, e, x1, y1, x3, y3);
                    } else if (x0 < x1 and x3 > x1) { // two segments across x, down-right
                        handleClippedEdge(scanline, x, e, x0, y0, x1, y1);
                        handleClippedEdge(scanline, x, e, x1, y1, x3, y3);
                    } else if (x3 < x1 and x0 > x1) { // two segments across x, down-left
                        handleClippedEdge(scanline, x, e, x0, y0, x1, y1);
                        handleClippedEdge(scanline, x, e, x1, y1, x3, y3);
                    } else if (x0 < x2 and x3 > x2) { // two segments across x+1, down-right
                        handleClippedEdge(scanline, x, e, x0, y0, x2, y2);
                        handleClippedEdge(scanline, x, e, x2, y2, x3, y3);
                    } else if (x3 < x2 and x0 > x2) { // two segments across x+1, down-left
                        handleClippedEdge(scanline, x, e, x0, y0, x2, y2);
                        handleClippedEdge(scanline, x, e, x2, y2, x3, y3);
                    } else { // one segment
                        handleClippedEdge(scanline, x, e, x0, y0, x3, y3);
                    }
                }
            }
        }
    }
}

/// The edge passed in here does not cross the vertical line at x or the
/// vertical line at x+1 (i.e. it has already been clipped to those).
fn handleClippedEdge(
    scanline: []f32,
    x: u32,
    e: *ActiveEdge,
    x0_start: f32,
    y0_start: f32,
    x1_start: f32,
    y1_start: f32,
) void {
    var x0 = x0_start;
    var y0 = y0_start;
    var x1 = x1_start;
    var y1 = y1_start;
    if (y0 == y1) return;
    assert(y0 < y1);
    assert(e.sy <= e.ey);
    if (y0 > e.ey) return;
    if (y1 < e.sy) return;
    if (y0 < e.sy) {
        x0 += (x1 - x0) * (e.sy - y0) / (y1 - y0);
        y0 = e.sy;
    }
    if (y1 > e.ey) {
        x1 += (x1 - x0) * (e.ey - y1) / (y1 - y0);
        y1 = e.ey;
    }

    const xf: f32 = @floatFromInt(x);

    if (x0 == xf)
        assert(x1 <= xf + 1)
    else if (x0 == xf + 1)
        assert(x1 >= xf)
    else if (x0 <= xf)
        assert(x1 <= xf)
    else if (x0 >= xf + 1)
        assert(x1 >= xf + 1)
    else {
        assert(x1 >= xf);
        assert(x1 <= xf + 1);
    }

    if (x0 <= xf and x1 <= xf) {
        scanline[x] += e.direction * (y1 - y0);
    } else if (x0 >= xf + 1 and x1 >= xf + 1) {
        // Do nothing.
    } else {
        assert(x0 >= xf);
        assert(x0 <= xf + 1);
        assert(x1 >= xf);
        assert(x1 <= xf + 1);
        // coverage = 1 - average x position
        scanline[x] += e.direction * (y1 - y0) * (1 - ((x0 - xf) + (x1 - xf)) / 2);
    }
}

fn coverageIndex(bytes: []const u8, coverage_table: u32, glyph: GlyphIndex) ?u32 {
    const coverage_format = readInt(u16, bytes[coverage_table..][0..2], .big);
    switch (coverage_format) {
        1 => {
            const glyph_count = readInt(u16, bytes[coverage_table + 2 ..][0..2], .big);

            // Binary search.
            var l: u32 = 0;
            var r: u32 = glyph_count - 1;
            const needle = @intFromEnum(glyph);
            while (l <= r) {
                const glyph_array = coverage_table + 4;
                const m = (l + r) >> 1;
                const glyph_id = readInt(u16, bytes[glyph_array + 2 * m ..][0..2], .big);
                const straw = glyph_id;
                if (needle < straw) {
                    if (m == 0) break;
                    r = m - 1;
                } else if (needle > straw) {
                    l = m + 1;
                } else {
                    return m;
                }
            }
        },
        2 => {
            const range_count = readInt(u16, bytes[coverage_table + 2 ..][0..2], .big);
            const range_array = coverage_table + 4;

            // Binary search.
            var l: u32 = 0;
            var r: u32 = range_count - 1;
            const needle = @intFromEnum(glyph);
            while (l <= r) {
                const m = (l + r) >> 1;
                const range_record = range_array + 6 * m;
                const straw_start = readInt(u16, bytes[range_record..][0..2], .big);
                const straw_end = readInt(u16, bytes[range_record + 2 ..][0..2], .big);
                if (needle < straw_start) {
                    if (m == 0) break;
                    r = m - 1;
                } else if (needle > straw_end) {
                    l = m + 1;
                } else {
                    const start_coverage_index = readInt(u16, bytes[range_record + 4 ..][0..2], .big);
                    return start_coverage_index + needle - straw_start;
                }
            }
        },
        else => {},
    }
    return null;
}

fn glyphClass(bytes: []const u8, class_def_table: u32, glyph: GlyphIndex) u32 {
    const glyph_int = @intFromEnum(glyph);
    const class_def_format = readInt(u16, bytes[class_def_table..][0..2], .big);
    switch (class_def_format) {
        1 => {
            const start_glyph_id = readInt(u16, bytes[class_def_table + 2 ..][0..2], .big);
            const glyph_count = readInt(u16, bytes[class_def_table + 4 ..][0..2], .big);
            const class_def1_value_array = class_def_table + 6;

            if (glyph_int >= start_glyph_id and glyph_int < start_glyph_id + glyph_count)
                return readInt(u16, bytes[class_def1_value_array + 2 * (glyph_int - start_glyph_id) ..][0..2], .big);
        },
        2 => {
            const class_range_count = readInt(u16, bytes[class_def_table + 2 ..][0..2], .big);
            const class_range_records = class_def_table + 4;

            // Binary search.
            var l: u32 = 0;
            var r: u32 = class_range_count - 1;
            while (l <= r) {
                const m = (l + r) >> 1;
                const class_range_record = class_range_records + 6 * m;
                const straw_start = readInt(u16, bytes[class_range_record..][0..2], .big);
                const straw_end = readInt(u16, bytes[class_range_record + 2 ..][0..2], .big);
                if (glyph_int < straw_start) {
                    if (m == 0) break;
                    r = m - 1;
                } else if (glyph_int > straw_end) {
                    l = m + 1;
                } else {
                    return readInt(u16, bytes[class_range_record + 4 ..][0..2], .big);
                }
            }
        },
        else => return std.math.maxInt(u32), // Unsupported definition type, return an error.
    }

    // "All glyphs not assigned to a class fall into class 0". (OpenType spec)
    return 0;
}
