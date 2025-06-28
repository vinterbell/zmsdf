const std = @import("std");
const zmsdf = @import("zmsdf");
const builtin = @import("builtin");

pub fn main() !void {
    var dbg_allocator: std.heap.DebugAllocator(.{}) = .init;
    defer _ = dbg_allocator.deinit();
    const allocator = dbg_allocator.allocator();

    var args = try std.process.argsWithAllocator(allocator);
    defer args.deinit();

    _ = args.next(); // skip the first argument (the program name)

    const size = 64;

    const mode = std.meta.stringToEnum(enum { msdf, sdf }, args.next() orelse "msdf") orelse .msdf;
    const char: u21 = try std.unicode.utf8Decode(args.next() orelse "A");

    const font_path = args.next() orelse switch (builtin.os.tag) {
        .windows => "C:\\Windows\\Fonts\\arial.ttf",
        .linux => "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        .macos => "/Library/Fonts/Arial.ttf",
        else => {
            std.log.err("Please specify a font file path as an argument or use a default one for your OS.");
            return error.InvalidFontPath;
        },
    };

    const font_bytes = std.fs.cwd().readFileAlloc(
        allocator,
        font_path,
        std.math.maxInt(usize),
    ) catch |err| {
        std.log.err("Failed to read font file: {s}", .{font_path});
        return err;
    };
    defer allocator.free(font_bytes);

    const ttf: zmsdf.TrueType = try .load(font_bytes);

    var shape = zmsdf.Shape.init(allocator);
    defer shape.deinit();

    const glyph_a = ttf.codepointGlyphIndex(char) orelse {
        std.log.err("Glyph not found in font: {u}", .{char});
        return error.GlyphNotFound;
    };
    try zmsdf.font.loadGlyph(
        &ttf,
        allocator,
        &shape,
        glyph_a,
        .em_normalized,
    );
    try shape.normalize();
    try zmsdf.edgeColoringSimple(&shape, 3.0, 0, allocator);

    var my_bitmap: zmsdf.Bitmap(f32, 3) = try .init(allocator, size, size);
    defer my_bitmap.deinit();

    var config: zmsdf.MSDFGeneratorConfig = .default;
    config.generator.allocator = allocator;

    try zmsdf.generateMSDF(my_bitmap.ref(), &shape, .init(
        .init(.init(size, size), .init(0.125, 0.125)),
        .fromRange(.symmetrical(0.125)),
    ), config);

    // stb stuff
    {
        const stb = @import("stb");

        stb.init(allocator);
        defer stb.deinit();

        var img = try stb.Image.createEmpty(
            my_bitmap.data.width,
            my_bitmap.data.height,
            4,
            .{
                .bytes_per_component = 4,
                .bytes_per_row = my_bitmap.data.width * 4 * 4,
            },
        );
        defer img.deinit();

        const float_data_src = my_bitmap.data.slice();
        const float_data_dst = std.mem.bytesAsSlice(f32, img.data);
        for (0..my_bitmap.data.height) |y| {
            const row_dst = my_bitmap.data.height - y - 1;
            const row_src = y;
            for (0..my_bitmap.data.width) |x| {
                const src_index = (row_dst * my_bitmap.data.width + x) * 3;
                const dst_index = (row_src * my_bitmap.data.width + x) * 4;

                const medianed = zmsdf.util.median(
                    f64,
                    float_data_src[src_index + 0],
                    float_data_src[src_index + 1],
                    float_data_src[src_index + 2],
                );

                if (mode == .sdf) {
                    float_data_dst[dst_index + 0] = @floatCast(medianed);
                    float_data_dst[dst_index + 1] = @floatCast(medianed);
                    float_data_dst[dst_index + 2] = @floatCast(medianed);
                    float_data_dst[dst_index + 3] = 1.0;
                } else {
                    float_data_dst[dst_index + 0] = float_data_src[src_index + 0];
                    float_data_dst[dst_index + 1] = float_data_src[src_index + 1];
                    float_data_dst[dst_index + 2] = float_data_src[src_index + 2];
                    float_data_dst[dst_index + 3] = 1.0;
                }
            }
        }

        try img.writeToFile("glyph.hdr", .hdr);
    }
}
