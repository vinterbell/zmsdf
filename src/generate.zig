const std = @import("std");
const zmsdf = @import("root.zig");

/// The configuration of the MSDF error correction pass.
pub const ErrorCorrectionConfig = struct {
    mode: Mode,
    distance_check_mode: DistanceCheckMode,
    min_deviation_ratio: f64,
    min_improve_ratio: f64,

    pub const default_min_deviation_ratio: f64 = 0.01;
    pub const default_min_improve_ratio: f64 = 0.5;

    pub const Mode = enum {
        /// Skips error correction pass.
        disabled,
        /// Corrects all discontinuities of the distance field regardless if edges are adversely affected.
        indiscriminate,
        /// Corrects artifacts at edges and other discontinuous distances only if it does not affect edges or corners.
        edge_priority,
        /// Only corrects artifacts at edges.
        edge_only,
    };

    /// Configuration of whether to use an algorithm that computes the exact shape distance at the positions of suspected artifacts. This algorithm can be much slower.
    pub const DistanceCheckMode = enum {
        /// Never computes exact shape distance.
        do_not_check_distance,
        /// Only computes exact shape distance at edges. Provides a good balance between speed and precision.
        check_distance_at_edge,
        /// Computes and compares the exact shape distance for each suspected artifact.
        always_check_distance,
    };

    pub const default: ErrorCorrectionConfig = .{
        .mode = .edge_priority,
        .distance_check_mode = .check_distance_at_edge,
        .min_deviation_ratio = default_min_deviation_ratio,
        .min_improve_ratio = default_min_improve_ratio,
    };
};

/// The configuration of the distance field generator algorithm.
pub const GeneratorConfig = struct {
    allocator: std.mem.Allocator,
    /// Specifies whether to use the version of the algorithm that supports overlapping contours with the same winding. May be set to false to improve performance when no such contours are present.
    overlap_support: bool,

    /// Must set allocator before using this config.
    pub const default: GeneratorConfig = .{
        .allocator = undefined, // Must be set by the caller.
        .overlap_support = true,
    };
};

pub const MSDFGeneratorConfig = struct {
    generator: GeneratorConfig,
    error_correction: ErrorCorrectionConfig,

    /// Must set allocator in `generator` before using this config.
    pub const default: MSDFGeneratorConfig = .{
        .generator = .default,
        .error_correction = .default,
    };
};

pub fn DistancePixelConversion(
    comptime DistanceType: type,
) type {
    return struct {
        const Self = @This();

        pub const BitmapRefType = switch (DistanceType) {
            f64 => zmsdf.BitmapRef(f32, 1),
            zmsdf.MultiDistance => zmsdf.BitmapRef(f32, 3),
            zmsdf.MultiAndTrueDistance => zmsdf.BitmapRef(f32, 4),
            else => @compileError("Unsupported DistanceType " ++ @typeName(DistanceType)),
        };

        mapping: zmsdf.DistanceMapping,

        pub fn getPixelFloat(self: Self, distance: f32, out_pixels: []f32) void {
            out_pixels[0] = @floatCast(self.mapping.of(distance));
        }

        pub fn getPixelMultiDistance(self: Self, distance: zmsdf.MultiDistance, out_pixels: []f32) void {
            out_pixels[0] = @floatCast(self.mapping.of(distance.r));
            out_pixels[1] = @floatCast(self.mapping.of(distance.g));
            out_pixels[2] = @floatCast(self.mapping.of(distance.b));
        }

        pub fn getPixelMultiAndTrueDistance(self: Self, distance: zmsdf.MultiAndTrueDistance, out_pixels: []f32) void {
            out_pixels[0] = @floatCast(self.mapping.of(distance.r));
            out_pixels[1] = @floatCast(self.mapping.of(distance.g));
            out_pixels[2] = @floatCast(self.mapping.of(distance.b));
            out_pixels[3] = @floatCast(self.mapping.of(distance.a));
        }

        pub fn getPixel(self: Self, distance: anytype, out_pixels: []f32) void {
            switch (DistanceType) {
                f64 => self.getPixelFloat(distance, out_pixels),
                zmsdf.MultiDistance => self.getPixelMultiDistance(distance, out_pixels),
                zmsdf.MultiAndTrueDistance => self.getPixelMultiAndTrueDistance(distance, out_pixels),
                else => @compileError("Unsupported DistanceType " ++ @typeName(DistanceType)),
            }
        }
    };
}

pub fn generateDistanceField(
    comptime ContourCombiner: type,
    output: DistancePixelConversion(ContourCombiner.DistanceType).BitmapRefType,
    shape: *const zmsdf.Shape,
    transformation: zmsdf.SDFTransformation,
    temp: std.mem.Allocator,
) !void {
    const distance_pixel_conversion: DistancePixelConversion(ContourCombiner.DistanceType) = .{
        .mapping = transformation.distance_mapping,
    };
    {
        var distance_finder: zmsdf.ShapeDistanceFinder(ContourCombiner) = try .init(temp, shape);
        defer distance_finder.deinit();
        var right_to_left = false;

        for (0..output.height) |y| {
            const row = if (shape.inverse_y_axis) output.height - y - 1 else y;
            for (0..output.width) |col| {
                const u32_col: u32 = @intCast(col);
                const x = if (right_to_left) output.width - u32_col - 1 else u32_col;
                const f32_x: f32 = @floatFromInt(x);
                const f32_y: f32 = @floatFromInt(row);
                const p = transformation.unproject(.init(f32_x + 0.5, f32_y + 0.5));
                const distance = distance_finder.distance(p);
                const pixel = output.get(x, @intCast(row));
                distance_pixel_conversion.getPixel(distance, pixel);
            }
            right_to_left = !right_to_left;
        }
    }
}

pub fn generateSDF(
    output: zmsdf.BitmapRef(f32, 1),
    shape: *const zmsdf.Shape,
    transformation: zmsdf.SDFTransformation,
    config: GeneratorConfig,
) !void {
    if (config.generator.overlap_support) {
        try generateDistanceField(
            zmsdf.OverlappingTrueDistanceContourCombiner,
            output,
            shape,
            transformation,
            config.generator.allocator,
        );
    } else {
        try generateDistanceField(
            zmsdf.SimpleTrueDistanceContourCombiner,
            output,
            shape,
            transformation,
            config.generator.allocator,
        );
    }

    // TODO: Implement error correction for MSDF
}

pub fn generateMSDF(
    output: zmsdf.BitmapRef(f32, 3),
    shape: *const zmsdf.Shape,
    transformation: zmsdf.SDFTransformation,
    config: MSDFGeneratorConfig,
) !void {
    if (config.generator.overlap_support) {
        try generateDistanceField(
            zmsdf.OverlappingMultiDistanceContourCombiner,
            output,
            shape,
            transformation,
            config.generator.allocator,
        );
    } else {
        try generateDistanceField(
            zmsdf.SimpleMultiDistanceContourCombiner,
            output,
            shape,
            transformation,
            config.generator.allocator,
        );
    }

    // TODO: Implement error correction for MSDF
}

pub fn generateMTSDF(
    output: zmsdf.BitmapRef(f32, 4),
    shape: *const zmsdf.Shape,
    transformation: zmsdf.SDFTransformation,
    config: MSDFGeneratorConfig,
) !void {
    if (config.generator.overlap_support) {
        try generateDistanceField(
            zmsdf.OverlappingMultiAndTrueDistanceContourCombiner,
            output,
            shape,
            transformation,
            config.generator.allocator,
        );
    } else {
        try generateDistanceField(
            zmsdf.SimpleMultiAndTrueDistanceContourCombiner,
            output,
            shape,
            transformation,
            config.generator.allocator,
        );
    }

    // TODO: Implement error correction for MTSDF
}
