const std = @import("std");

const generate = @import("generate.zig");
pub const ErrorCorrectionConfig = generate.ErrorCorrectionConfig;
pub const GeneratorConfig = generate.GeneratorConfig;
pub const MSDFGeneratorConfig = generate.MSDFGeneratorConfig;
pub const generateMSDF = generate.generateMSDF;

pub const Vector2 = @import("Vector2.zig");
const bitmap = @import("bitmap.zig");
pub const Bitmap = bitmap.Bitmap;
pub const BitmapConstRef = bitmap.BitmapConstRef;
pub const BitmapRef = bitmap.BitmapRef;
pub const Contour = @import("Contour.zig");

const contour_combiners = @import("contour_combiners.zig");

pub const SimpleContourCombiner = contour_combiners.SimpleContourCombiner;
pub const SimpleTrueDistanceContourCombiner = contour_combiners.SimpleTrueDistanceContourCombiner;
pub const SimplePerpendicularDistanceContourCombiner = contour_combiners.SimplePerpendicularDistanceContourCombiner;
pub const SimpleMultiDistanceContourCombiner = contour_combiners.SimpleMultiDistanceContourCombiner;
pub const SimpleMultiAndTrueDistanceContourCombiner = contour_combiners.SimpleMultiAndTrueDistanceContourCombiner;

pub const OverlappingContourCombiner = contour_combiners.OverlappingContourCombiner;
pub const OverlappingTrueDistanceContourCombiner = contour_combiners.OverlappingTrueDistanceContourCombiner;
pub const OverlappingPerpendicularDistanceContourCombiner = contour_combiners.OverlappingPerpendicularDistanceContourCombiner;
pub const OverlappingMultiDistanceContourCombiner = contour_combiners.OverlappingMultiDistanceContourCombiner;
pub const OverlappingMultiAndTrueDistanceContourCombiner = contour_combiners.OverlappingMultiAndTrueDistanceContourCombiner;

pub const DistanceMapping = @import("DistanceMapping.zig");
pub const EdgeSegment = @import("EdgeSegment.zig");

const edge_coloring = @import("edge_coloring.zig");
pub const edgeColoringSimple = edge_coloring.edgeColoringSimple;
pub const edgeColoringInkTrap = edge_coloring.edgeColoringInkTrap;
pub const edgeColoringByDistance = edge_coloring.edgeColoringByDistance;

const edge_selectors = @import("edge_selectors.zig");
pub const MultiAndTrueDistance = edge_selectors.MultiAndTrueDistance;
pub const MultiAndTrueDistanceSelector = edge_selectors.MultiAndTrueDistanceSelector;
pub const MultiDistance = edge_selectors.MultiDistance;
pub const MultiDistanceSelector = edge_selectors.MultiDistanceSelector;
pub const PerpendicularDistanceSelector = edge_selectors.PerpendicularDistanceSelector;
pub const PerpendicularDistanceSelectorBase = edge_selectors.PerpendicularDistanceSelectorBase;
pub const TrueDistanceSelector = edge_selectors.TrueDistanceSelector;

pub const Projection = @import("Projection.zig");

pub const Range = @import("Range.zig");
pub const Scanline = @import("Scanline.zig");
pub const SDFTransformation = @import("SDFTransformation.zig");
pub const SignedDistance = @import("SignedDistance.zig");
pub const Shape = @import("Shape.zig");
pub const ShapeDistanceFinder = @import("shape_distance_finder.zig").ShapeDistanceFinder;

pub const EdgeColor = packed struct(u8) {
    red_channel: bool,
    green_channel: bool,
    blue_channel: bool,
    _padding: u5 = 0,

    pub const black: EdgeColor = .{
        .red_channel = false,
        .green_channel = false,
        .blue_channel = false,
    };
    pub const red: EdgeColor = .{
        .red_channel = true,
        .green_channel = false,
        .blue_channel = false,
    };
    pub const green: EdgeColor = .{
        .red_channel = false,
        .green_channel = true,
        .blue_channel = false,
    };
    pub const yellow: EdgeColor = .{
        .red_channel = true,
        .green_channel = true,
        .blue_channel = false,
    };
    pub const blue: EdgeColor = .{
        .red_channel = false,
        .green_channel = false,
        .blue_channel = true,
    };
    pub const magenta: EdgeColor = .{
        .red_channel = true,
        .green_channel = false,
        .blue_channel = true,
    };
    pub const cyan: EdgeColor = .{
        .red_channel = false,
        .green_channel = true,
        .blue_channel = true,
    };
    pub const white: EdgeColor = .{
        .red_channel = true,
        .green_channel = true,
        .blue_channel = true,
    };

    pub fn toInt(self: EdgeColor) u8 {
        return @bitCast(self);
    }

    pub fn fromInt(value: u8) EdgeColor {
        var color: EdgeColor = @bitCast(value);
        color._padding = 0; // Ensure padding is zeroed
        return color;
    }

    pub fn eql(self: EdgeColor, other: EdgeColor) bool {
        return self.red_channel == other.red_channel and
            self.green_channel == other.green_channel and
            self.blue_channel == other.blue_channel;
    }
};

pub const Bounds = struct {
    left: f64,
    bottom: f64,
    right: f64,
    top: f64,

    pub const zero: Bounds = .{
        .left = 0,
        .bottom = 0,
        .right = 0,
        .top = 0,
    };
};

pub const Polarity = enum(i32) {
    pos = 1,
    neg = -1,
    zero = 0,

    pub inline fn of(x: anytype) Polarity {
        const result = std.math.sign(x);
        if (result > 0) {
            return .pos;
        } else if (result < 0) {
            return .neg;
        } else {
            return .zero;
        }
    }

    pub fn asInt(self: Polarity) i32 {
        return @intFromEnum(self);
    }

    pub fn asFloat(self: Polarity) f64 {
        return @floatFromInt(self.asInt());
    }

    pub inline fn invert(self: Polarity) Polarity {
        return switch (self) {
            .pos => .neg,
            .neg => .pos,
            .zero => .zero,
        };
    }

    pub inline fn isPositive(self: Polarity) bool {
        return self.asInt() > 0;
    }

    pub inline fn isNegative(self: Polarity) bool {
        return self.asInt() < 0;
    }
};

comptime {
    std.testing.refAllDeclsRecursive(@This());
}

// external
pub const TrueType = @import("external/TrueType.zig");
pub const font = @import("import_font.zig");

pub const util = @import("util.zig");
