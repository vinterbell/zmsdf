const std = @import("std");

pub fn Bitmap(comptime T: type, comptime Channels: u8) type {
    return struct {
        const Self = @This();

        allocator: std.mem.Allocator,
        data: BitmapRef(T, Channels),

        pub fn init(
            allocator: std.mem.Allocator,
            width: u32,
            height: u32,
        ) !Self {
            const data = try allocator.alloc(
                T,
                width * height * Channels,
            );
            return .{
                .allocator = allocator,
                .data = .{
                    .data = data,
                    .width = width,
                    .height = height,
                },
            };
        }

        pub fn deinit(self: *Self) void {
            self.allocator.free(self.data.data);
        }

        pub fn ref(self: *Self) BitmapRef(T, Channels) {
            return self.data;
        }

        pub fn constRef(self: *const Self) BitmapConstRef(T, Channels) {
            return self.data.constRef();
        }
    };
}

pub fn BitmapRef(comptime T: type, comptime Channels: u8) type {
    return struct {
        const Self = @This();

        data: []T,
        width: u32,
        height: u32,

        pub const zero: Self = .{
            .data = &.{},
            .width = 0,
            .height = 0,
        };

        pub fn slice(self: *const Self) []T {
            return self.data[0 .. self.width * self.height * @as(usize, Channels)];
        }

        pub fn byteSlice(self: *const Self) []u8 {
            return std.mem.bytesAsSlice(u8, self.slice());
        }

        pub fn getPixel(self: *const Self, x: u32, y: u32) *[Channels]T {
            return self.data[(y * self.width + x) * Channels ..][0..Channels];
        }

        pub fn get(self: *const Self, x: u32, y: u32) []T {
            return self.data[(y * self.width + x) * Channels ..];
        }

        pub fn constRef(self: *const Self) BitmapConstRef(T, Channels) {
            return .{
                .data = self.data,
                .width = self.width,
                .height = self.height,
            };
        }
    };
}

pub fn BitmapConstRef(comptime T: type, comptime Channels: u8) type {
    return struct {
        const Self = @This();

        data: []const T,
        width: u32,
        height: u32,

        pub const zero: Self = .{
            .data = &.{},
            .width = 0,
            .height = 0,
        };

        pub fn slice(self: *const Self) []const T {
            return self.data[0 .. self.width * self.height * @as(usize, Channels)];
        }

        pub fn byteSlice(self: *const Self) []const u8 {
            return std.mem.bytesAsSlice(u8, self.slice());
        }

        pub fn getPixel(self: *const Self, x: u32, y: u32) [Channels]T {
            return self.data[(y * self.width + x) * Channels ..][0..Channels].*;
        }

        pub fn get(self: *const Self, x: u32, y: u32) []const T {
            return self.data[(y * self.width + x) * Channels ..];
        }
    };
}
