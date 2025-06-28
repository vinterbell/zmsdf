//! Represents a range between two real values.
//! For example, the range of representable signed distances.
const Range = @This();

lower: f64,
upper: f64,

pub fn symmetrical(width: f64) Range {
    return .{ .lower = -0.5 * width, .upper = 0.5 * width };
}

pub fn init(lower: f64, upper: f64) Range {
    return .{ .lower = lower, .upper = upper };
}

pub fn multiply(self: Range, factor: f64) Range {
    return .{ .lower = self.lower * factor, .upper = self.upper * factor };
}

pub fn divide(self: Range, divisor: f64) Range {
    return .{ .lower = self.lower / divisor, .upper = self.upper / divisor };
}
