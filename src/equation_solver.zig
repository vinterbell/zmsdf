const std = @import("std");

const cpi = 3.14159265358979323846;

pub fn solveQuadratic(a: f64, b: f64, c: f64) std.BoundedArray(f64, 2) {
    if (a == 0 or @abs(b) > 1e12 * @abs(a)) {
        // a == 0, b == 0 -> no solution
        if (b == 0) {
            if (c == 0)
                return .{}; // 0 == 0
            return .{};
        }
        // return .{-c / b, 0};
        return .{
            .buffer = .{ -c / b, 0 },
            .len = 1,
        };
    }

    // d = b^2 - 4ac
    const discr = b * b - 4 * a * c;
    if (discr > 0) {
        const sqrt_dscr = @sqrt(discr);
        return .{
            .buffer = .{
                (-b + sqrt_dscr) / (2.0 * a),
                (-b - sqrt_dscr) / (2.0 * a),
            },
            .len = 2,
        };
    } else if (discr == 0) {
        return .{
            .buffer = .{ -b / (2.0 * a), 0 },
            .len = 1,
        };
    } else {
        return .{}; // No real solutions
    }
}

pub fn solveCubic(a: f64, b: f64, c: f64, d: f64) std.BoundedArray(f64, 3) {
    if (a != 0) {
        const bn = b / a;
        if (@abs(bn) < 1e6) { // Above this ratio, the numerical error gets larger than if we treated a as zero
            return solveCubicNormed(bn, c / a, d / a);
        }
    }
    const results = solveQuadratic(b, c, d);
    return .{
        .buffer = .{ results.buffer[0], results.buffer[1], 0.0 },
        .len = results.len,
    };
}

fn solveCubicNormed(_a: f64, b: f64, c: f64) std.BoundedArray(f64, 3) {
    var a = _a;
    const a2 = a * a;
    const q = 1.0 / 9.0 * (a2 - 3.0 * b);
    const r = 1.0 / 54.0 * (a * (2.0 * a2 - 9.0 * b) + 27 * c);
    const r2 = r * r;
    const q3 = q * q * q;
    a *= 1.0 / 3.0;

    if (r2 < q3) {
        var t = r / @sqrt(q3);
        if (t < -1) t = -1;
        if (t > 1) t = 1;
        t = std.math.acos(t);
        const q_negated = -2.0 * @sqrt(q);
        return .{
            .buffer = .{
                q_negated * @cos(1.0 / 3.0 * t) - a,
                q_negated * @cos(1.0 / 3.0 * (t + 2.0 * cpi)) - a,
                q_negated * @cos(1.0 / 3.0 * (t - 2.0 * cpi)) - a,
            },
            .len = 3,
        };
    } else {
        const u = @as(f64, if (r < 0) 1 else -1) * std.math.pow(f64, @abs(r) + @sqrt(r2 - q3), 1.0 / 3.0);
        const v = if (u == 0) 0 else q / u;
        const sol_one = (u + v) - a;
        if (u == v or @abs(u - v) < 1e-12 * @abs(u + v)) {
            const sol_two = -0.5 * (u + v) - a;
            return .{
                .buffer = .{ sol_one, sol_two, 0.0 },
                .len = 2,
            };
        }
        return .{
            .buffer = .{ sol_one, 0.0, 0.0 },
            .len = 1,
        };
    }
}
