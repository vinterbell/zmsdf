const std = @import("std");

pub fn build(b: *std.Build) !void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const zon = @import("build.zig.zon");
    const version = std.SemanticVersion.parse(zon.version) catch |err| {
        std.log.err("Failed to parse version: {}", .{err});
        return err;
    };

    const build_example = b.option(bool, "build-example", "Build the example executable") orelse false;

    // new rewrite
    const zmsdf = b.addModule("zmsdf", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });

    if (build_example) {
        const ztb_dep = b.lazyDependency("stb", .{
            .target = target,
            .optimize = optimize,
        }) orelse return;
        const example_mod = b.createModule(.{
            .root_source_file = b.path("src/example/main.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "zmsdf", .module = zmsdf },
                .{ .name = "stb", .module = ztb_dep.module("zstb") },
            },
        });
        const example_exe = b.addExecutable(.{
            .name = "example",
            .root_module = example_mod,
            .optimize = optimize,
            .version = version,
        });
        b.installArtifact(example_exe);

        const run_example = b.addRunArtifact(example_exe);
        const run_step = b.step("run-example", "Run the example");
        run_step.dependOn(&run_example.step);
        run_step.dependOn(b.getInstallStep());

        if (b.args) |args| run_example.addArgs(args);
    }

    const test_exe = b.addTest(.{
        .root_module = zmsdf,
        .target = target,
        .optimize = optimize,
    });
    b.installArtifact(test_exe);

    const test_step = b.step("test", "Run the tests");
    test_step.dependOn(&b.addRunArtifact(test_exe).step);
}
