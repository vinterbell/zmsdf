const std = @import("std");
const zmsdf = @import("root.zig");

const util = @import("util.zig");

pub const edge_length_precision = 4;

/// Assigns colors to edges of the shape in accordance to the multi-channel distance field technique.
/// May split some edges if necessary.
/// angleThreshold specifies the maximum angle (in radians) to be considered a corner, for example 3 (~172 degrees).
/// Values below 1/2 PI will be treated as the external angle.
pub fn edgeColoringSimple(
    shape: *zmsdf.Shape,
    angle_threshold: f64,
    in_seed: u64,
    temp: std.mem.Allocator,
) !void {
    var seed = in_seed;
    const cross_threshold = std.math.sin(angle_threshold);
    var color: zmsdf.EdgeColor = initColor(&seed);
    var corners: std.ArrayListUnmanaged(i32) = .empty;
    defer corners.deinit(temp);

    for (shape.contours.items) |*contour| {
        if (contour.edges.items.len == 0) continue;
        { // Identify corners
            corners.clearRetainingCapacity();
            var prev_direction = contour.edges.getLast().direction(1.0);
            for (contour.edges.items, 0..) |edge, index| {
                if (isCorner(
                    prev_direction.normalize(false),
                    edge.direction(0).normalize(false),
                    cross_threshold,
                )) {
                    try corners.append(temp, @intCast(index));
                }
                prev_direction = edge.direction(1);
            }
        }

        // smooth contour
        if (corners.items.len == 0) {
            switchColor(&color, &seed);
            for (contour.edges.items) |*edge| {
                edge.color = color;
            }
        } else if (corners.items.len == 1) {
            // "teardrop" case
            switchColor(&color, &seed);
            var colors: [3]zmsdf.EdgeColor = .{
                color,
                .white,
                .white,
            };
            switchColor(&color, &seed);
            colors[2] = color;
            const corner: usize = @intCast(corners.items[0]);
            if (contour.edges.items.len >= 3) {
                const m: usize = @intCast(contour.edges.items.len);
                for (0..m) |i| {
                    contour.edges.items[(corner + i) % m].color = colors[
                        1 + symmetricalTrichotomy(
                            @intCast(i),
                            @intCast(m),
                        )
                    ];
                }
            } else if (contour.edges.items.len >= 1) {
                // Less than three edge segments for three colors => edges must be split
                var parts: [7]?zmsdf.EdgeSegment = @splat(null);
                const part_set_1 = contour.edges.items[0].splitInThirds();
                parts[0 + 3 * corner] = part_set_1[0];
                parts[1 + 3 * corner] = part_set_1[1];
                parts[2 + 3 * corner] = part_set_1[2];
                if (contour.edges.items.len >= 2) {
                    const part_set_2 = contour.edges.items[1].splitInThirds();
                    parts[3 - 3 * corner] = part_set_2[0];
                    parts[4 - 3 * corner] = part_set_2[1];
                    parts[5 - 3 * corner] = part_set_2[2];

                    parts[1].?.color = colors[0];
                    parts[3].?.color = colors[1];
                    parts[5].?.color = colors[2];

                    parts[0].?.color = colors[0];
                    parts[2].?.color = colors[1];
                    parts[4].?.color = colors[2];
                } else {
                    parts[0].?.color = colors[0];
                    parts[1].?.color = colors[1];
                    parts[2].?.color = colors[2];
                }
                contour.edges.clearRetainingCapacity();
                for (&parts) |part| {
                    try contour.edges.append(contour.allocator, part orelse break);
                }
            }
        } else {
            // multiple corners
            const corner_count: usize = corners.items.len;
            var spline: usize = 0;
            const start: usize = @intCast(corners.items[0]);
            const m: usize = contour.edges.items.len;
            switchColor(&color, &seed);
            const initial_color: zmsdf.EdgeColor = color;
            for (0..m) |i| {
                const index = (start + i) % m;
                if (spline + 1 < corner_count and corners.items[spline + 1] == index) {
                    spline += 1;
                    switchColorBanned(&color, &seed, if (spline == corner_count - 1) initial_color else .black);
                }
                contour.edges.items[index].color = color;
            }
        }
    }
}

/// The alternative "ink trap" coloring strategy is designed for better results with typefaces
/// that use ink traps as a design feature. It guarantees that even if all edges that are shorter than
/// both their neighboring edges are removed, the coloring remains consistent with the established rules.
pub fn edgeColoringInkTrap(
    shape: *zmsdf.Shape,
    angle_threshold: f64,
    seed: u64,
) void {
    _ = shape;
    _ = angle_threshold;
    _ = seed;
}

/// The alternative coloring by distance tries to use different colors for edges that are close together.
/// This should theoretically be the best strategy on average. However, since it needs to compute the
/// distance between all pairs of edges, and perform a graph optimization task, it is much slower
/// than the rest.
pub fn edgeColoringByDistance(
    shape: *zmsdf.Shape,
    angle_threshold: f64,
    seed: u64,
) void {
    _ = shape;
    _ = angle_threshold;
    _ = seed;
}

// For each position < n, this function will return -1, 0, or 1,
// depending on whether the position is closer to the beginning, middle, or end, respectively.
// It is guaranteed that the output will be balanced in that the total for positions 0 through n-1 will be zero.
fn symmetricalTrichotomy(position: i32, n: i32) usize {
    const position_float: f64 = @floatFromInt(position);
    const n_float: f64 = @floatFromInt(n);
    const value = (3 + 2.875 * position_float / (n_float - 1) - 1.4375 + 0.5);
    const value_int: i32 = @intFromFloat(@trunc(value));
    return @intCast(value_int - 3);
}

fn isCorner(
    a_dir: zmsdf.Vector2,
    b_dir: zmsdf.Vector2,
    cross_threshold: f64,
) bool {
    return zmsdf.Vector2.dot(a_dir, b_dir) <= 0 or
        @abs(zmsdf.Vector2.cross(a_dir, b_dir)) > cross_threshold;
}

fn estimateLength(
    edge: zmsdf.EdgeSegment,
) f64 {
    var length: f64 = 0.0;
    var prev = edge.point(0.0);
    for (0..edge_length_precision) |i| {
        const cur = edge.point(@as(f64, @floatFromInt(i + 1)) / @as(f64, @floatFromInt(edge_length_precision)));
        length += (cur.subtract(prev)).length();
        prev = cur;
    }
    return length;
}

fn seedExtract2(seed: *u64) i32 {
    const v: i32 = @intCast(seed.* & 1);
    seed.* >>= 1;
    return v;
}

fn seedExtract3(seed: *u64) i32 {
    const v: i21 = @intCast(seed.* % 3);
    seed.* /= 3;
    return v;
}

fn initColor(seed: *u64) zmsdf.EdgeColor {
    const colors: [3]zmsdf.EdgeColor = .{ .cyan, .magenta, .yellow };
    return colors[@intCast(seedExtract3(seed))];
}

fn switchColor(
    color: *zmsdf.EdgeColor,
    seed: *u64,
) void {
    const shifted = @as(u8, color.toInt() << @as(u3, @intCast(1 + seedExtract2(seed))));
    color.* = .fromInt(shifted | (shifted >> 3));
}

fn switchColorBanned(
    color: *zmsdf.EdgeColor,
    seed: *u64,
    banned: zmsdf.EdgeColor,
) void {
    var combined: zmsdf.EdgeColor = .{
        .red_channel = (color.red_channel and banned.red_channel),
        .green_channel = (color.green_channel and banned.green_channel),
        .blue_channel = (color.blue_channel and banned.blue_channel),
    };
    if (combined.eql(.red) or combined.eql(.green) or combined.eql(.blue)) {
        color.* = .fromInt(combined.toInt() ^ zmsdf.EdgeColor.white.toInt());
    } else {
        switchColor(color, seed);
    }
}

// EDGE COLORING BY DISTANCE - EXPERIMENTAL IMPLEMENTATION - WORK IN PROGRESS
// const max_recolor_steps = 16;
// const edge_distance_precision = 16;

// fn edgeToEdgeDistance(a: *const zmsdf.EdgeSegment, b: *const zmsdf.EdgeSegment, precision: u32) f64 {
//     if (a.point(0).eql(b.point(0)) or
//         a.point(0).eql(b.point(1)) or
//         a.point(1).eql(b.point(0)) or
//         a.point(1).eql(b.point(1)))
//     {
//         return 0.0;
//     }

//     const i_fac: f64 = 1.0 / @as(f64, @floatFromInt(precision));
//     var min_distance = b.point(0).subtract(a.point(0)).length();
//     for (0..precision + 1) |i| {
//         const t: f64 = i_fac * @as(f64, @floatFromInt(i));
//         const d = @abs(a.signedDistance(b.point(t), t).distance.distance);
//         min_distance = @min(min_distance, d);
//     }
//     for (0..precision + 1) |i| {
//         const t: f64 = i_fac * @as(f64, @floatFromInt(i));
//         const d = @abs(b.signedDistance(a.point(t), t).distance.distance);
//         min_distance = @min(min_distance, d);
//     }
//     return min_distance;
// }

// fn splineToSplineDistance(a_segments: []const zmsdf.EdgeSegment, b_segments: []const zmsdf.EdgeSegment, precision: u32) f64 {
//     var min_distance: f64 = util.f64_max;
//     for (a_segments) |a| {
//         for (b_segments) |b| {
//             const d = edgeToEdgeDistance(&a, &b, precision);
//             min_distance = @min(min_distance, d);
//         }
//     }
//     return min_distance;
// }

// inline fn invertInt(value: anytype) @TypeOf(value) {
//     return if (value == 0) return 1 else 0;
// }

// fn colorSecondDegreeGraph(coloring: []u8, edge_matrix: []const []const u8, vertex_count: usize, seed_constant: u64) void {
//     var seed = seed_constant;
//     for (0..vertex_count) |i| {
//         var possible_colors: u8 = 7;
//         for (0..i) |j| {
//             if (edge_matrix[i][j] != 0) {
//                 possible_colors &= ~(1 << @intCast(coloring[j]));
//             }
//         }
//         var color: zmsdf.EdgeColor = .black;
//         switch (possible_colors) {
//             1 => color = .fromInt(0),
//             2 => color = .fromInt(1),
//             3 => color = .fromInt(seedExtract2(&seed)), // 0 or 1
//             4 => color = .fromInt(2),
//             5 => color = .fromInt((@as(u8, invertInt(seedExtract2(&seed))) << 1)), // 2 or 0
//             6 => color = .fromInt(seedExtract2(&seed) + 1), // 1 or 2
//             7 => color = .fromInt((seedExtract3(&seed) + i) % 3), // 0 or 1 or 2
//             else => unreachable,
//         }
//         coloring[i] = color.toInt();
//     }
// }

// fn vertexPossibleColors(coloring: []const u8, edge_vector: []const u8, vertex_count: usize) u8 {
//     var used_colors: u8 = 0;
//     for (0..vertex_count) |i| {
//         if (edge_vector[i] != 0) {
//             used_colors |= 1 << @intCast(coloring[i]);
//         }
//     }
//     return 7 & ~used_colors;
// }

// static void uncolorSameNeighbors(std::queue<int> &uncolored, int *coloring, const int *const *edgeMatrix, int vertex, int vertexCount) {
//     for (int i = vertex+1; i < vertexCount; ++i) {
//         if (edgeMatrix[vertex][i] && coloring[i] == coloring[vertex]) {
//             coloring[i] = -1;
//             uncolored.push(i);
//         }
//     }
//     for (int i = 0; i < vertex; ++i) {
//         if (edgeMatrix[vertex][i] && coloring[i] == coloring[vertex]) {
//             coloring[i] = -1;
//             uncolored.push(i);
//         }
//     }
// }

// static bool tryAddEdge(int *coloring, int *const *edgeMatrix, int vertexCount, int vertexA, int vertexB, int *coloringBuffer) {
//     static const int FIRST_POSSIBLE_COLOR[8] = { -1, 0, 1, 0, 2, 2, 1, 0 };
//     edgeMatrix[vertexA][vertexB] = 1;
//     edgeMatrix[vertexB][vertexA] = 1;
//     if (coloring[vertexA] != coloring[vertexB])
//         return true;
//     int bPossibleColors = vertexPossibleColors(coloring, edgeMatrix[vertexB], vertexCount);
//     if (bPossibleColors) {
//         coloring[vertexB] = FIRST_POSSIBLE_COLOR[bPossibleColors];
//         return true;
//     }
//     memcpy(coloringBuffer, coloring, sizeof(int)*vertexCount);
//     std::queue<int> uncolored;
//     {
//         int *coloring = coloringBuffer;
//         coloring[vertexB] = FIRST_POSSIBLE_COLOR[7&~(1<<coloring[vertexA])];
//         uncolorSameNeighbors(uncolored, coloring, edgeMatrix, vertexB, vertexCount);
//         int step = 0;
//         while (!uncolored.empty() && step < MAX_RECOLOR_STEPS) {
//             int i = uncolored.front();
//             uncolored.pop();
//             int possibleColors = vertexPossibleColors(coloring, edgeMatrix[i], vertexCount);
//             if (possibleColors) {
//                 coloring[i] = FIRST_POSSIBLE_COLOR[possibleColors];
//                 continue;
//             }
//             do {
//                 coloring[i] = step++%3;
//             } while (edgeMatrix[i][vertexA] && coloring[i] == coloring[vertexA]);
//             uncolorSameNeighbors(uncolored, coloring, edgeMatrix, i, vertexCount);
//         }
//     }
//     if (!uncolored.empty()) {
//         edgeMatrix[vertexA][vertexB] = 0;
//         edgeMatrix[vertexB][vertexA] = 0;
//         return false;
//     }
//     memcpy(coloring, coloringBuffer, sizeof(int)*vertexCount);
//     return true;
// }

// static int cmpDoublePtr(const void *a, const void *b) {
//     return sign(**reinterpret_cast<const double *const *>(a)-**reinterpret_cast<const double *const *>(b));
// }

// void edgeColoringByDistance(Shape &shape, double angleThreshold, unsigned long long seed) {

//     std::vector<EdgeSegment *, Allocator<EdgeSegment *>> edgeSegments;
//     std::vector<int, Allocator<int>> splineStarts;

//     double crossThreshold = sin(angleThreshold);
//     std::vector<int, Allocator<int>> corners;
//     for (std::vector<Contour, Allocator<Contour>>::iterator contour = shape.contours.begin(); contour != shape.contours.end(); ++contour)
//         if (!contour->edges.empty()) {
//             // Identify corners
//             corners.clear();
//             Vector2 prevDirection = contour->edges.back()->direction(1);
//             int index = 0;
//             for (std::vector<EdgeHolder, Allocator<EdgeHolder>>::const_iterator edge = contour->edges.begin(); edge != contour->edges.end(); ++edge, ++index) {
//                 if (isCorner(prevDirection.normalize(), (*edge)->direction(0).normalize(), crossThreshold))
//                     corners.push_back(index);
//                 prevDirection = (*edge)->direction(1);
//             }

//             splineStarts.push_back((int) edgeSegments.size());
//             // Smooth contour
//             if (corners.empty())
//                 for (std::vector<EdgeHolder, Allocator<EdgeHolder>>::iterator edge = contour->edges.begin(); edge != contour->edges.end(); ++edge)
//                     edgeSegments.push_back(&**edge);
//             // "Teardrop" case
//             else if (corners.size() == 1) {
//                 int corner = corners[0];
//                 if (contour->edges.size() >= 3) {
//                     int m = (int) contour->edges.size();
//                     for (int i = 0; i < m; ++i) {
//                         if (i == m/2)
//                             splineStarts.push_back((int) edgeSegments.size());
//                         if (symmetricalTrichotomy(i, m))
//                             edgeSegments.push_back(&*contour->edges[(corner+i)%m]);
//                         else
//                             contour->edges[(corner+i)%m]->color = WHITE;
//                     }
//                 } else if (contour->edges.size() >= 1) {
//                     // Less than three edge segments for three colors => edges must be split
//                     EdgeSegment *parts[7] = { };
//                     contour->edges[0]->splitInThirds(parts[0+3*corner], parts[1+3*corner], parts[2+3*corner]);
//                     if (contour->edges.size() >= 2) {
//                         contour->edges[1]->splitInThirds(parts[3-3*corner], parts[4-3*corner], parts[5-3*corner]);
//                         edgeSegments.push_back(parts[0]);
//                         edgeSegments.push_back(parts[1]);
//                         parts[2]->color = parts[3]->color = WHITE;
//                         splineStarts.push_back((int) edgeSegments.size());
//                         edgeSegments.push_back(parts[4]);
//                         edgeSegments.push_back(parts[5]);
//                     } else {
//                         edgeSegments.push_back(parts[0]);
//                         parts[1]->color = WHITE;
//                         splineStarts.push_back((int) edgeSegments.size());
//                         edgeSegments.push_back(parts[2]);
//                     }
//                     contour->edges.clear();
//                     for (int i = 0; parts[i]; ++i)
//                         contour->edges.push_back(EdgeHolder(parts[i]));
//                 }
//             }
//             // Multiple corners
//             else {
//                 int cornerCount = (int) corners.size();
//                 int spline = 0;
//                 int start = corners[0];
//                 int m = (int) contour->edges.size();
//                 for (int i = 0; i < m; ++i) {
//                     int index = (start+i)%m;
//                     if (spline+1 < cornerCount && corners[spline+1] == index) {
//                         splineStarts.push_back((int) edgeSegments.size());
//                         ++spline;
//                     }
//                     edgeSegments.push_back(&*contour->edges[index]);
//                 }
//             }
//         }
//     splineStarts.push_back((int) edgeSegments.size());

//     int segmentCount = (int) edgeSegments.size();
//     int splineCount = (int) splineStarts.size()-1;
//     if (!splineCount)
//         return;

//     std::vector<double, Allocator<double>> distanceMatrixStorage(splineCount*splineCount);
//     std::vector<double *, Allocator<double *>> distanceMatrix(splineCount);
//     for (int i = 0; i < splineCount; ++i)
//         distanceMatrix[i] = &distanceMatrixStorage[i*splineCount];
//     const double *distanceMatrixBase = &distanceMatrixStorage[0];

//     for (int i = 0; i < splineCount; ++i) {
//         distanceMatrix[i][i] = -1;
//         for (int j = i+1; j < splineCount; ++j) {
//             double dist = splineToSplineDistance(&edgeSegments[0], splineStarts[i], splineStarts[i+1], splineStarts[j], splineStarts[j+1], EDGE_DISTANCE_PRECISION);
//             distanceMatrix[i][j] = dist;
//             distanceMatrix[j][i] = dist;
//         }
//     }

//     std::vector<const double *, Allocator<const double *>> graphEdgeDistances;
//     graphEdgeDistances.reserve(splineCount*(splineCount-1)/2);
//     for (int i = 0; i < splineCount; ++i)
//         for (int j = i+1; j < splineCount; ++j)
//             graphEdgeDistances.push_back(&distanceMatrix[i][j]);
//     int graphEdgeCount = (int) graphEdgeDistances.size();
//     if (!graphEdgeDistances.empty())
//         qsort(&graphEdgeDistances[0], graphEdgeDistances.size(), sizeof(const double *), &cmpDoublePtr);

//     std::vector<int, Allocator<int>> edgeMatrixStorage(splineCount*splineCount);
//     std::vector<int *, Allocator<int *>> edgeMatrix(splineCount);
//     for (int i = 0; i < splineCount; ++i)
//         edgeMatrix[i] = &edgeMatrixStorage[i*splineCount];
//     int nextEdge = 0;
//     for (; nextEdge < graphEdgeCount && !*graphEdgeDistances[nextEdge]; ++nextEdge) {
//         int elem = (int) (graphEdgeDistances[nextEdge]-distanceMatrixBase);
//         int row = elem/splineCount;
//         int col = elem%splineCount;
//         edgeMatrix[row][col] = 1;
//         edgeMatrix[col][row] = 1;
//     }

//     std::vector<int, Allocator<int>> coloring(2*splineCount);
//     colorSecondDegreeGraph(&coloring[0], &edgeMatrix[0], splineCount, seed);
//     for (; nextEdge < graphEdgeCount; ++nextEdge) {
//         int elem = (int) (graphEdgeDistances[nextEdge]-distanceMatrixBase);
//         tryAddEdge(&coloring[0], &edgeMatrix[0], splineCount, elem/splineCount, elem%splineCount, &coloring[splineCount]);
//     }

//     const EdgeColor colors[3] = { YELLOW, CYAN, MAGENTA };
//     int spline = -1;
//     for (int i = 0; i < segmentCount; ++i) {
//         if (splineStarts[spline+1] == i)
//             ++spline;
//         edgeSegments[i]->color = colors[coloring[spline]];
//     }
// }

// namespace zmsdf {

// struct EdgeColoringInkTrapCorner {
//     int index;
//     double prevEdgeLengthEstimate;
//     bool minor;
//     EdgeColor color;
// };

// void edgeColoringInkTrap(Shape &shape, double angleThreshold, unsigned long long seed) {
//     typedef EdgeColoringInkTrapCorner Corner;
//     double crossThreshold = sin(angleThreshold);
//     EdgeColor color = initColor(seed);
//     std::vector<Corner, Allocator<Corner>> corners;
//     for (std::vector<Contour, Allocator<Contour>>::iterator contour = shape.contours.begin(); contour != shape.contours.end(); ++contour) {
//         if (contour->edges.empty())
//             continue;
//         double splineLength = 0;
//         { // Identify corners
//             corners.clear();
//             Vector2 prevDirection = contour->edges.back()->direction(1);
//             int index = 0;
//             for (std::vector<EdgeHolder, Allocator<EdgeHolder>>::const_iterator edge = contour->edges.begin(); edge != contour->edges.end(); ++edge, ++index) {
//                 if (isCorner(prevDirection.normalize(), (*edge)->direction(0).normalize(), crossThreshold)) {
//                     Corner corner = { index, splineLength };
//                     corners.push_back(corner);
//                     splineLength = 0;
//                 }
//                 splineLength += estimateEdgeLength(*edge);
//                 prevDirection = (*edge)->direction(1);
//             }
//         }

//         // Smooth contour
//         if (corners.empty()) {
//             switchColor(color, seed);
//             for (std::vector<EdgeHolder, Allocator<EdgeHolder>>::iterator edge = contour->edges.begin(); edge != contour->edges.end(); ++edge)
//                 (*edge)->color = color;
//         }
//         // "Teardrop" case
//         else if (corners.size() == 1) {
//             EdgeColor colors[3];
//             switchColor(color, seed);
//             colors[0] = color;
//             colors[1] = WHITE;
//             switchColor(color, seed);
//             colors[2] = color;
//             int corner = corners[0].index;
//             if (contour->edges.size() >= 3) {
//                 int m = (int) contour->edges.size();
//                 for (int i = 0; i < m; ++i)
//                     contour->edges[(corner+i)%m]->color = colors[1+symmetricalTrichotomy(i, m)];
//             } else if (contour->edges.size() >= 1) {
//                 // Less than three edge segments for three colors => edges must be split
//                 EdgeSegment *parts[7] = { };
//                 contour->edges[0]->splitInThirds(parts[0+3*corner], parts[1+3*corner], parts[2+3*corner]);
//                 if (contour->edges.size() >= 2) {
//                     contour->edges[1]->splitInThirds(parts[3-3*corner], parts[4-3*corner], parts[5-3*corner]);
//                     parts[0]->color = parts[1]->color = colors[0];
//                     parts[2]->color = parts[3]->color = colors[1];
//                     parts[4]->color = parts[5]->color = colors[2];
//                 } else {
//                     parts[0]->color = colors[0];
//                     parts[1]->color = colors[1];
//                     parts[2]->color = colors[2];
//                 }
//                 contour->edges.clear();
//                 for (int i = 0; parts[i]; ++i)
//                     contour->edges.push_back(EdgeHolder(parts[i]));
//             }
//         }
//         // Multiple corners
//         else {
//             int cornerCount = (int) corners.size();
//             int majorCornerCount = cornerCount;
//             if (cornerCount > 3) {
//                 corners.begin()->prevEdgeLengthEstimate += splineLength;
//                 for (int i = 0; i < cornerCount; ++i) {
//                     if (
//                         corners[i].prevEdgeLengthEstimate > corners[(i+1)%cornerCount].prevEdgeLengthEstimate &&
//                         corners[(i+1)%cornerCount].prevEdgeLengthEstimate < corners[(i+2)%cornerCount].prevEdgeLengthEstimate
//                     ) {
//                         corners[i].minor = true;
//                         --majorCornerCount;
//                     }
//                 }
//             }
//             EdgeColor initialColor = BLACK;
//             for (int i = 0; i < cornerCount; ++i) {
//                 if (!corners[i].minor) {
//                     --majorCornerCount;
//                     switchColor(color, seed, EdgeColor(!majorCornerCount*initialColor));
//                     corners[i].color = color;
//                     if (!initialColor)
//                         initialColor = color;
//                 }
//             }
//             for (int i = 0; i < cornerCount; ++i) {
//                 if (corners[i].minor) {
//                     EdgeColor nextColor = corners[(i+1)%cornerCount].color;
//                     corners[i].color = EdgeColor((color&nextColor)^WHITE);
//                 } else
//                     color = corners[i].color;
//             }
//             int spline = 0;
//             int start = corners[0].index;
//             color = corners[0].color;
//             int m = (int) contour->edges.size();
//             for (int i = 0; i < m; ++i) {
//                 int index = (start+i)%m;
//                 if (spline+1 < cornerCount && corners[spline+1].index == index)
//                     color = corners[++spline].color;
//                 contour->edges[index]->color = color;
//             }
//         }
//     }
// }

// }
