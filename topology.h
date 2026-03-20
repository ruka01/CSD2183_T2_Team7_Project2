#pragma once
#include <vector>

struct Vertex;
class Ring;

// ---------------------------------------------------------------------------
// Topology checks
//
// Before applying any collapse A→B→C→D → A→E→D we must verify:
//   (a) The new edge A→E does not intersect any existing edge in any ring.
//   (b) The new edge E→D does not intersect any existing edge in any ring.
//   (c) Edges that SHARE an endpoint with A→E or E→D are allowed to touch
//       at that shared endpoint but must not otherwise cross.
//
// A naive O(n) check is correct and sufficient to get the project working.
// Replace with a spatial index (R-tree / grid) for large inputs.
// ---------------------------------------------------------------------------

// Returns true if segment (p1→p2) and segment (q1→q2) properly intersect
// (i.e. cross in their interiors, not just at shared endpoints).
bool segments_intersect(double p1x, double p1y, double p2x, double p2y,
                        double q1x, double q1y, double q2x, double q2y);

// Check whether inserting vertex E (at ex, ey) between A and D (so that
// the ring gains edges A→E and E→D and loses edges A→B, B→C, C→D) would
// create any intersection with any edge in `rings`.
//
// `collapse_B` is the vertex being collapsed (B in A→B→C→D); it and C are
// about to be removed.  Pass it so we can skip the edges being deleted.
//
// Returns true if the collapse is topologically valid (no intersections).
bool topology_valid(Vertex* collapse_B,
                    const std::vector<Ring*>& rings);
