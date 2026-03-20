#include "topology.h"
#include "ring.h"
#include <cmath>
#include <vector>

// ---------------------------------------------------------------------------
// Cross product of vectors (b-a) and (c-a)  — 2D "z" component
// ---------------------------------------------------------------------------
static double cross2d(double ax, double ay,
                      double bx, double by,
                      double cx, double cy) {
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
}

// ---------------------------------------------------------------------------
// segments_intersect
//
// Uses the standard cross-product orientation test.
// Returns true only for *proper* intersections (interior × interior).
// Collinear / endpoint-touching cases return false so that adjacent
// edges that share a vertex are not flagged as intersecting.
// ---------------------------------------------------------------------------
bool segments_intersect(double p1x, double p1y, double p2x, double p2y,
                        double q1x, double q1y, double q2x, double q2y) {
    double d1 = cross2d(q1x, q1y, q2x, q2y, p1x, p1y);
    double d2 = cross2d(q1x, q1y, q2x, q2y, p2x, p2y);
    double d3 = cross2d(p1x, p1y, p2x, p2y, q1x, q1y);
    double d4 = cross2d(p1x, p1y, p2x, p2y, q2x, q2y);

    // Proper intersection: the two segments straddle each other
    if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
        ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0)))
        return true;

    // Collinear cases — treat as non-intersecting for shared-endpoint edges
    // (the spec only forbids proper crossings for topological validity)
    return false;
}

// ---------------------------------------------------------------------------
// topology_valid
//
// Naive O(n) scan.  For each ring, walk every edge and test it against
// the two new edges A→E and E→D.
//
// Edges we must SKIP (they are being removed or share an endpoint):
//   A→B, B→C, C→D   (deleted by the collapse)
//   The reverse edges that share A or D as an endpoint are allowed to
//   touch at the shared endpoint — segments_intersect already handles
//   this because it returns false for endpoint-only touches.
// ---------------------------------------------------------------------------
bool topology_valid(Vertex* B, const std::vector<Ring*>& rings) {
    Vertex* A = B->prev;
    Vertex* C = B->next;
    Vertex* D = C->next;

    double ex = B->ex, ey = B->ey;

    // The two new edges
    double ae_x1 = A->x, ae_y1 = A->y, ae_x2 = ex,   ae_y2 = ey;
    double ed_x1 = ex,   ed_y1 = ey,   ed_x2 = D->x, ed_y2 = D->y;

    for (Ring* ring : rings) {
        if (!ring->head) continue;
        Vertex* cur = ring->head;
        do {
            Vertex* nxt = cur->next;

            // Skip the three edges being deleted
            if ((cur == A && nxt == B) ||
                (cur == B && nxt == C) ||
                (cur == C && nxt == D)) {
                cur = nxt;
                continue;
            }

            double ex1 = cur->x, ey1 = cur->y;
            double ex2 = nxt->x, ey2 = nxt->y;

            if (segments_intersect(ae_x1, ae_y1, ae_x2, ae_y2,
                                   ex1, ey1, ex2, ey2))
                return false;

            if (segments_intersect(ed_x1, ed_y1, ed_x2, ed_y2,
                                   ex1, ey1, ex2, ey2))
                return false;

            cur = nxt;
        } while (cur != ring->head);
    }
    return true;
}
