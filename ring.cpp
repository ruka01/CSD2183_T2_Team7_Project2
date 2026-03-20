#include "ring.h"
#include <cassert>
#include <cmath>
#include <stdexcept>

// ---------------------------------------------------------------------------
// Internal geometry helpers (file-local)
// ---------------------------------------------------------------------------

// Signed area of triangle (ax,ay)→(bx,by)→(cx,cy)
static double tri_signed_area(double ax, double ay,
                               double bx, double by,
                               double cx, double cy) {
    return 0.5 * ((bx - ax) * (cy - ay) - (cx - ax) * (by - ay));
}

// ---------------------------------------------------------------------------
// Ring::append
// ---------------------------------------------------------------------------
void Ring::append(double x, double y) {
    Vertex* v = new Vertex();
    v->x      = x;
    v->y      = y;
    v->ring_id = id;

    if (!head) {
        // First vertex: point to itself
        head    = v;
        v->prev = v;
        v->next = v;
    } else {
        // Insert before head (i.e. at the tail of the sequence)
        Vertex* tail = head->prev;
        tail->next   = v;
        v->prev      = tail;
        v->next      = head;
        head->prev   = v;
    }
    ++size;
}

// ---------------------------------------------------------------------------
// Ring::clear
// ---------------------------------------------------------------------------
void Ring::clear() {
    if (!head) return;
    Vertex* cur = head;
    do {
        Vertex* nxt = cur->next;
        delete cur;
        cur = nxt;
    } while (cur != head);
    head = nullptr;
    size = 0;
}

// ---------------------------------------------------------------------------
// Ring::signed_area  —  shoelace formula
// ---------------------------------------------------------------------------
double Ring::signed_area() const {
    if (!head || size < 3) return 0.0;
    double area = 0.0;
    Vertex* cur = head;
    do {
        area += cur->x * cur->next->y - cur->next->x * cur->y;
        cur = cur->next;
    } while (cur != head);
    return area * 0.5;
}

// ---------------------------------------------------------------------------
// Ring::compute_candidate
//
// Roles: A = B->prev, B = B, C = B->next, D = B->next->next
//
// APSC formula: find E on the area-preserving constraint line
//   area(AED) = S_quad = area(A,B,C,D)
// that minimises the areal displacement |tri(A,B,E)| + |tri(C,D,E)|.
//
// The constraint is linear in E, defining a line in 2D.
// Both |tri(A,B,E)| and |tri(C,D,E)| are |linear in t| along that line,
// so the optimal t is at one of the two sign-change breakpoints.
// ---------------------------------------------------------------------------
bool Ring::compute_candidate(Vertex* B) {
    assert(B && B->ring_id == id);

    // A collapse removes 2 vertices (B and C) and inserts 1 (E), net -1.
    // We must never let a ring fall below 3 vertices (a triangle is the
    // minimum valid simple polygon).
    if (size <= 3) {
        B->areal_displacement = 1e18;
        B->invalid = true;
        return false;
    }

    Vertex* A = B->prev;
    Vertex* C = B->next;
    Vertex* D = C->next;

    // ------------------------------------------------------------------
    // APSC formula (Kronenfeld et al. 2020)
    //
    // We want E on segment A→D such that the signed area of triangle AED
    // equals the signed area of quadrilateral ABCD.
    //
    // Signed area of quad ABCD (shoelace over 4 vertices):
    //   S = tri_signed_area(A,B,C) + tri_signed_area(A,C,D)
    //
    // Signed area of triangle AED with E = A + t*(D-A):
    //   area(AED) = 0.5 * cross(AD, AE)
    //             = 0.5 * |AD|^2 * t  ... no, let's derive carefully.
    //
    // Let AD = (dx, dy) = (D-A).
    // E = A + t*AD  =>  AE = t*AD
    // area(AED) = 0.5 * (AD x AE) ... but AE = t*AD so cross = 0.
    //
    // Correct derivation: E is NOT constrained to lie on AD.
    // Instead E lies on the LINE through A and D (extended if needed),
    // i.e. E = A + t*(D - A) for any real t.
    //
    // area(triangle A, E, D):
    //   = 0.5 * ((E-A) x (D-A))
    //   = 0.5 * (t*(D-A)) x (D-A)
    //   = 0  ...  that's still zero because E-A is parallel to D-A.
    //
    // The paper places E on the PERPENDICULAR bisector of AD, not on AD.
    // Specifically (Kronenfeld eq. 4):
    //
    //   E = midpoint(A,D) + h * perp_unit(AD)
    //
    // where perp_unit(AD) is the unit vector perpendicular to AD, and h is
    // chosen so area(AED) == S_quad.
    //
    //   area(AED) = 0.5 * |AD| * |h|  (base |AD|, height |h|)
    //   => h = 2 * S_quad / |AD|
    //   (sign of h determines which side of AD E falls on)
    //
    // Areal displacement = area of symmetric difference between
    // polyline A→B→C→D and A→E→D, which equals
    //   |area(ABE)| + |area(BCE)| + |area(CDE)|
    // but more simply: since area is preserved,
    //   displacement = |area(ABE) + area(BCD_new)| summed over left/right
    // The cleanest formula: compute absolute areas of sub-triangles of the
    // symmetric difference.  For a convex quad the symmetric difference has
    // two triangular regions L and R (see Figure 2 of the spec).
    //
    //   disp = |area(A,B,E)| + |area(C,D,E)|
    //          (the two "leftover" triangles when the quad is split by AE and ED)
    // ------------------------------------------------------------------

    double S_quad = tri_signed_area(A->x, A->y, B->x, B->y, C->x, C->y)
                  + tri_signed_area(A->x, A->y, C->x, C->y, D->x, D->y);

    // The constraint for area preservation: area(AED) = S_quad
    // Expanded: 0.5 * ((ex - A.x)*(D.y - A.y) - (ey - A.y)*(D.x - A.x)) = S_quad
    // => p*ex + q*ey = RHS   where p = D.y-A.y, q = -(D.x-A.x)
    double p   = D->y - A->y;
    double q   = -(D->x - A->x);
    double RHS = 2.0 * S_quad + p * A->x + q * A->y;
    double pq2 = p*p + q*q;   // = |AD|^2

    if (pq2 < 1e-15) {
        B->areal_displacement = 1e18;
        B->invalid = true;
        return false;
    }

    // Parameterize all valid E points on the constraint line:
    //   E(t) = E0 + t * (tx, ty)
    // where E0 = foot from origin to constraint line, tangent = (q,-p)/|pq|
    double pq_len = std::sqrt(pq2);
    double tx = q / pq_len;
    double ty = -p / pq_len;
    double t0 = RHS / pq2;
    double e0x = t0 * p;
    double e0y = t0 * q;

    // f(t) = |tri(A,B,E(t))| + |tri(C,D,E(t))|
    // Each tri is linear in t: tri(A,B,E(t)) = c1_0 + k1*t
    // f(t) is piecewise linear; minimum is at one of the two breakpoints
    // (where each term changes sign) or at t=0.
    double k1 = 0.5 * ((B->x - A->x) * ty - (B->y - A->y) * tx);
    double k2 = 0.5 * ((D->x - C->x) * ty - (D->y - C->y) * tx);
    double c1  = tri_signed_area(A->x, A->y, B->x, B->y, e0x, e0y);
    double c2  = tri_signed_area(C->x, C->y, D->x, D->y, e0x, e0y);

    auto f = [&](double t) {
        return std::abs(c1 + k1*t) + std::abs(c2 + k2*t);
    };

    // Breakpoints where each |linear| term changes sign
    double best_t = 0.0;
    double best_f = f(0.0);

    if (std::abs(k1) > 1e-12) {
        double tb = -c1 / k1;
        double fb = f(tb);
        if (fb < best_f) { best_f = fb; best_t = tb; }
    }
    if (std::abs(k2) > 1e-12) {
        double tb = -c2 / k2;
        double fb = f(tb);
        if (fb < best_f) { best_f = fb; best_t = tb; }
    }

    B->ex = e0x + best_t * tx;
    B->ey = e0y + best_t * ty;
    B->areal_displacement = best_f;
    return true;
}

// ---------------------------------------------------------------------------
// Ring::apply_collapse
//
// Before calling this you must have:
//   1. Verified B->invalid == false
//   2. Passed the topology check for the new edges A→E and E→D
//
// What this does:
//   - Creates new vertex E at (B->ex, B->ey)
//   - Inserts E between A and D in the linked list
//   - Removes (deletes) B and C
//   - Updates ring size
//   - Returns pointer to E so the caller can re-enqueue affected candidates
// ---------------------------------------------------------------------------
Vertex* Ring::apply_collapse(Vertex* B) {
    Vertex* A = B->prev;
    Vertex* C = B->next;
    Vertex* D = C->next;

    // Build replacement vertex E
    Vertex* E   = new Vertex();
    E->x        = B->ex;
    E->y        = B->ey;
    E->ring_id  = id;
    E->pq_index = -1;

    // Relink: A ↔ E ↔ D
    A->next = E;
    E->prev = A;
    E->next = D;
    D->prev = E;

    // If head pointed to B or C, redirect it
    if (head == B || head == C) head = E;

    // Free the two removed vertices (pq_index already cleared by caller)
    B->pq_index = -1;
    C->pq_index = -1;
    delete B;
    delete C;

    size -= 1; // removed 2 (B, C), inserted 1 (E)  →  net -1

    return E;
}

// ---------------------------------------------------------------------------
// Ring::collect_vertices
// ---------------------------------------------------------------------------
void Ring::collect_vertices(std::vector<Vertex*>& out) const {
    if (!head) return;
    Vertex* cur = head;
    int vid = 0;
    do {
        cur->vertex_id = vid++;
        out.push_back(cur);
        cur = cur->next;
    } while (cur != head);
}
