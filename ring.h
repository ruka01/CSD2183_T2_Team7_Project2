#pragma once
#include <cmath>
#include <vector>

// ---------------------------------------------------------------------------
// Vertex  — one node in a ring's circular doubly-linked list.
//
// Each vertex "owns" one candidate collapse: the sequence
//   prev → this → next → next->next   (roles: A → B → C → D)
// where *this* is B (the first vertex removed in a collapse).
//
// pq_index lets the priority queue do O(log n) decrease-key by storing
// where in the heap array this vertex currently lives.  Set to -1 when
// the vertex is not in the queue (e.g. it has been removed).
// ---------------------------------------------------------------------------
struct Vertex {
    double x, y;

    Vertex* prev = nullptr;
    Vertex* next = nullptr;

    int ring_id   = 0;
    int vertex_id = 0;   // renumbered on output, not fixed during simplification

    // Priority-queue bookkeeping
    int    pq_index          = -1;   // position in heap array; -1 = not in PQ
    double areal_displacement = 1e18; // cost of collapsing A→B→C→D

    // Position of the replacement vertex E (computed when cost is evaluated)
    double ex = 0.0, ey = 0.0;

    // True if this candidate was evaluated and found topologically invalid.
    // We use lazy invalidation: keep it in the PQ but skip it when popped.
    bool invalid = false;
};

// ---------------------------------------------------------------------------
// Ring  — a closed polygon ring stored as a circular doubly-linked list.
//
// Ownership: the Ring allocates all its Vertex nodes and frees them in the
// destructor.  Do not delete Vertex pointers from outside the Ring.
// ---------------------------------------------------------------------------
class Ring {
public:
    int  id;          // 0 = exterior, 1+ = interior holes
    bool is_exterior; // true iff id == 0

    // Pointer to any one vertex in the ring (used as the entry point).
    // After collapses this still points to a live vertex.
    Vertex* head = nullptr;

    int size = 0;     // current vertex count (updated on every insert/remove)

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    explicit Ring(int ring_id)
        : id(ring_id), is_exterior(ring_id == 0) {}

    ~Ring() { clear(); }

    // Append a new vertex (x, y) at the tail of the circular list.
    // Call this in input order; after all appends the list is automatically
    // circular (last->next == head).
    void append(double x, double y);

    // -----------------------------------------------------------------------
    // Geometry
    // -----------------------------------------------------------------------

    // Signed area via the shoelace formula.
    // Positive for CCW rings (exterior), negative for CW rings (holes).
    double signed_area() const;

    // -----------------------------------------------------------------------
    // Collapse operation
    // -----------------------------------------------------------------------

    // Given that vertex B is the candidate to collapse (roles A→B→C→D),
    // compute the position of area-preserving replacement vertex E and store
    // it in B->ex, B->ey.  Also compute and store B->areal_displacement.
    //
    // Returns false if the collapse is geometrically degenerate (e.g. A == D).
    bool compute_candidate(Vertex* B);

    // Apply the collapse: remove B and C from the linked list and insert a
    // new vertex E at their position.  Returns a pointer to the new vertex E
    // so the caller can update the priority queue.
    //
    // Preconditions:
    //   - B->invalid == false
    //   - topology check has already passed
    Vertex* apply_collapse(Vertex* B);

    // -----------------------------------------------------------------------
    // Output helpers
    // -----------------------------------------------------------------------

    // Write all vertices of this ring into the output vector in linked-list
    // order, renumbering vertex_id from 0.
    void collect_vertices(std::vector<Vertex*>& out) const;

private:
    void clear(); // free all allocated Vertex nodes
};
