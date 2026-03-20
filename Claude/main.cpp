#include "ring.h"
#include "priority_queue.h"
#include "topology.h"
#include "io.h"
#include <iostream>
#include <vector>
#include <numeric>

// ---------------------------------------------------------------------------
// total_vertex_count  —  sum of all ring sizes
// ---------------------------------------------------------------------------
static int total_vertex_count(const std::vector<Ring*>& rings) {
    int n = 0;
    for (Ring* r : rings) n += r->size;
    return n;
}

// ---------------------------------------------------------------------------
// init_queue
//
// For every ring, compute the candidate cost for every vertex B (treating
// B as the first interior vertex of A→B→C→D) and push B into the PQ.
//
// Vertices in rings with fewer than 4 vertices cannot be collapsed and are
// not added to the queue.
// ---------------------------------------------------------------------------
static void init_queue(const std::vector<Ring*>& rings, CollapseQueue& pq) {
    for (Ring* ring : rings) {
        if (ring->size < 4) continue;
        Vertex* cur = ring->head;
        do {
            // cur plays the role of B
            if (ring->compute_candidate(cur))
                pq.push(cur);
            cur = cur->next;
        } while (cur != ring->head);
    }
}

// ---------------------------------------------------------------------------
// update_neighbours
//
// After collapsing B (which produced vertex E), up to 4 candidate entries
// in the PQ become stale.  Re-evaluate and update them:
//
//   The 4 affected B-roles after inserting E between A and D:
//     1. E itself     (sequence A → E → D → D->next)
//     2. A            (sequence A->prev → A → E → D)
//     3. A->prev      (sequence A->prev->prev → A->prev → A → E)
//     4. D            (sequence E → D → D->next → D->next->next)
//
//   Additionally, some old entries (for the deleted vertices B and C) are
//   already gone; their pq_index was set to -1 inside apply_collapse via
//   the delete operation, but the PQ may still hold a pointer to freed
//   memory if we used lazy deletion.  To avoid this, we must mark them
//   invalid BEFORE calling apply_collapse (see main loop below).
// ---------------------------------------------------------------------------
static void update_neighbours(Vertex* E, const std::vector<Ring*>& rings,
                               CollapseQueue& pq) {
    Ring* ring = nullptr;
    for (Ring* r : rings)
        if (r->id == E->ring_id) { ring = r; break; }
    if (!ring) return;

    // The four affected B-positions after inserting E between A and D:
    //   E itself, A (=E->prev), A->prev, D (=E->next)
    // On small rings these can alias each other — use a set to deduplicate.
    Vertex* raw[4] = {
        E,
        E->prev,
        E->prev->prev,
        E->next
    };

    // Deduplicate: on a 3-vertex ring, E->prev->prev == E->next etc.
    for (int i = 0; i < 4; ++i) {
        Vertex* v = raw[i];
        if (!v || v->ring_id != ring->id) continue;

        // Skip if this pointer appeared earlier in the array
        bool seen = false;
        for (int j = 0; j < i; ++j)
            if (raw[j] == v) { seen = true; break; }
        if (seen) continue;

        // If ring is too small, invalidate and skip
        if (ring->size <= 3) {
            v->invalid = true;
            if (v->pq_index >= 0) {
                v->areal_displacement = 1e18;
                pq.update(v);
            }
            continue;
        }

        bool ok = ring->compute_candidate(v);
        if (!ok) {
            v->invalid = true;
            if (v->pq_index >= 0) pq.update(v);
            continue;
        }

        if (v->pq_index >= 0) {
            pq.update(v);
        } else {
            pq.push(v);
        }
    }
}

// ---------------------------------------------------------------------------
// simplify  —  the main APSC loop
// ---------------------------------------------------------------------------
static double simplify(std::vector<Ring*>& rings, int target_n) {
    CollapseQueue pq;
    init_queue(rings, pq);

    double total_areal_displacement = 0.0;

    while (total_vertex_count(rings) > target_n && !pq.empty()) {
        Vertex* B = pq.pop_best();
        if (!B) break;

        // Re-validate: the ring may have shrunk since B was enqueued
        Ring* ring = nullptr;
        for (Ring* r : rings)
            if (r->id == B->ring_id) { ring = r; break; }
        if (!ring || ring->size <= 3) continue;

        // Topology check for the proposed collapse
        if (!topology_valid(B, rings)) {
            B->invalid = true;
            // Do not re-enqueue — mark and move on.
            // The affected neighbours will be re-evaluated from other collapses.
            continue;
        }

        // Accumulate areal displacement before the vertices are deleted
        total_areal_displacement += B->areal_displacement;

        // Grab C before apply_collapse frees it, and remove it from PQ
        Vertex* C = B->next;
        B->invalid = true;
        C->invalid = true;
        // Clear C's pq_index so the heap slot is treated as stale
        if (C->pq_index >= 0) {
            // Mark the heap slot invalid; pop_best will skip it
            // We can't physically remove from mid-heap without decrease-key,
            // so set cost to sentinel — pop_best skips invalid entries anyway
            C->areal_displacement = 1e18;
            pq.update(C);
            // After update C floats toward the bottom; mark invalid so pop skips it
        }

        // Apply the collapse: returns the new vertex E
        Vertex* E = ring->apply_collapse(B);

        // Update PQ for the four affected neighbours
        update_neighbours(E, rings, pq);
    }

    return total_areal_displacement;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: ./simplify <input_file.csv> <target_vertices>\n";
        return 1;
    }

    std::string input_file   = argv[1];
    int         target_verts = std::stoi(argv[2]);

    // Read input
    std::vector<Ring*> rings;
    try {
        rings = read_csv(input_file);
    } catch (const std::exception& e) {
        std::cerr << "Error reading input: " << e.what() << "\n";
        return 1;
    }

    // Record input area before any changes
    double input_signed_area = 0.0;
    for (Ring* r : rings) input_signed_area += r->signed_area();

    int current_n = total_vertex_count(rings);
    std::cerr << "Input: " << current_n << " vertices across "
              << rings.size() << " ring(s)\n";

    // Run simplification if needed
    double total_areal_displacement = 0.0;
    if (current_n > target_verts)
        total_areal_displacement = simplify(rings, target_verts);

    std::cerr << "Output: " << total_vertex_count(rings) << " vertices\n";

    // Write result to after.csv
    const std::string output_file = "after.csv";
    write_csv(rings, input_signed_area, total_areal_displacement, output_file);
    std::cerr << "Saved simplified polygon to " << output_file << "\n";

    // Clean up
    for (Ring* r : rings) delete r;

    return 0;
}
