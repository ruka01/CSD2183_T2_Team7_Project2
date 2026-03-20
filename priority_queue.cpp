#include "priority_queue.h"
#include "ring.h"
#include <cassert>
#include <stdexcept>

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

bool CollapseQueue::less_than(int i, int j) const {
    return heap_[i]->areal_displacement < heap_[j]->areal_displacement;
}

void CollapseQueue::swap_entries(int i, int j) {
    std::swap(heap_[i], heap_[j]);
    heap_[i]->pq_index = i;
    heap_[j]->pq_index = j;
}

// ---------------------------------------------------------------------------
// sift_up — used after push or after cost decreases
// ---------------------------------------------------------------------------
void CollapseQueue::sift_up(int i) {
    while (i > 0) {
        int p = parent(i);
        if (less_than(i, p)) {
            swap_entries(i, p);
            i = p;
        } else {
            break;
        }
    }
}

// ---------------------------------------------------------------------------
// sift_down — used after pop or after cost increases
// ---------------------------------------------------------------------------
void CollapseQueue::sift_down(int i) {
    int n = static_cast<int>(heap_.size());
    while (true) {
        int smallest = i;
        int l = left(i), r = right(i);
        if (l < n && less_than(l, smallest)) smallest = l;
        if (r < n && less_than(r, smallest)) smallest = r;
        if (smallest == i) break;
        swap_entries(i, smallest);
        i = smallest;
    }
}

// ---------------------------------------------------------------------------
// push
// ---------------------------------------------------------------------------
void CollapseQueue::push(Vertex* v) {
    assert(v->pq_index == -1 && "Vertex already in queue");
    v->pq_index = static_cast<int>(heap_.size());
    heap_.push_back(v);
    sift_up(v->pq_index);
}

// ---------------------------------------------------------------------------
// update  —  call after v->areal_displacement has been recomputed
// ---------------------------------------------------------------------------
void CollapseQueue::update(Vertex* v) {
    assert(v->pq_index >= 0 && "Vertex not in queue");
    int i = v->pq_index;
    // Try both directions; only one will actually move the node
    sift_up(i);
    sift_down(v->pq_index); // pq_index may have changed after sift_up
}

// ---------------------------------------------------------------------------
// pop_best
//
// Skips entries that are marked invalid or have cost == 1e18.
// Lazy deletion: invalid entries stay in the heap until they bubble to
// the top, then get discarded here.
// ---------------------------------------------------------------------------
Vertex* CollapseQueue::pop_best() {
    while (!heap_.empty()) {
        // Move root to back, shrink heap
        int last = static_cast<int>(heap_.size()) - 1;
        swap_entries(0, last);
        Vertex* v = heap_.back();
        heap_.pop_back();
        v->pq_index = -1;

        if (!heap_.empty()) sift_down(0);

        // Skip dead entries
        if (v->invalid || v->areal_displacement >= 1e17) continue;

        return v;
    }
    return nullptr; // queue exhausted
}
