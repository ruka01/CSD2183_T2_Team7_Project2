#pragma once
#include <vector>
#include <functional>

// ---------------------------------------------------------------------------
// CollapseQueue  —  min-heap keyed on Vertex::areal_displacement.
//
// Design notes:
//   - We store raw Vertex* pointers in the heap array.
//   - Each Vertex tracks its own position in the heap via pq_index.
//     This lets us do O(log n) decrease-key (re-cost after a neighbour
//     is updated) instead of rebuilding the heap from scratch.
//   - Vertices with areal_displacement == 1e18 or invalid == true are
//     "dead" entries; they are skipped when popped.
//
// API:
//   push(v)        — insert v (v->areal_displacement must already be set)
//   update(v)      — call after recomputing v->areal_displacement; fixes heap
//   pop_best()     — remove and return the vertex with lowest valid cost,
//                    or nullptr if the queue is empty / all entries invalid
//   empty()        — true if no valid entries remain
// ---------------------------------------------------------------------------

struct Vertex; // forward declaration (defined in ring.h)

class CollapseQueue {
public:
    CollapseQueue() = default;

    // Insert a vertex that is not yet in the queue.
    void push(Vertex* v);

    // Re-heapify after v->areal_displacement has changed.
    // v must already be in the queue (v->pq_index != -1).
    void update(Vertex* v);

    // Remove and return the cheapest valid candidate.
    // Returns nullptr when the queue is exhausted.
    Vertex* pop_best();

    bool empty() const { return heap_.empty(); }
    int  size()  const { return static_cast<int>(heap_.size()); }

private:
    std::vector<Vertex*> heap_;

    // Standard binary-heap helpers
    // Parent/child index arithmetic (0-based)
    static int parent(int i) { return (i - 1) / 2; }
    static int left(int i)   { return 2 * i + 1;   }
    static int right(int i)  { return 2 * i + 2;   }

    // Move node at index i upward until heap property restored
    void sift_up(int i);

    // Move node at index i downward until heap property restored
    void sift_down(int i);

    // Swap two heap entries, keeping pq_index consistent
    void swap_entries(int i, int j);

    bool less_than(int i, int j) const; // true if heap_[i] < heap_[j] by cost
};
