#include "ring.h"
#include "priority_queue.h"
#include "topology.h"
#include "io.h"
#include "view_template.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <filesystem>

// ---------------------------------------------------------------------------
// total_vertex_count
// ---------------------------------------------------------------------------
static int total_vertex_count(const std::vector<Ring*>& rings) {
    int n = 0;
    for (Ring* r : rings) n += r->size;
    return n;
}

// ---------------------------------------------------------------------------
// init_queue
// ---------------------------------------------------------------------------
static void init_queue(const std::vector<Ring*>& rings, CollapseQueue& pq) {
    for (Ring* ring : rings) {
        if (ring->size < 4) continue;
        Vertex* cur = ring->head;
        do {
            if (ring->compute_candidate(cur)) pq.push(cur);
            cur = cur->next;
        } while (cur != ring->head);
    }
}

// ---------------------------------------------------------------------------
// update_neighbours
// ---------------------------------------------------------------------------
static void update_neighbours(Vertex* E, const std::vector<Ring*>& rings,
    CollapseQueue& pq) {
    Ring* ring = nullptr;
    for (Ring* r : rings)
        if (r->id == E->ring_id) { ring = r; break; }
    if (!ring) return;

    Vertex* raw[4] = { E, E->prev, E->prev->prev, E->next };
    for (int i = 0; i < 4; ++i) {
        Vertex* v = raw[i];
        if (!v || v->ring_id != ring->id) continue;
        bool seen = false;
        for (int j = 0; j < i; ++j) if (raw[j] == v) { seen = true; break; }
        if (seen) continue;

        if (ring->size <= 3) {
            v->invalid = true;
            if (v->pq_index >= 0) { v->areal_displacement = 1e18; pq.update(v); }
            continue;
        }
        bool ok = ring->compute_candidate(v);
        if (!ok) { v->invalid = true; if (v->pq_index >= 0) pq.update(v); continue; }
        if (v->pq_index >= 0) pq.update(v); else pq.push(v);
    }
}

// ---------------------------------------------------------------------------
// simplify
// ---------------------------------------------------------------------------
static double simplify(std::vector<Ring*>& rings, int target_n) {
    CollapseQueue pq;
    init_queue(rings, pq);
    double total = 0.0;

    while (total_vertex_count(rings) > target_n && !pq.empty()) {
        Vertex* B = pq.pop_best();
        if (!B) break;

        Ring* ring = nullptr;
        for (Ring* r : rings) if (r->id == B->ring_id) { ring = r; break; }
        if (!ring || ring->size <= 3) continue;

        if (!topology_valid(B, rings)) { B->invalid = true; continue; }

        total += B->areal_displacement;

        Vertex* C = B->next;
        B->invalid = C->invalid = true;
        if (C->pq_index >= 0) { C->areal_displacement = 1e18; pq.update(C); }

        update_neighbours(ring->apply_collapse(B), rings, pq);
    }
    return total;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: ./simplify <input_file.csv> <target_vertices>\n";
        return 1;
    }

    std::string input_file = argv[1];
    int         target_verts = std::stoi(argv[2]);

    std::vector<Ring*> rings;
    try { rings = read_csv(input_file); }
    catch (const std::exception& e) {
        std::cerr << "Error reading input: " << e.what() << "\n"; return 1;
    }

    double input_signed_area = 0.0;
    for (Ring* r : rings) input_signed_area += r->signed_area();

    int current_n = total_vertex_count(rings);
    std::cerr << "Input: " << current_n << " vertices across "
        << rings.size() << " ring(s)\n";

    double total_areal_displacement = 0.0;
    if (current_n > target_verts)
        total_areal_displacement = simplify(rings, target_verts);

    std::cerr << "Output: " << total_vertex_count(rings) << " vertices\n";

    // -----------------------------------------------------------------------
    // Create output folder named after the input file stem
    // -----------------------------------------------------------------------
    namespace fs = std::filesystem;

    fs::path    input_path(input_file);
    std::string folder_name = input_path.stem().string();
    fs::path    folder(folder_name);
    fs::create_directories(folder);

    // Copy input CSV into folder
    fs::copy_file(input_path, folder / "input.csv",
        fs::copy_options::overwrite_existing);

    // Write output CSV into folder
    std::string dest_output = (folder / "output.csv").string();
    write_csv(rings, input_signed_area, total_areal_displacement, dest_output);
    std::cerr << "Saved output to " << dest_output << "\n";

    // -----------------------------------------------------------------------
    // Read both CSVs back as strings, escape for JS template literals,
    // then substitute into the HTML template and write View.html.
    // The file is fully self-contained — no fetch() needed, works from file://.
    // -----------------------------------------------------------------------
    auto read_str = [](const std::string& path) -> std::string {
        std::ifstream f(path);
        std::ostringstream ss;
        ss << f.rdbuf();
        return ss.str();
        };

    auto escape_js = [](const std::string& s) -> std::string {
        std::string out;
        out.reserve(s.size());
        for (char c : s) {
            if (c == '\\') out += "\\\\";
            else if (c == '`')  out += "\\`";
            else                out += c;
        }
        return out;
        };

    std::string input_data = escape_js(read_str((folder / "input.csv").string()));
    std::string output_data = escape_js(read_str(dest_output));

    // Replace the two markers in the template
    std::string html = VIEW_HTML_TEMPLATE;
    auto replace_marker = [&](const std::string& marker, const std::string& val) {
        size_t pos = html.find(marker);
        if (pos != std::string::npos)
            html.replace(pos, marker.size(), val);
        };
    replace_marker("%%INPUT_CSV%%", input_data);
    replace_marker("%%OUTPUT_CSV%%", output_data);

    std::string view_path = (folder / "View.html").string();
    std::ofstream out(view_path);
    if (!out.is_open())
        std::cerr << "Warning: could not write View.html\n";
    else {
        out << html;
        std::cerr << "Wrote View.html to " << view_path << "\n";
    }

    std::cerr << "\nOpen " << folder_name << "/View.html in your browser.\n";

    for (Ring* r : rings) delete r;
    return 0;
}