#include "io.h"
#include "ring.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <map>
#include <vector>
#include <cstdio>
#include <cmath>

// ---------------------------------------------------------------------------
// read_csv
// ---------------------------------------------------------------------------
std::vector<Ring*> read_csv(const std::string& filename) {
    std::ifstream f(filename);
    if (!f.is_open())
        throw std::runtime_error("Cannot open input file: " + filename);

    std::string line;
    // Skip header line
    if (!std::getline(f, line))
        throw std::runtime_error("Empty input file");

    // Use a map so rings are built in ring_id order regardless of CSV row order
    std::map<int, Ring*> ring_map;

    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string tok;

        // Parse: ring_id,vertex_id,x,y
        int ring_id, vertex_id;
        double x, y;

        try {
            std::getline(ss, tok, ','); ring_id   = std::stoi(tok);
            std::getline(ss, tok, ','); vertex_id = std::stoi(tok);
            std::getline(ss, tok, ','); x         = std::stod(tok);
            std::getline(ss, tok, ','); y         = std::stod(tok);
        } catch (...) {
            throw std::runtime_error("Malformed CSV line: " + line);
        }

        if (ring_map.find(ring_id) == ring_map.end())
            ring_map[ring_id] = new Ring(ring_id);

        ring_map[ring_id]->append(x, y);
        (void)vertex_id; // input vertex_id not stored; we renumber on output
    }

    // Convert map to vector in ring_id order (0, 1, 2, ...)
    std::vector<Ring*> rings;
    rings.reserve(ring_map.size());
    for (auto& [id, ring] : ring_map)
        rings.push_back(ring);

    if (rings.empty())
        throw std::runtime_error("No rings found in input");

    return rings;
}

// ---------------------------------------------------------------------------
// write_csv
// ---------------------------------------------------------------------------
void write_csv(const std::vector<Ring*>& rings, double input_signed_area,
               double total_areal_displacement,
               const std::string& output_file) {
    FILE* f = fopen(output_file.c_str(), "w");
    if (!f)
        throw std::runtime_error("Cannot open output file: " + output_file);

    fprintf(f, "ring_id,vertex_id,x,y\n");

    double output_signed_area = 0.0;

    for (Ring* ring : rings) {
        std::vector<Vertex*> verts;
        ring->collect_vertices(verts);

        for (Vertex* v : verts)
            fprintf(f, "%d,%d,%g,%g\n", ring->id, v->vertex_id, v->x, v->y);

        output_signed_area += ring->signed_area();
    }

    fclose(f);

    // Three summary lines to stdout (scientific notation, 6 decimal places)
    printf("Total signed area in input: %.6e\n",  input_signed_area);
    printf("Total signed area in output: %.6e\n", output_signed_area);
    printf("Total areal displacement: %.6e\n",    total_areal_displacement);
}
