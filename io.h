#pragma once
#include <string>
#include <vector>

class Ring;

// ---------------------------------------------------------------------------
// CSV I/O
// ---------------------------------------------------------------------------

// Parse the input CSV file.  Returns a vector of Ring* in ring_id order.
// Ring 0 is always the exterior ring.
// Throws std::runtime_error on malformed input.
std::vector<Ring*> read_csv(const std::string& filename);

// Write the simplified polygon to `output_file` in CSV format.
// Also prints the three summary lines to stdout.
//
// `input_signed_area`  — total signed area of the input (computed before
//                        any simplification so it stays fixed).
// `rings`              — the simplified rings (in ring_id order).
void write_csv(const std::vector<Ring*>& rings, double input_signed_area,
               double total_areal_displacement,
               const std::string& output_file);
