#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

static constexpr double EPS = 1e-9;

struct Point
{
    double x{};
    double y{};
};

struct Ring
{
    int ring_id{};
    std::vector<Point> pts;
};

struct Polygon
{
    std::vector<Ring> rings;
};

struct CollapseCandidate
{
    bool valid{ false };
    int ring_idx{ -1 };
    int start_idx{ -1 }; // A index in A,B,C,D
    Point e{};
    double displacement{ std::numeric_limits<double>::infinity() };
};

double cross(const Point& a, const Point& b)
{
    return a.x * b.y - a.y * b.x;
}

Point operator+(const Point& a, const Point& b)
{
    return { a.x + b.x, a.y + b.y };
}

Point operator-(const Point& a, const Point& b)
{
    return { a.x - b.x, a.y - b.y };
}

Point operator*(const Point& a, double s)
{
    return { a.x * s, a.y * s };
}

double dot(const Point& a, const Point& b)
{
    return a.x * b.x + a.y * b.y;
}

double norm2(const Point& a)
{
    return dot(a, a);
}

double orient(const Point& a, const Point& b, const Point& c)
{
    return cross(b - a, c - a);
}

bool almostEqual(double a, double b, double eps = EPS)
{
    return std::fabs(a - b) <= eps;
}

bool samePoint(const Point& a, const Point& b, double eps = EPS)
{
    return almostEqual(a.x, b.x, eps) && almostEqual(a.y, b.y, eps);
}

double signedArea(const std::vector<Point>& poly)
{
    const int n = static_cast<int>(poly.size());
    if (n < 3)
        return 0.0;

    double s = 0.0;
    for (int i = 0; i < n; ++i)
    {
        const Point& p = poly[i];
        const Point& q = poly[(i + 1) % n];
        s += cross(p, q);
    }
    return 0.5 * s;
}

double totalSignedArea(const Polygon& poly)
{
    double s = 0.0;
    for (const Ring& r : poly.rings)
        s += signedArea(r.pts);
    return s;
}

int totalVertices(const Polygon& poly)
{
    int cnt = 0;
    for (const Ring& r : poly.rings)
        cnt += static_cast<int>(r.pts.size());
    return cnt;
}

bool onSegment(const Point& a, const Point& b, const Point& p)
{
    if (std::fabs(orient(a, b, p)) > EPS)
        return false;

    return p.x >= std::min(a.x, b.x) - EPS &&
           p.x <= std::max(a.x, b.x) + EPS &&
           p.y >= std::min(a.y, b.y) - EPS &&
           p.y <= std::max(a.y, b.y) + EPS;
}

bool segmentsIntersectInclusive(const Point& a, const Point& b,
                                const Point& c, const Point& d)
{
    const double o1 = orient(a, b, c);
    const double o2 = orient(a, b, d);
    const double o3 = orient(c, d, a);
    const double o4 = orient(c, d, b);

    if (((o1 > EPS && o2 < -EPS) || (o1 < -EPS && o2 > EPS)) &&
        ((o3 > EPS && o4 < -EPS) || (o3 < -EPS && o4 > EPS)))
        return true;

    if (std::fabs(o1) <= EPS && onSegment(a, b, c)) return true;
    if (std::fabs(o2) <= EPS && onSegment(a, b, d)) return true;
    if (std::fabs(o3) <= EPS && onSegment(c, d, a)) return true;
    if (std::fabs(o4) <= EPS && onSegment(c, d, b)) return true;

    return false;
}

bool segmentsProperlyIntersectOrOverlap(const Point& a, const Point& b,
                                        const Point& c, const Point& d)
{
    if (!segmentsIntersectInclusive(a, b, c, d))
        return false;

    // Allow touching only at shared endpoints
    const bool shareEndpoint =
        samePoint(a, c) || samePoint(a, d) || samePoint(b, c) || samePoint(b, d);

    if (shareEndpoint)
    {
        // If they are collinear and overlap more than endpoint-touch, reject.
        if (std::fabs(orient(a, b, c)) <= EPS && std::fabs(orient(a, b, d)) <= EPS)
        {
            // Collinear: conservative rejection unless just single-endpoint touch.
            return true;
        }
        return false;
    }

    return true;
}

bool pointInRing(const std::vector<Point>& ring, const Point& p)
{
    // Ray casting
    bool inside = false;
    const int n = static_cast<int>(ring.size());
    for (int i = 0, j = n - 1; i < n; j = i++)
    {
        const Point& a = ring[j];
        const Point& b = ring[i];

        const bool intersect = ((a.y > p.y) != (b.y > p.y)) &&
                               (p.x < (b.x - a.x) * (p.y - a.y) / ((b.y - a.y) + 1e-30) + a.x);
        if (intersect)
            inside = !inside;
    }
    return inside;
}

Polygon readCSV(const std::string& filename)
{
    std::ifstream fin(filename);
    if (!fin)
    {
        throw std::runtime_error("Failed to open input file: " + filename);
    }

    std::string line;
    std::getline(fin, line); // header

    std::map<int, Ring> ringMap;

    while (std::getline(fin, line))
    {
        if (line.empty())
            continue;

        std::stringstream ss(line);
        std::string tok;
        std::vector<std::string> fields;

        while (std::getline(ss, tok, ','))
            fields.push_back(tok);

        if (fields.size() != 4)
            continue;

        int ring_id = std::stoi(fields[0]);
        // int vertex_id = std::stoi(fields[1]); // not needed
        double x = std::stod(fields[2]);
        double y = std::stod(fields[3]);

        ringMap[ring_id].ring_id = ring_id;
        ringMap[ring_id].pts.push_back({ x, y });
    }

    Polygon poly;
    for (auto& [id, ring] : ringMap)
        poly.rings.push_back(ring);

    return poly;
}

void writeCSVToStream(std::ostream& os, const Polygon& poly)
{
    os << "ring_id,vertex_id,x,y\n";
    for (int r = 0; r < static_cast<int>(poly.rings.size()); ++r)
    {
        const auto& pts = poly.rings[r].pts;
        for (int i = 0; i < static_cast<int>(pts.size()); ++i)
        {
            os << r << "," << i << ","
               << std::setprecision(15) << pts[i].x << ","
               << std::setprecision(15) << pts[i].y << "\n";
        }
    }
}

double triangleAreaAbs(const Point& a, const Point& b, const Point& c)
{
    return 0.5 * std::fabs(orient(a, b, c));
}

double collapseDisplacementApprox(const Point& a, const Point& b, const Point& c,
                                  const Point& d, const Point& e)
{
    // Conservative local estimate of areal displacement.
    // Uses triangulation of the closed chain A-B-C-D-E.
    return triangleAreaAbs(e, a, b) +
           triangleAreaAbs(e, b, c) +
           triangleAreaAbs(e, c, d);
}

Point computeAreaPreservingBasePoint(const Point& a, const Point& b,
                                     const Point& c, const Point& d)
{
    // Preserve signed area of A-B-C-D by replacing B,C with E in A-E-D.
    // Need:
    // cross(A,B) + cross(B,C) + cross(C,D) = cross(A,E) + cross(E,D)
    // => cross(E, D-A) = cross(A,B) + cross(B,C) + cross(C,D)
    const Point v = d - a;
    const double k = cross(a, b) + cross(b, c) + cross(c, d);
    const double vv = norm2(v);

    if (vv <= EPS)
        return a;

    // One solution E0 satisfying cross(E0, v) = k:
    // E0 = (k * vy / |v|^2, -k * vx / |v|^2)
    return { (k * v.y) / vv, (-k * v.x) / vv };
}

Point linePointAt(const Point& base, const Point& dir, double t)
{
    return base + dir * t;
}

bool ringHasSelfIntersection(const std::vector<Point>& pts)
{
    const int n = static_cast<int>(pts.size());
    if (n < 3)
        return true;

    for (int i = 0; i < n; ++i)
    {
        Point a = pts[i];
        Point b = pts[(i + 1) % n];

        for (int j = i + 1; j < n; ++j)
        {
            // Skip same edge and adjacent edges in cyclic polygon
            if (j == i) continue;
            if ((j + 1) % n == i) continue;
            if ((i + 1) % n == j) continue;

            Point c = pts[j];
            Point d = pts[(j + 1) % n];

            if (segmentsProperlyIntersectOrOverlap(a, b, c, d))
                return true;
        }
    }
    return false;
}

bool ringsIntersect(const std::vector<Point>& a, const std::vector<Point>& b)
{
    const int na = static_cast<int>(a.size());
    const int nb = static_cast<int>(b.size());

    for (int i = 0; i < na; ++i)
    {
        Point a1 = a[i];
        Point a2 = a[(i + 1) % na];
        for (int j = 0; j < nb; ++j)
        {
            Point b1 = b[j];
            Point b2 = b[(j + 1) % nb];

            if (segmentsProperlyIntersectOrOverlap(a1, a2, b1, b2))
                return true;
        }
    }

    // Containment check
    if (!a.empty() && pointInRing(b, a[0])) return true;
    if (!b.empty() && pointInRing(a, b[0])) return true;

    return false;
}

bool topologyValid(const Polygon& poly)
{
    // Ring simplicity
    for (const Ring& r : poly.rings)
    {
        if (static_cast<int>(r.pts.size()) < 3)
            return false;
        if (ringHasSelfIntersection(r.pts))
            return false;
    }

    // Pairwise ring non-intersection
    for (int i = 0; i < static_cast<int>(poly.rings.size()); ++i)
    {
        for (int j = i + 1; j < static_cast<int>(poly.rings.size()); ++j)
        {
            if (ringsIntersect(poly.rings[i].pts, poly.rings[j].pts))
                return false;
        }
    }

    // Hole containment heuristic: every hole point must remain inside exterior ring.
    if (poly.rings.empty())
        return false;

    const auto& exterior = poly.rings[0].pts;
    for (int i = 1; i < static_cast<int>(poly.rings.size()); ++i)
    {
        if (poly.rings[i].pts.empty())
            return false;
        if (!pointInRing(exterior, poly.rings[i].pts[0]))
            return false;
    }

    return true;
}

std::vector<Point> applyCollapseToRing(const std::vector<Point>& pts, int startIdx, const Point& e)
{
    const int n = static_cast<int>(pts.size());
    std::vector<Point> out;
    out.reserve(n - 1);

    const int A = startIdx;
    const int B = (A + 1) % n;
    const int C = (A + 2) % n;
    const int D = (A + 3) % n;

    for (int i = 0; i < n; ++i)
    {
        if (i == B || i == C)
            continue;

        if (i == A)
        {
            out.push_back(pts[A]);
            out.push_back(e);
            continue;
        }

        // Skip D here? No, keep D in normal loop.
        out.push_back(pts[i]);
    }

    // The above duplicates ordering incorrectly if A is not 0 in cyclic order.
    // Rebuild in cyclic order from A.
    std::vector<Point> cyc;
    cyc.reserve(n - 1);
    cyc.push_back(pts[A]);
    cyc.push_back(e);
    int i = D;
    while (i != A)
    {
        cyc.push_back(pts[i]);
        i = (i + 1) % n;
    }

    return cyc;
}

bool ringOrientationMatchesOriginal(const std::vector<Point>& original,
                                    const std::vector<Point>& simplified)
{
    const double a0 = signedArea(original);
    const double a1 = signedArea(simplified);
    if (std::fabs(a0) <= EPS || std::fabs(a1) <= EPS)
        return false;
    return (a0 > 0 && a1 > 0) || (a0 < 0 && a1 < 0);
}

CollapseCandidate bestCollapseForRing(const Polygon& poly, int ringIdx)
{
    CollapseCandidate best;
    const auto& pts = poly.rings[ringIdx].pts;
    const int n = static_cast<int>(pts.size());

    if (n <= 3)
        return best;

    for (int i = 0; i < n; ++i)
    {
        const Point& a = pts[i];
        const Point& b = pts[(i + 1) % n];
        const Point& c = pts[(i + 2) % n];
        const Point& d = pts[(i + 3) % n];

        Point v = d - a;
        double vv = norm2(v);
        if (vv <= EPS)
            continue;

        Point dir = { v.x / std::sqrt(vv), v.y / std::sqrt(vv) };
        Point base = computeAreaPreservingBasePoint(a, b, c, d);

        // Search along the area-preserving line E = base + t * dir
        // using a bounded 1D sampling + local refinement.
        const double ta = dot(a, dir);
        const double tb = dot(b, dir);
        const double tc = dot(c, dir);
        const double td = dot(d, dir);

        double lo = std::min(std::min(ta, tb), std::min(tc, td));
        double hi = std::max(std::max(ta, tb), std::max(tc, td));
        const double span = std::max(1.0, hi - lo);

        lo -= 3.0 * span;
        hi += 3.0 * span;

        const int SAMPLES = 80;
        double bestT = 0.0;
        double bestCost = std::numeric_limits<double>::infinity();

        for (int s = 0; s <= SAMPLES; ++s)
        {
            double t = lo + (hi - lo) * static_cast<double>(s) / SAMPLES;
            Point e = linePointAt(base, dir, t);
            double cost = collapseDisplacementApprox(a, b, c, d, e);

            if (cost < bestCost)
            {
                bestCost = cost;
                bestT = t;
            }
        }

        // Local ternary refinement around best sample
        double left = std::max(lo, bestT - span);
        double right = std::min(hi, bestT + span);

        for (int it = 0; it < 60; ++it)
        {
            double m1 = left + (right - left) / 3.0;
            double m2 = right - (right - left) / 3.0;

            Point e1 = linePointAt(base, dir, m1);
            Point e2 = linePointAt(base, dir, m2);

            double c1 = collapseDisplacementApprox(a, b, c, d, e1);
            double c2 = collapseDisplacementApprox(a, b, c, d, e2);

            if (c1 < c2)
                right = m2;
            else
                left = m1;
        }

        Point e = linePointAt(base, dir, 0.5 * (left + right));

        // Avoid nearly duplicate vertices
        if (samePoint(e, a) || samePoint(e, d))
            continue;

        std::vector<Point> newRingPts = applyCollapseToRing(pts, i, e);

        if (static_cast<int>(newRingPts.size()) < 3)
            continue;

        // Area check for this ring
        double oldArea = signedArea(pts);
        double newArea = signedArea(newRingPts);
        if (std::fabs(oldArea - newArea) > 1e-6)
            continue;

        // Orientation check
        if (!ringOrientationMatchesOriginal(pts, newRingPts))
            continue;

        Polygon trial = poly;
        trial.rings[ringIdx].pts = std::move(newRingPts);

        if (!topologyValid(trial))
            continue;

        double disp = collapseDisplacementApprox(a, b, c, d, e);
        if (disp < best.displacement)
        {
            best.valid = true;
            best.ring_idx = ringIdx;
            best.start_idx = i;
            best.e = e;
            best.displacement = disp;
        }
    }

    return best;
}

Polygon simplifyPolygon(const Polygon& input, int targetVertices, double& totalDisplacement)
{
    Polygon poly = input;
    totalDisplacement = 0.0;

    while (totalVertices(poly) > targetVertices)
    {
        CollapseCandidate best;

        for (int r = 0; r < static_cast<int>(poly.rings.size()); ++r)
        {
            CollapseCandidate cand = bestCollapseForRing(poly, r);
            if (cand.valid && cand.displacement < best.displacement)
                best = cand;
        }

        if (!best.valid)
            break;

        auto& pts = poly.rings[best.ring_idx].pts;
        std::vector<Point> updated = applyCollapseToRing(pts, best.start_idx, best.e);
        pts = std::move(updated);
        totalDisplacement += best.displacement;
    }

    return poly;
}

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        std::cerr << "Usage: ./simplify <input_file.csv> <target_vertices>\n";
        return 1;
    }

    try
    {
        const std::string filename = argv[1];
        const int targetVertices = std::stoi(argv[2]);

        if (targetVertices < 3)
        {
            std::cerr << "target_vertices must be at least 3\n";
            return 1;
        }

        Polygon input = readCSV(filename);
        if (!topologyValid(input))
        {
            std::cerr << "Input polygon is not topologically valid.\n";
            return 1;
        }

        const double inputArea = totalSignedArea(input);

        double totalDisplacement = 0.0;
        Polygon output = simplifyPolygon(input, targetVertices, totalDisplacement);

        const double outputArea = totalSignedArea(output);

        // Write simplified polygon to after.csv
        std::ofstream fout("after.csv");
        if (!fout)
        {
            std::cerr << "Failed to create after.csv\n";
            return 1;
        }

        writeCSVToStream(fout, output);
        fout.close();

        // Print summary to terminal
        std::cout << std::scientific << std::setprecision(15);
        std::cout << "Saved simplified polygon to after.csv\n";
        std::cout << "Total signed area in input: " << inputArea << "\n";
        std::cout << "Total signed area in output: " << outputArea << "\n";
        std::cout << "Total areal displacement: " << totalDisplacement << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}