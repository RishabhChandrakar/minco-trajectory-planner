#ifndef HSTAR_2D_FINAL_HPP
#define HSTAR_2D_FINAL_HPP

#include <vector>
#include <queue>
#include <cmath>
#include <limits>
#include <cstring>
#include <iostream>
#include <algorithm>

using namespace std;

// ... (Vec2 and Vec2i structs remain the same) ...
struct Vec2 {
    double x, y;
    Vec2() : x(0), y(0) {}
    Vec2(double _x, double _y) : x(_x), y(_y) {}
    Vec2 operator+(const Vec2& v) const { return Vec2(x + v.x, y + v.y); }
    Vec2 operator-(const Vec2& v) const { return Vec2(x - v.x, y - v.y); }
    Vec2 operator*(double s) const { return Vec2(x * s, y * s); }
    double norm() const { return sqrt(x*x + y*y); }
};

struct Vec2i {
    int x, y;
    Vec2i() : x(0), y(0) {}
    Vec2i(int _x, int _y) : x(_x), y(_y) {}
    bool operator==(const Vec2i& v) const { return x==v.x && y==v.y; }
    bool operator!=(const Vec2i& v) const { return !(*this==v); }
};

struct HNode {
    int gx, gy;
    int hbin;
    double cx, cy;
    double heading;
    double g, f;
    int parent;
    int8_t state;

    HNode() : gx(0), gy(0), hbin(0), cx(0), cy(0), heading(0),
          g(numeric_limits<double>::infinity()), f(numeric_limits<double>::infinity()),
          parent(-1), state(0) {}
};

struct MotionPrimitive {
    double body_dx, body_dy;
    double dheading;
    double cost;
    MotionPrimitive(double dx, double dy, double dh, double c)
        : body_dx(dx), body_dy(dy), dheading(dh), cost(c) {}
};

class HStarPlanner {
private:
    int max_x, max_y;
    double resolution;
    int heading_bins;
    double heading_res;
    uint8_t* occupancy;
    HNode* nodes;
    size_t node_count;
    vector<MotionPrimitive> motions;
    double min_turn_radius;
    double tie_breaker;

public:
    HStarPlanner() : max_x(0), max_y(0), resolution(1.0), occupancy(nullptr), nodes(nullptr) {}
    ~HStarPlanner() { cleanup(); }

    void init(int grid_x, int grid_y, double res = 1.0, int h_bins = 16, double turn_radius = 2.0) {
        cleanup();
        max_x = grid_x;
        max_y = grid_y;
        resolution = res;
        heading_bins = max(4, h_bins);
        heading_res = 2.0 * M_PI / heading_bins;
        min_turn_radius = max(0.5, turn_radius);

        size_t grid_size = (size_t)max_x * max_y;
        occupancy = new uint8_t[grid_size];
        memset(occupancy, 0, grid_size);

        // FIX: Ensure correct node count allocation
        node_count = grid_size * heading_bins;
        nodes = new HNode[node_count];
        initNodes();
        buildMotionPrimitives();
    }

    // ================= MAP LOADING (FIXED) =================
    void loadMap(const vector<uint8_t>& map, int w, int h) {
        // Ensure we don't overflow if map size mismatches init size
        int limits_x = min(w, max_x);
        int limits_y = min(h, max_y);
        
        for (int y = 0; y < limits_y; y++) {
            for (int x = 0; x < limits_x; x++) {
                if (map[y * w + x]) {
                    // Use the corrected cellIndex
                    occupancy[cellIndex(x, y)] = 1; 
                }
            }
        }
    }

    bool isOccupied(int gx, int gy) const {
        if (!inBounds(gx, gy)) return true;
        return occupancy[cellIndex(gx, gy)] != 0;
    }

    // ================= SEARCH (FIXED) =================
    vector<Vec2> search(int sx, int sy, double sheading,
                        int gx, int gy, double gheading = 0.0,
                        int max_iter = 200000)
    {
        if (!inBounds(sx, sy) || !inBounds(gx, gy)) return {};
        if (isOccupied(sx, sy) || isOccupied(gx, gy)) return {};

        sheading = normalizeAngle(sheading);
        int sbin = angleToBin(sheading);
        resetNodes();

        size_t start_idx = nodeIndex(sx, sy, sbin);
        HNode& start = nodes[start_idx];
        start.g = 0.0;
        start.f = heuristic(sx, sy, gx, gy);
        start.state = 1;

        struct Q { size_t idx; double f; };
        auto cmp = [](const Q& a, const Q& b) { return a.f > b.f; };
        priority_queue<Q, vector<Q>, decltype(cmp)> pq(cmp);
        pq.push({start_idx, start.f});

        int iterations = 0;

        while (!pq.empty() && iterations++ < max_iter) {
            auto q = pq.top(); pq.pop();
            HNode& cur = nodes[q.idx];

            if (cur.state == -1) continue;
            if (cur.gx == gx && cur.gy == gy) return reconstructPath(q.idx);

            cur.state = -1;

            double cos_h = cos(cur.heading);
            double sin_h = sin(cur.heading);

            for (const auto& m : motions) {
                double world_dx = m.body_dx * cos_h - m.body_dy * sin_h;
                double world_dy = m.body_dx * sin_h + m.body_dy * cos_h;

                double next_cx = cur.cx + world_dx;
                double next_cy = cur.cy + world_dy;
                double next_heading = normalizeAngle(cur.heading + m.dheading);

                int ngx = static_cast<int>(round(next_cx / resolution));
                int ngy = static_cast<int>(round(next_cy / resolution));

                if (!inBounds(ngx, ngy)) continue;

                // --- FIX 2: LINE COLLISION CHECK ---
                // We must check if the path FROM (cur.gx, cur.gy) TO (ngx, ngy) is clear.
                // Just checking endpoints causes "tunneling" through walls.
                if (!isLineFree(cur.gx, cur.gy, ngx, ngy)) continue;

                int hbin = angleToBin(next_heading);
                size_t ni = nodeIndex(ngx, ngy, hbin);
                HNode& nb = nodes[ni];

                if (nb.state == -1) continue;

                double ng = cur.g + m.cost;
                if (nb.state != 1 || ng < nb.g) {
                    nb.parent = static_cast<int>(q.idx);
                    nb.g = ng;
                    nb.f = ng + tie_breaker * heuristic(ngx, ngy, gx, gy);
                    nb.state = 1;
                    pq.push({ni, nb.f});
                }
            }
        }
        return {};
    }

    vector<Vec2i> searchGridPath(int sx, int sy, double sheading, int gx, int gy) {
        vector<Vec2> world = search(sx, sy, sheading, gx, gy);
        if (world.empty()) return {};

        vector<Vec2i> sparse;
        for (const auto& p : world) sparse.push_back(worldToGrid(p));
        return densifyGridPath(sparse);
    }

    Vec2i worldToGrid(const Vec2& w) const {
        return Vec2i((int)round(w.x / resolution), (int)round(w.y / resolution));
    }

private:
    // --- FIX 1: CORRECT INDEXING ---
    inline size_t cellIndex(int x, int y) const {
        // Must match Row-Major logic: y * width + x
        return (size_t)y * max_x + x;
    }

    inline size_t nodeIndex(int x, int y, int h) const {
        // Must match cellIndex logic for consistency
        return ((size_t)y * max_x + x) * heading_bins + h;
    }

    // --- NEW HELPER: BRESENHAM LINE CHECK ---
    bool isLineFree(int x0, int y0, int x1, int y1) const {
        int dx = abs(x1 - x0), dy = abs(y1 - y0);
        int sx = x0 < x1 ? 1 : -1;
        int sy = y0 < y1 ? 1 : -1;
        int err = dx - dy;

        int x = x0, y = y0;
        
        while (true) {
            if (isOccupied(x, y)) return false;
            if (x == x1 && y == y1) break;
            
            int e2 = 2 * err;
            if (e2 > -dy) { err -= dy; x += sx; }
            if (e2 < dx) { err += dx; y += sy; }
        }
        return true;
    }

    // UNCHANGED HELPERS
    vector<Vec2i> densifyGridPath(const vector<Vec2i>& s) const {
        vector<Vec2i> d;
        if (s.empty()) return d;
        d.push_back(s[0]);
        for (size_t i = 1; i < s.size(); i++) {
            int x0 = s[i - 1].x, y0 = s[i - 1].y;
            int x1 = s[i].x, y1 = s[i].y;
            int dx = abs(x1 - x0), dy = abs(y1 - y0);
            int sx = x0 < x1 ? 1 : -1, sy = y0 < y1 ? 1 : -1;
            int err = dx - dy;
            while (x0 != x1 || y0 != y1) {
                int e2 = 2 * err;
                if (e2 > -dy) { err -= dy; x0 += sx; }
                if (e2 < dx) { err += dx; y0 += sy; }
                d.emplace_back(x0, y0);
            }
        }
        return d;
    }

    void initNodes() {
        size_t i = 0;
        // Keep loops consistent with memory layout (Y then X is often better for caching if row-major)
        for (int y = 0; y < max_y; y++) {
            for (int x = 0; x < max_x; x++) {
                for (int h = 0; h < heading_bins; h++) {
                    // Recalculate index to be safe, or just increment if order matches nodeIndex
                    size_t idx = nodeIndex(x, y, h);
                    nodes[idx].gx = x;
                    nodes[idx].gy = y;
                    nodes[idx].cx = x * resolution;
                    nodes[idx].cy = y * resolution;
                    nodes[idx].heading = h * heading_res;
                    nodes[idx].parent = -1;
                    nodes[idx].state = 0;
                    nodes[idx].g = nodes[idx].f = numeric_limits<double>::infinity();
                }
            }
        }
    }

    void resetNodes() {
        for (size_t i = 0; i < node_count; i++) {
            nodes[i].g = nodes[i].f = numeric_limits<double>::infinity();
            nodes[i].parent = -1;
            nodes[i].state = 0;
        }
    }

    void buildMotionPrimitives() {
        motions.clear();
        double step = resolution;

        // Straight
        motions.emplace_back(step, 0.0, 0.0, step);

        // Arcs
        double rad = max(0.5, min_turn_radius);
        double ang1 = heading_res;
        double arc_len1 = rad * ang1;
        double dx1 = rad * sin(ang1);
        double dy1 = rad * (1.0 - cos(ang1));

        motions.emplace_back(dx1, dy1, ang1, arc_len1);
        motions.emplace_back(dx1, -dy1, -ang1, arc_len1);

        double ang2 = 2.0 * heading_res;
        double arc_len2 = rad * ang2;
        double dx2 = rad * sin(ang2);
        double dy2 = rad * (1.0 - cos(ang2));
        motions.emplace_back(dx2, dy2, ang2, arc_len2);
        motions.emplace_back(dx2, -dy2, -ang2, arc_len2);

        // REMOVED LARGE STEPS (2.0*step, 3.0*step)
        // These are dangerous for densification. If you really need them for speed,
        // you rely on 'isLineFree' to catch collisions. 
        // I left them out to guarantee safety, but you can add them back 
        // now that 'isLineFree' is implemented.
        
        for (MotionPrimitive& m : motions) {
            m.cost = max(1e-6, m.cost / resolution);
        }
    }

    vector<Vec2> reconstructPath(size_t idx) const {
        vector<Vec2> p;
        while (idx != (size_t)-1) {
            p.emplace_back(nodes[idx].cx, nodes[idx].cy);
            idx = nodes[idx].parent;
        }
        reverse(p.begin(), p.end());
        return p;
    }

    inline bool inBounds(int x, int y) const {
        return x >= 0 && y >= 0 && x < max_x && y < max_y;
    }
    inline double normalizeAngle(double a) const {
        while (a < 0) a += 2 * M_PI;
        while (a >= 2 * M_PI) a -= 2 * M_PI;
        return a;
    }
    inline int angleToBin(double a) const {
        return (int)round(normalizeAngle(a) / heading_res) % heading_bins;
    }
    inline double heuristic(int ax, int ay, int bx, int by) const {
        return hypot(ax - bx, ay - by);
    }
    void cleanup() {
        if(occupancy) delete[] occupancy; occupancy = nullptr;
        if(nodes) delete[] nodes; nodes = nullptr;
    }
};

#endif