#ifndef FINAL_ASTAR_HPP
#define FINAL_ASTAR_HPP

#include <vector>
#include <queue>
#include <cmath>
#include <limits>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <Eigen/Dense> // <--- ADDED EIGEN

using namespace std;

// REMOVED custom 'struct Vec2i' to prevent redefinition errors.
// We are now natively using Eigen::Vector2i (which is Eigen::Matrix<int, 2, 1>).

// ================= A* NODE =================
struct ANode {
    int gx, gy;
    double g, f;
    int parent;   
    int8_t state; 

    ANode() : gx(0), gy(0), 
              g(numeric_limits<double>::infinity()), 
              f(numeric_limits<double>::infinity()),
              parent(-1), state(0) {}
};

// ================= A* PLANNER =================
class AStarPlanner {
private:
    int max_x, max_y;
    uint8_t* occupancy;
    ANode* nodes;
    size_t node_count;
    double tie_breaker;

    // 8-connected grid motions
    const int dx[8] = {1, -1, 0, 0, 1, 1, -1, -1};
    const int dy[8] = {0, 0, 1, -1, 1, -1, 1, -1};
    const double dcost[8] = {1.0, 1.0, 1.0, 1.0, 1.41421, 1.41421, 1.41421, 1.41421};

public:
    AStarPlanner() : max_x(0), max_y(0), occupancy(nullptr), nodes(nullptr), tie_breaker(1.001) {}
    
    ~AStarPlanner() { cleanup(); }

    void init(int grid_x, int grid_y) {
        cleanup();
        max_x = grid_x;
        max_y = grid_y;

        node_count = (size_t)max_x * max_y;
        
        occupancy = new uint8_t[node_count];
        memset(occupancy, 0, node_count);

        nodes = new ANode[node_count];
        initNodes();
    }

    void loadMap(const vector<uint8_t>& map, int w, int h) {
        int limits_x = min(w, max_x);
        int limits_y = min(h, max_y);
        
        for (int y = 0; y < limits_y; y++) {
            for (int x = 0; x < limits_x; x++) {
                if (map[y * w + x]) {
                    occupancy[cellIndex(x, y)] = 1; 
                }
            }
        }
    }

    bool isOccupied(int gx, int gy) const {
        if (!inBounds(gx, gy)) return true;
        return occupancy[cellIndex(gx, gy)] != 0;
    }

    // ================= SEARCH (Now returns Eigen vectors) =================
    vector<Eigen::Vector2i> searchGridPath(int sx, int sy, int gx, int gy, int max_iter = 500000) {
        if (!inBounds(sx, sy) || !inBounds(gx, gy)) {
            cerr << "[AStar] Start or Goal out of bounds!" << endl;
            return {};
        }
        if (isOccupied(sx, sy)) {
            auto new_start = findNearestFree(sx, sy);
            if (new_start.x() < 0) {
                cerr << "[AStar] No free cell near start!" << endl;
                return {};
            }
            cerr << "[AStar] Start projected to free cell: "
                << new_start.x() << "," << new_start.y() << endl;
            sx = new_start.x();
            sy = new_start.y();
        }

        if (isOccupied(gx, gy)) {
            auto new_goal = findNearestFree(gx, gy);
            if (new_goal.x() < 0) {
                cerr << "[AStar] No free cell near goal!" << endl;
                return {};
            }
            cerr << "[AStar] Goal projected to free cell: "
                << new_goal.x() << "," << new_goal.y() << endl;
            gx = new_goal.x();
            gy = new_goal.y();
        }

        resetNodes();

        size_t start_idx = cellIndex(sx, sy);
        ANode& start = nodes[start_idx];
        start.g = 0.0;
        start.f = heuristic(sx, sy, gx, gy);
        start.state = 1; 

        struct Q { size_t idx; double f; };
        auto cmp =[](const Q& a, const Q& b) { return a.f > b.f; };
        priority_queue<Q, vector<Q>, decltype(cmp)> pq(cmp);
        
        pq.push({start_idx, start.f});

        int iterations = 0;

        while (!pq.empty() && iterations++ < max_iter) {
            auto q = pq.top(); 
            pq.pop();
            
            ANode& cur = nodes[q.idx];

            if (cur.state == -1) continue;
            
            if (cur.gx == gx && cur.gy == gy) {
                vector<Eigen::Vector2i> sparse_path = reconstructPath(q.idx);
                return densifyGridPath(sparse_path);
            }

            cur.state = -1;

            for (int i = 0; i < 8; ++i) {
                int nx = cur.gx + dx[i];
                int ny = cur.gy + dy[i];

                if (!inBounds(nx, ny) || isOccupied(nx, ny)) continue;

                if (i >= 4) { 
                    if (isOccupied(cur.gx, ny) && isOccupied(nx, cur.gy)) continue; 
                }

                size_t n_idx = cellIndex(nx, ny);
                ANode& neighbor = nodes[n_idx];

                if (neighbor.state == -1) continue; 

                double tentative_g = cur.g + dcost[i];

                if (neighbor.state != 1 || tentative_g < neighbor.g) {
                    neighbor.parent = static_cast<int>(q.idx);
                    neighbor.g = tentative_g;
                    neighbor.f = tentative_g + tie_breaker * heuristic(nx, ny, gx, gy);
                    neighbor.state = 1; 
                    pq.push({n_idx, neighbor.f});
                }
            }
        }
        cerr << "[AStar] Path not found within max iterations." << endl;
        return {};
    }

private:
    inline size_t cellIndex(int x, int y) const { return (size_t)y * max_x + x; }
    inline bool inBounds(int x, int y) const { return x >= 0 && y >= 0 && x < max_x && y < max_y; }
    inline double heuristic(int ax, int ay, int bx, int by) const { return std::hypot(ax - bx, ay - by); }

    void initNodes() {
        for (int y = 0; y < max_y; ++y) {
            for (int x = 0; x < max_x; ++x) {
                size_t idx = cellIndex(x, y);
                nodes[idx].gx = x;
                nodes[idx].gy = y;
                nodes[idx].parent = -1;
                nodes[idx].state = 0;
                nodes[idx].g = nodes[idx].f = numeric_limits<double>::infinity();
            }
        }
    }

    void resetNodes() {
        for (size_t i = 0; i < node_count; ++i) {
            nodes[i].parent = -1;
            nodes[i].state = 0;
            nodes[i].g = nodes[i].f = numeric_limits<double>::infinity();
        }
    }

    // Now uses Eigen::Vector2i
    vector<Eigen::Vector2i> reconstructPath(size_t idx) const {
        vector<Eigen::Vector2i> path;
        while (idx != (size_t)-1) {
            path.emplace_back(nodes[idx].gx, nodes[idx].gy);
            idx = nodes[idx].parent;
        }
        reverse(path.begin(), path.end());
        return path;
    }

    // Now uses Eigen::Vector2i and .x() / .y() accessors
    vector<Eigen::Vector2i> densifyGridPath(const vector<Eigen::Vector2i>& s) const {
        vector<Eigen::Vector2i> d;
        if (s.empty()) return d;
        d.push_back(s[0]);
        for (size_t i = 1; i < s.size(); i++) {
            int x0 = s[i - 1].x(), y0 = s[i - 1].y();
            int x1 = s[i].x(), y1 = s[i].y();
            int dx = abs(x1 - x0), dy = abs(y1 - y0);
            int sx = x0 < x1 ? 1 : -1, sy = y0 < y1 ? 1 : -1;
            int err = dx - dy;
            while (x0 != x1 || y0 != y1) {
                int e2 = 2 * err;
                if (e2 > -dy) { err -= dy; x0 += sx; }
                if (e2 < dx) { err += dx; y0 += sy; }
                if (x0 != x1 || y0 != y1) d.emplace_back(x0, y0);
            }
            d.emplace_back(x1, y1);
        }
        
        // Safely remove duplicates using a custom lambda for Eigen vectors
        d.erase(unique(d.begin(), d.end(),[](const Eigen::Vector2i& a, const Eigen::Vector2i& b) {
            return a.x() == b.x() && a.y() == b.y();
        }), d.end());
        
        return d;
    }

    Eigen::Vector2i findNearestFree(int sx, int sy) const {

        if (!inBounds(sx, sy))
            return Eigen::Vector2i(-1, -1);

        if (!isOccupied(sx, sy))
            return Eigen::Vector2i(sx, sy);

        std::vector<uint8_t> visited(node_count, 0);
        std::queue<Eigen::Vector2i> q;

        q.emplace(sx, sy);
        visited[cellIndex(sx, sy)] = 1;

        const int dx4[4] = {1, -1, 0, 0};
        const int dy4[4] = {0, 0, 1, -1};

        while (!q.empty()) {
            auto p = q.front();
            q.pop();

            for (int i = 0; i < 4; ++i) {
                int nx = p.x() + dx4[i];
                int ny = p.y() + dy4[i];

                if (!inBounds(nx, ny))
                    continue;

                size_t idx = cellIndex(nx, ny);
                if (visited[idx])
                    continue;

                visited[idx] = 1;

                if (!isOccupied(nx, ny)) {
                    return Eigen::Vector2i(nx, ny);
                }

                q.emplace(nx, ny);
            }
        }

        return Eigen::Vector2i(-1, -1);  // no free cell found
    }
    
    void cleanup() {
            if (occupancy) {
                delete[] occupancy;
                occupancy = nullptr;
            }
            if (nodes) {
                delete[] nodes;
                nodes = nullptr;
            }
    }

};

#endif
