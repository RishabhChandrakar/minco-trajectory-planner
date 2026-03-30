#ifndef CORRIDOR_BUILDER_HPP
#define CORRIDOR_BUILDER_HPP

#include <vector>
#include <cstdint>
#include <Eigen/Dense>
#include <algorithm>
#include <iostream>

#include "polyhedrongenerator_class.hpp"
#include "final_hstar.hpp"
#include "final_astar.hpp"

#include "final_obstacle_inflate.hpp"

#include <random>
#include <fstream>   
#include <string>



class CorridorBuilder
{
public:
    using Vec2i = Eigen::Vector2i;
    using Vec2d = Eigen::Vector2d;

    CorridorBuilder(int width,
                    int height,
                    double resolution,
                    double cellSizeCm,
                    double droneRadiusCm)
        : width_(width),
          height_(height),
          resolution_(resolution),
          cellSizeCm_(cellSizeCm),
          droneRadiusCm_(droneRadiusCm),
          poly_obj_(width, height)
    {
        poly_obj_._resolution = resolution_;
        poly_obj_.paramSet(false, false, true,
                           width_, height_, resolution_,
                           100, 50);

        planner_.init(width_, height_, resolution_, 16, 1.0);

        // 4. Initialize the planner memory (allocate nodes)
        astar.init(width_, height_);


    }

    // ----------------------------
    // Main entry point
    // ----------------------------
    bool build(const std::vector<uint8_t>& map,
               int sx, int sy,
               int gx, int gy)
    {
        if ((int)map.size() != width_ * height_) {
            std::cerr << "[CorridorBuilder] Map size mismatch\n";
            return false;
        }

        // Inflate obstacles
        map_inflated_ = map;

        // Load into polytope generator
        poly_obj_.loadMapData(map_inflated_);
        
        // Load into planner
        //planner_.loadMap(map_inflated_, width_, height_);

        if (planner_.isOccupied(sx, sy) || planner_.isOccupied(gx, gy)) {
           std::cerr << "[CorridorBuilder] Start or goal in obstacle\n";
           return false;
        }

        // Plan path
        //double sh = 0.0, gh = 0.0;
        //std::vector<Vec2> path =
        //    planner_.search(sx, sy, sh, gx, gy, gh, 200000);

        //if (path.empty()) {
        //    std::cerr << "[CorridorBuilder] No path found\n";
        //    return false;
        //}

        // Convert to dense grid path
        //std::vector<Vec2i> path_grid = densifyPath(path);
        


        // 5. Load the collision map
        astar.loadMap(map_inflated_, width_, height_);
        std::vector<Vec2i> path_grid = astar.searchGridPath(sx, sy, gx, gy);


        if (path_grid.empty()) {
            std::cerr << "[CorridorBuilder] No path found\n";
            return false;
        }

        logPathToFile(path_grid, "grid_path_log.csv"); 

        

        // Generate polytopes along path
        generateCorridor(path_grid);

        return true;
    }


    void logPathToFile(const std::vector<Vec2i>& path, const std::string& filename) {

        // 1. Create an output file stream
        std::ofstream outFile(filename);

        // 2. Check if the file opened successfully
        if (!outFile.is_open()) {
            std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
            return;
        }

        // 3. Header (Optional: helps if you import into Excel/Python later)
        outFile << "x,y\n";

        // 4. Loop through the vector and write coordinates
        for (const auto& point : path) {
            // Format: x, y (using '\n' is faster than std::endl)
            outFile << point.x() << "," << point.y() << "\n";
        }

        // 5. Close the file
        outFile.close();
        std::cout << "Path successfully logged to " << filename << std::endl;
    }

    // ----------------------------
    // Accessors
    // ----------------------------
    // corridor_builder.hpp -- replace your previous getInequalities implementation
    std::vector<std::pair<Eigen::MatrixXd, Eigen::VectorXd>>
    getInequalities() const
    {
        std::vector<std::pair<Eigen::MatrixXd, Eigen::VectorXd>> out;

        // poly_obj_.getCorridor() returns std::vector<Polytope>
        for (const auto& poly : poly_obj_.getCorridor()) {

            const auto& planes = poly.planes;   // std::vector<Eigen::Vector3d>

            Eigen::MatrixXd A( (int)planes.size(), 2 );
            Eigen::VectorXd b( (int)planes.size() );

            for (size_t i = 0; i < planes.size(); ++i) {

                
                const Eigen::Vector3d& p = planes[i];
                    
                A(i, 0) = p(0);   // a
                A(i, 1) = p(1);   // b
                b(i)    = -p(2);  // ax + by + c <= 0  →  ax + by <= -c
            }

            out.emplace_back(std::move(A), std::move(b));
        }

        return out;
    }




private:
    // ----------------------------
    // Path densification
    // ----------------------------
    // corridor_builder.hpp -- densifyPath replacement
    std::vector<Vec2i> densifyPath(const std::vector<Vec2>& path)
    {
        std::vector<Vec2i> dense; // if Vec2i is not appropriate, we will use auto below

        for (size_t i = 1; i < path.size(); ++i) {
            // use auto to match whatever worldToGrid() returns
            auto a = planner_.worldToGrid(path[i - 1]);
            auto b = planner_.worldToGrid(path[i]);

            // Access as fields (your planner's Vec2i likely has .x and .y fields)
            int ax = a.x; int ay = a.y;
            int bx = b.x; int by = b.y;

            int dx = std::abs(bx - ax);
            int dy = std::abs(by - ay);
            int sx = (ax < bx) ? 1 : -1;
            int sy = (ay < by) ? 1 : -1;
            int err = dx - dy;

            int x = ax;
            int y = ay;

            while (true) {
                dense.emplace_back(Eigen::Vector2i(x, y)); // store as Eigen::Vector2i
                if (x == bx && y == by) break;
                int e2 = 2 * err;
                if (e2 > -dy) { err -= dy; x += sx; }
                if (e2 <  dx) { err += dx; y += sy; }
            }
        }

        // Remove duplicates (lambda updated for Vector2i)
        dense.erase(
            std::unique(dense.begin(), dense.end(),
                [](const Eigen::Vector2i& p1, const Eigen::Vector2i& p2) {
                    return p1.x() == p2.x() && p1.y() == p2.y();
                }),
            dense.end());

        return dense;
    }
    // ----------------------------
    // Gap-Free Shortcut Logic
    // ----------------------------
    bool overlapPolytope(const polyhedronGenerator::Polytope& p1, const polyhedronGenerator::Polytope& p2) {
        // ROBUST CHECK: We check if any physical grid cell of p1 is strictly inside p2's walls.
        // Checking discrete grid cells (rather than continuous vertices) guarantees a physical, navigable tunnel.
        for (size_t i = 0; i < p1.cluster_x_idx.size(); ++i) {
            double px = static_cast<double>(p1.cluster_x_idx[i]);
            double py = static_cast<double>(p1.cluster_y_idx[i]);
            bool inside = true;
            for (const auto& plane : p2.planes) {
                // ax + by + c <= 0 (epsilon of 0.01 handles floating point rounding)
                if (plane.x() * px + plane.y() * py + plane.z() > 0.01) { 
                    inside = false; 
                    break; 
                }
            }
            if (inside) return true; // Found guaranteed safe physical overlap!
        }
        return false;
    }

    void shortCutCorridor() {
        auto& hpolys = poly_obj_.getCorridorMutable();
        const int M = static_cast<int>(hpolys.size());
        if (M <= 2) return; // Need at least 3 to bypass a middleman

        std::vector<polyhedronGenerator::Polytope> htemp = hpolys;
        hpolys.clear();

        std::vector<int> keep_indices;
        keep_indices.push_back(0); // Always keep the starting polytope

        int current = 0;
        while (current < M - 1) {
            int next_to_keep = current + 1; // Fallback: keep adjacent
            
            // Forward Greedy Search: Look as far ahead as possible
            for (int j = M - 1; j > current + 1; --j) {
                // Check overlap both ways for maximum safety
                if (overlapPolytope(htemp[current], htemp[j]) || overlapPolytope(htemp[j], htemp[current])) {
                    next_to_keep = j;
                    break;
                }
            }
            keep_indices.push_back(next_to_keep);
            current = next_to_keep;
        }

        // Rebuild optimized corridor
        for (int idx : keep_indices) {
            hpolys.push_back(htemp[idx]);
        }
        std::cerr << "[SFC Optimizer] Reduced from " << M << " to " << hpolys.size() << " polytopes.\n";
    }


    // ----------------------------
    // Corridor generation
    // ----------------------------
    void generateCorridor(const std::vector<Vec2i>& path)
    {
        poly_obj_.clearCorridor();
        if (path.empty()) return;

        // --- TUNING PARAMETERS ---
        const int progress = 10;  // How far ahead to grab target 'b'
        const int range = 9;
        Eigen::Vector2i last(-1, -1);
        
        for (size_t i = 0; i < path.size(); ++i) {
            const Eigen::Vector2i& grid = path[i];
            if (grid.x() == last.x() && grid.y() == last.y()) continue;

            if (poly_obj_.getCorridor().size() >= 2 && poly_obj_.isInsideLastSecondPolytope(grid)) {
                poly_obj_.getCorridorMutable().pop_back();
            }

            if (poly_obj_.getCorridor().empty() || poly_obj_.isOutsideLatestPolytope(grid)) {
                
                // 1. Define 'a' and 'b'
                const Eigen::Vector2i& a = grid;
                size_t b_idx = std::min(i + progress, path.size() - 1);
                const Eigen::Vector2i& b = path[b_idx];

                // 2. Define the local 'bd' bounding box (safely clamped to global map dimensions)
                int minX = std::max(0, std::min(a.x(), b.x()) - range);
                int maxX = std::min(width_ - 1, std::max(a.x(), b.x()) + range);
                int minY = std::max(0, std::min(a.y(), b.y()) - range);
                int maxY = std::min(height_ - 1, std::max(a.y(), b.y()) + range);

                // 3. Apply bounds to generator
                poly_obj_.setLocalBoundingBox(minX, maxX, minY, maxY);

                // 4. Generate shape
                std::vector<int> sx = { a.x() };
                std::vector<int> sy = { a.y() };
                poly_obj_.generateAndStorePolytope(sx, sy);
            }

            last = grid;
        }
        // 5. Delete redundant polytopes to speed up GCOPTER
        //shortCutCorridor();

        // 6. Reset bounds to full map
        poly_obj_.setLocalBoundingBox(0, width_ - 1, 0, height_ - 1);
    }

private:
    int width_, height_;
    double resolution_;
    double cellSizeCm_;
    double droneRadiusCm_;

    std::vector<uint8_t> map_inflated_;

    polyhedronGenerator poly_obj_;
    HStarPlanner planner_;

    // 3. Instantiate the AStarPlanner
    AStarPlanner astar;
};

#endif // CORRIDOR_BUILDER_HPP