#ifndef POLYHEDRON_GENERATOR_CLASS_HPP
#define POLYHEDRON_GENERATOR_CLASS_HPP

#include <iostream>
#include <Eigen/Dense>
#include <limits>
#include <string>
#include <fstream>
#include <math.h>
#include <cmath>
#include <vector>
#include <cstring>
#include <iomanip>

#include "quickhull_class.hpp"
#include "linear_eq.hpp"

using namespace std;
using namespace Eigen;

class polyhedronGenerator {

public:

    
    double _resolution;

    double _map_lower_x=0.0;
    double _map_lower_y=0.0;


    struct PolyInequality
    {
        Eigen::MatrixXd A;   // N × 2 (or N × 3 if extended to 3D)
        Eigen::VectorXd b;   // N × 1
    };

    struct Polytope
    {
        // Plane representation: ax + by + c <= 0
        std::vector<Eigen::Vector3d> planes;

        // Geometric attributes
        Eigen::Vector2d center;
        Eigen::Vector2d seed;

        std::vector<int> cluster_x_idx;     // grid X indices of convex region
        std::vector<int> cluster_y_idx;     // grid Y indices of convex region


        std::vector<Eigen::Vector2d> convex_vertices;

        PolyInequality inequality;


        void appendPlane(const Eigen::Vector3d& p) { planes.push_back(p); }

        void setCenter(const Eigen::Vector2d& c) { center = c; }

        void setSeed(const Eigen::Vector2d& s) { seed = s; }

        void addClusterIndices(const std::vector<int>& xs,const std::vector<int>& ys)
        {
            cluster_x_idx = xs;   // copy entire list
            cluster_y_idx = ys;
        }

        // Store continuous polygon vertices (optional)
        void addConvexVertex(const Eigen::Vector2d& v)
        {
            convex_vertices.push_back(v);
        }

        // NEW: Assign inequality form A·x <= b

        void setInequality(const Eigen::MatrixXd& A_in,const Eigen::VectorXd& b_in)
        {
            inequality.A = A_in;
            inequality.b = b_in;
        }
    };

    vector<Polytope> polytope_corridor;       // to store the data of each ploygon



    
    polyhedronGenerator(int max_x_id, int max_y_id) {
        _max_x_id = max_x_id;
        _max_y_id = max_y_id;
        _grid_num = max_x_id * max_y_id;
        
        // Initialize Eigen matrices
        map_data = MatrixXi::Zero(_max_x_id, _max_y_id);
        use_data = MatrixXi::Zero(_max_x_id, _max_y_id);
        invalid_data = MatrixXi::Zero(_max_x_id, _max_y_id);
        inside_data = MatrixXi::Zero(_max_x_id, _max_y_id);
        
        // Initialize vertex indices
        vertex_idx = VectorXi::Zero(8);
        vertex_idx_lst = VectorXi::Zero(8);
        
        inf_step = 1;
        itr_inflate_max = 1000;
        _cluster_buffer_size = 50000;
        _candidate_buffer_size = 10000;
    }
    
    ~polyhedronGenerator() {
        // Eigen handles memory management automatically
    }

    void mapClear() {
        map_data.setZero();
    }

    void flagClear() {
        use_data.setZero();
        invalid_data.setZero();
        inside_data.setZero();
    }
    //adding 
    // Called before each polytope generation to set the local BD box
    void setLocalBoundingBox(int minx, int maxx, int miny, int maxy) {
        limit_min_x = std::max(0,             minx);
        limit_max_x = std::min(_max_x_id - 1, maxx);
        limit_min_y = std::max(0,             miny);
        limit_max_y = std::min(_max_y_id - 1, maxy);
    }

    void paramSet(bool is_gpu_on_stage_1, bool is_gpu_on_stage_2, bool is_cluster_on,
                  const int & max_x_id, const int & max_y_id, double resolution, 
                  double itr_inflate_max_, double itr_cluster_max_) {
        _resolution = resolution;
        _max_x_id = max_x_id;
        _max_y_id = max_y_id;
        _grid_num = max_x_id * max_y_id;

        if(is_cluster_on) {
            itr_inflate_max = itr_inflate_max_;
            itr_cluster_max = itr_cluster_max_;
        }
        else {
            itr_inflate_max = 1000;
            itr_cluster_max = 0;   
        }

        _cluster_buffer_size = 50000;
        _candidate_buffer_size = 10000;
        _cluster_buffer_size_square = _candidate_buffer_size * _candidate_buffer_size;

        // Reinitialize matrices with new sizes
        map_data = MatrixXi::Zero(max_x_id, max_y_id);
        use_data = MatrixXi::Zero(max_x_id, max_y_id);
        invalid_data = MatrixXi::Zero(max_x_id, max_y_id);
        inside_data = MatrixXi::Zero(max_x_id, max_y_id);

        // Initialize buffer matrices
        active_xy_id = MatrixXi::Zero(_candidate_buffer_size, 2);
        cluster_xy_id = MatrixXi::Zero(_cluster_buffer_size, 2);
        candidate_cell = MatrixXi::Zero(_candidate_buffer_size, 2);
        //Adding 
        limit_min_x = 0;
        limit_max_x = max_x_id - 1;
        limit_min_y = 0;
        limit_max_y = max_y_id - 1;
        
        inf_step = 1;
        mapClear();
    }
    void setObs(const int & id_x, const int & id_y) {
        if(id_x >= 0 && id_x < _max_x_id && id_y >= 0 && id_y < _max_y_id) {
            map_data(id_x, id_y) = 1;
        }
    }

    void setFr(const int & id_x, const int & id_y) {
        if(id_x >= 0 && id_x < _max_x_id && id_y >= 0 && id_y < _max_y_id) {
            map_data(id_x, id_y) = 0;
        }
    }

    void setVertexInitIndex(VectorXi& vertex_idx, int min_x, int min_y, int max_x, int max_y) {
        vertex_idx(0) = vertex_idx(1) = max_x;  
        vertex_idx(2) = vertex_idx(3) = min_x;  
        vertex_idx(4) = vertex_idx(7) = min_y;
        vertex_idx(5) = vertex_idx(6) = max_y;
    }

    void getGridsInRect(const VectorXi& vertex_idx, 
                         vector<int>& cube_grid_x, 
                         vector<int>& cube_grid_y) {
        int id_x, id_y;
        int min_x = vertex_idx(2);
        int max_x = vertex_idx(0);
        int min_y = vertex_idx(4);
        int max_y = vertex_idx(5);
        
        for(id_x = min_x; id_x <= max_x; id_x++) {   
            for(id_y = min_y; id_y <= max_y; id_y++) {
                cube_grid_x.push_back(id_x);
                cube_grid_y.push_back(id_y);
                if(map_data(id_x,id_y)!=1) 
                    {inside_data(id_x, id_y) = 2;}
            }
        }
    }

    void inflateX_n(VectorXi& vertex_idx) {
        if(vertex_idx(0 + 3) <= limit_min_x) return;     

        int id_x = vertex_idx(0 + 3) - 1;   
        for(int id_y = vertex_idx(4 + 3); id_y <= vertex_idx(4 + 2); id_y++) {
            if(map_data(id_x, id_y) == 1) {   
                return;
            }
        }

        vertex_idx(2) -= inf_step;
        vertex_idx(3) -= inf_step;
    }

    void inflateX_p(VectorXi& vertex_idx) {
        if(vertex_idx(0 + 0) >= limit_max_x) return;     

        int id_x = vertex_idx(0 + 0) + 1;   
        for(int id_y = vertex_idx(4 + 0); id_y <= vertex_idx(4 + 1); id_y++) {
            if(map_data(id_x, id_y) ==1) {   
                return;
            }
        }

        vertex_idx(0) += inf_step;
        vertex_idx(1) += inf_step;
    }

    void inflateY_n(VectorXi& vertex_idx) {
        if(vertex_idx(4 + 0) <= limit_min_y) return;     

        int id_y = vertex_idx(4 + 0) - 1;   
        for(int id_x = vertex_idx(0 + 3); id_x <= vertex_idx(0 + 0); id_x++) {
            if(map_data(id_x, id_y)==1) {   
                return;
            }
        }

        vertex_idx(4) -= inf_step;
        vertex_idx(7) -= inf_step;
    }

    void inflateY_p(VectorXi& vertex_idx) {
        if(vertex_idx(4 + 1) >= limit_max_y) return;     

        int id_y = vertex_idx(4 + 1) + 1;   
        for(int id_x = vertex_idx(0 + 2); id_x <= vertex_idx(0 + 1); id_x++) {
            if(map_data(id_x, id_y) ==1) {   
                return;
            }
        }

        vertex_idx(5) += inf_step;
        vertex_idx(6) += inf_step;
    }

    void cubeInflation_cpu(VectorXi& vertex_idx_lst, VectorXi& vertex_idx) {   
        int iter = 0;
        while(iter < itr_inflate_max) {   
            vertex_idx_lst = vertex_idx;

            inflateY_n(vertex_idx);
            inflateY_p(vertex_idx);
            inflateX_n(vertex_idx);
            inflateX_p(vertex_idx);

            bool is_inflate_conti = false;
            for(int vtx = 0; vtx < 4; vtx++) {
                if((vertex_idx_lst(vtx) != vertex_idx(vtx)) || 
                   (vertex_idx_lst(vtx + 4) != vertex_idx(vtx + 4))) {
                    is_inflate_conti = true;
                    break;
                }
            }

            if(is_inflate_conti == false)
                break;

            for(int vtx = 0; vtx < 4; vtx++) { 
                vertex_idx_lst(vtx + 0) = vertex_idx(vtx + 0);
                vertex_idx_lst(vtx + 4) = vertex_idx(vtx + 4);
            }

            iter++;
        }
        
        // cout << "  Inflation completed after " << iter << " iterations" << endl;
    }

    bool serialConvexTest2D(int can_x_index, int can_y_index,
                            int cluster_grid_num, int max_y_id,
                            const MatrixXi &cluster_xy_id,
                            const MatrixXi &inside_data,
                            const MatrixXi &map_data)
    {
        Vector2i index_1_med(can_x_index / 2, can_y_index / 2);
        Vector2i index_med, index_2;
        
        for (int i = cluster_grid_num - 1; i >= 0; i--)
        {
            index_2(0) = cluster_xy_id(i, 0);
            index_2(1) = cluster_xy_id(i, 1);

            if (inside_data(can_x_index, can_y_index) == 0)
            {
                index_med = index_1_med + (index_2/2);

                if (inside_data(index_med(0), index_med(1)) == 1)
                    continue;

                int x = can_x_index, y = can_y_index;
                int endX = index_2(0), endY = index_2(1);
                int dx = endX - x, dy = endY - y;
                int stepX = signum_cpu(dx), stepY = signum_cpu(dy);

                float tMaxX = intbound_cpu(0.5, dx);
                float tMaxY = intbound_cpu(0.5, dy);
                float tDeltaX = ((float)stepX) / dx;
                float tDeltaY = ((float)stepY) / dy;
                
                // Ray traversal
                while (true)
                {   
                    if (tMaxX < tMaxY)
                    {
                        x += stepX;
                        tMaxX += tDeltaX;
                    }
                    else if(tMaxY < tMaxX)
                    {
                        y += stepY;
                        tMaxY += tDeltaY;
                    }
                    else{
                        y += stepY;
                        x += stepX;
                        tMaxY += tDeltaY;
                        tMaxX += tDeltaX;
                    }

                    if (x == endX && y == endY)
                        break;
                    
                    if (inside_data(x, y) == 1)
                    {
                        continue;
                    }
                    else if (map_data(x, y) == 1)
                    {
                        return false;
                    }
                }
            }
        }

        return true;
    }
    
    void polytopeCluster_cpu(int& cluster_grid_num, int& active_grid_num)
    {   
        int candidate_grid_num;

        int itr_cluster_cnt = 0;
        while(itr_cluster_cnt < itr_cluster_max)
        {   
            candidate_grid_num = 0;
            for(int i = 0; i < active_grid_num; i++)
            {
                Eigen::Vector2i cur_cell(active_xy_id(i, 0), active_xy_id(i, 1));
                
                use_data(cur_cell.x(), cur_cell.y()) = 1;

                // Get all nearby voxels (8-connectivity in 2D)
                for(int dx = -1; dx < 2; dx++)
                { 
                    for(int dy = -1; dy < 2; dy++)
                    {   
                        if(dx == 0 && dy == 0) continue;
                        
                        Eigen::Vector2i nei_cell = cur_cell + Eigen::Vector2i(dx, dy);

                        // Bounds check
                        if(nei_cell.x() < limit_min_x || nei_cell.x() > limit_max_x
                        || nei_cell.y() < limit_min_y || nei_cell.y() > limit_max_y)
                            continue;

                        // Check if cell is valid candidate
                        if(map_data(nei_cell.x(), nei_cell.y()) == 1 
                        || use_data(nei_cell.x(), nei_cell.y()) == 1 
                        || invalid_data(nei_cell.x(), nei_cell.y()) == 1 
                        || inside_data(nei_cell.x(), nei_cell.y()) == 1)
                        {
                            continue;
                        }
                        else
                        {   
                            candidate_cell(candidate_grid_num, 0) = nei_cell.x();
                            candidate_cell(candidate_grid_num, 1) = nei_cell.y();

                            candidate_grid_num++;
                            use_data(nei_cell.x(), nei_cell.y()) = 1;
                        }
                    }
                }
            }

            if(candidate_grid_num == 0) 
                break;

            active_grid_num = 0;
            for(int i = 0; i < candidate_grid_num; i++)
            {
                // Test if candidate preserves convex hull property
                int cand_x = candidate_cell(i, 0);
                int cand_y = candidate_cell(i, 1);

                if(serialConvexTest2D(cand_x, cand_y, 
                                cluster_grid_num, _max_y_id, 
                                cluster_xy_id, inside_data, map_data))
                {   
                    cluster_xy_id(cluster_grid_num, 0) = cand_x;
                    cluster_xy_id(cluster_grid_num, 1) = cand_y;

                    active_xy_id(active_grid_num, 0) = cand_x;
                    active_xy_id(active_grid_num, 1) = cand_y;

                    cluster_grid_num++;
                    active_grid_num++;

                    inside_data(cand_x, cand_y) = 0;
                }
                else
                {   
                    invalid_data(cand_x, cand_y) = 1;
                }
            }
            
            if(active_grid_num == 0) 
                break;

            itr_cluster_cnt++; 
        }
        
        // cout << "  Clustering completed after " << itr_cluster_cnt << " iterations" << endl;
        // cout << "  Total cluster grid count: " << cluster_grid_num << endl;
    }

    void polygonGeneration(
        std::vector<int>& cluster_x_idx, std::vector<int>& cluster_y_idx)
    {   
        flagClear();
        
        if(cluster_x_idx.size() == 1)
        {
            // Single point
            Eigen::Vector2i single_point(cluster_x_idx[0], cluster_y_idx[0]);
            for(int vtx = 0; vtx < 4; vtx++)
            {   
                vertex_idx(vtx) = single_point.x();
                vertex_idx(vtx + 4) = single_point.y();
            }       
        }
        else
        {
            // Find bounding box
            Eigen::Vector2i min_corner(100000, 100000);
            Eigen::Vector2i max_corner(-100000, -100000);
            
            for(int i = 0; i < (int)cluster_x_idx.size(); i++)
            {
                Eigen::Vector2i point(cluster_x_idx[i], cluster_y_idx[i]);
                min_corner = min_corner.cwiseMin(point);
                max_corner = max_corner.cwiseMax(point);
            }

            // Initialize rectangle vertices
            setVertexInitIndex(vertex_idx, min_corner.x(), min_corner.y(), 
                            max_corner.x(), max_corner.y());
        }

        // Copy vertex indices
        for(int i = 0; i < 8; i++) {
            vertex_idx_lst(i) = vertex_idx(i);
        }

        // Rectangle inflation
        // cout << "  Starting cube inflation..." << endl;
        cubeInflation_cpu(vertex_idx_lst, vertex_idx);

        cluster_x_idx.clear(); 
        cluster_y_idx.clear();

        // Get grids within rectangle
        std::vector<int> rect_grid_x, rect_grid_y;
        getGridsInRect(vertex_idx, rect_grid_x, rect_grid_y);
        
        std::vector<Eigen::Vector2i> rect_outside_grids;

        if(rect_grid_x.size() == 1)
        {
            rect_outside_grids.push_back(Eigen::Vector2i(rect_grid_x[0], rect_grid_y[0]));
        }
        else
        {
            for(int i = 0; i < (int)rect_grid_x.size(); i++)
            {   
                Eigen::Vector2i grid_cell(rect_grid_x[i], rect_grid_y[i]);
                
                use_data(grid_cell.x(), grid_cell.y()) = 1;
                bool is_boundary=false;
                //int is_inside = 1;      
                
                // Check 8-neighborhood
                for(int dx = -1; dx < 2; dx++)
                {   
                    for(int dy = -1; dy < 2; dy++)
                    {   
                        if(dx == 0 && dy == 0)
                            continue;

                        Eigen::Vector2i neighbor = grid_cell + Eigen::Vector2i(dx, dy);

                        // Bounds check
                        if(neighbor.x() >= 0 && neighbor.x() < _max_x_id && 
                        neighbor.y() >= 0 && neighbor.y() < _max_y_id)
                        {
                            // is_inside *= inside_data(neighbor.x(), neighbor.y());
                            if(inside_data(neighbor.x(), neighbor.y())!=2){
                                is_boundary=true;
                                break;
                            }
                        }
                        else{
                            //is_inside = 0;   
                            is_boundary= true;
                            break;}
                    }
                    if(is_boundary) break;
                }

                if(is_boundary) // Boundary grid
                {
                    rect_outside_grids.push_back(grid_cell);
                }
            }
        }
        
        // Mark boundary grids
        for(const auto& grid : rect_outside_grids)
        {
            inside_data(grid.x(), grid.y()) = 0;
        }

        int init_cluster_grid_num = rect_outside_grids.size();
        int active_grid_num = init_cluster_grid_num;

        // Initialize cluster and active grids
        for(int i = 0; i < init_cluster_grid_num; i++)
        {
            cluster_xy_id(i, 0) = rect_outside_grids[i].x();
            cluster_xy_id(i, 1) = rect_outside_grids[i].y();

            active_xy_id(i, 0) = rect_outside_grids[i].x();
            active_xy_id(i, 1) = rect_outside_grids[i].y();
        }

        // Check if rectangle is degenerate
        Eigen::Vector2i vertex_diff_1(
            vertex_idx(3) - vertex_idx(1),
            vertex_idx(7) - vertex_idx(5)
        );

        if(vertex_diff_1.x() == 0 || vertex_diff_1.y() == 0)
        {   
            for(int i = 0; i < init_cluster_grid_num; i++)
            {
                cluster_x_idx.push_back(cluster_xy_id(i, 0));
                cluster_y_idx.push_back(cluster_xy_id(i, 1));
            }

            for (int i = 0; i < rect_grid_x.size(); ++i) {
                cluster_x_idx.push_back(rect_grid_x[i]);
                cluster_y_idx.push_back(rect_grid_y[i]);
            }

            return;
        }

        int cluster_grid_num = init_cluster_grid_num;

        // Expand cluster
        // cout << "  Starting polytope clustering..." << endl;

        // Extract final cluster

        
        polytopeCluster_cpu(cluster_grid_num, active_grid_num);

        // Add INSIDE rectangle points FIRST
        for (int i = 0; i < rect_grid_x.size(); ++i) {
            cluster_x_idx.push_back(rect_grid_x[i]);
            cluster_y_idx.push_back(rect_grid_y[i]);
        }

        
        for(int i = 0; i < cluster_grid_num; i++)
        {
            cluster_x_idx.push_back(cluster_xy_id(i, 0));
            cluster_y_idx.push_back(cluster_xy_id(i, 1));
        }

        // cout << "  [polyhedron_generator]{CPU} Finished polygon generation" << endl;
    }


    void printPolytopeInequalities() const
    {
        std::cout << "\n=============================================\n";
        std::cout << " POLYTOPE CORRIDOR INEQUALITY DEBUG OUTPUT\n";
        std::cout << " Total polytopes: " << polytope_corridor.size() << "\n";
        std::cout << "=============================================\n";

        // ============================================================
        // PRINT POLYTOPE CORRIDOR DETAILS
        // ============================================================
        for (size_t pid = 0; pid < polytope_corridor.size(); ++pid)
        {
            const Polytope& P = polytope_corridor[pid];

            std::cout << "\n=============================================\n";
            std::cout << " Polytope ID : " << pid << "\n";
            std::cout << "=============================================\n";

            // ----------------------------------------------------------
            // Seed & Center
            // ----------------------------------------------------------
            std::cout << " Seed   : (" << P.seed.x()   << ", " << P.seed.y()   << ")\n";
            std::cout << " Center : (" << P.center.x() << ", " << P.center.y() << ")\n";

            // ----------------------------------------------------------
            // Cluster grid indices (x,y pairs)
            // ----------------------------------------------------------
            size_t cluster_sz = std::min(P.cluster_x_idx.size(),
                                        P.cluster_y_idx.size());

            std::cout << "\n Cluster cells (" << cluster_sz << "):\n   ";
            std::cout<<P.cluster_x_idx.size()<<endl;

            for (size_t i = 0; i < cluster_sz; ++i)
            {
                std::cout << "("
                        << P.cluster_x_idx[i] << ", "
                        << P.cluster_y_idx[i] << ") ,";

                if ((i + 1) % 8 == 0) std::cout << "\n   ";
            }
            std::cout << "\n";

            // ----------------------------------------------------------
            // Convex polygon vertices (continuous coordinates)
            // ----------------------------------------------------------
            std::cout << "\n Convex vertices (" << P.convex_vertices.size() << "):\n";

            for (size_t i = 0; i < P.convex_vertices.size(); ++i)
            {
                const auto& v = P.convex_vertices[i];
                std::cout << "   [" << std::setw(2) << i << "]  ("
                        << std::fixed << std::setprecision(6)
                        << v.x() << ", " << v.y() << ")\n";
            }

            // ----------------------------------------------------------
            // Inequality representation: A * x <= b
            // ----------------------------------------------------------
            const auto& A = P.inequality.A;
            const auto& b = P.inequality.b;

            std::cout << "\n Inequality form:  A * x <= b\n";
            std::cout << "   A size : " << A.rows() << " x " << A.cols() << "\n";
            std::cout << "   b size : " << b.size() << "\n";

            std::cout << "\n Constraints:\n";
            for (int i = 0; i < b.size(); ++i)
            {
                std::cout << "   [" << std::setw(2) << i << "]  "
                        << std::fixed << std::setprecision(6)
                        << A(i, 0) << " * x  +  "
                        << A(i, 1) << " * y  <=  "
                        << b(i) << "\n";
            }

            // ----------------------------------------------------------
            // Plane representation: ax + by + c <= 0
            // ----------------------------------------------------------
            std::cout << "\n Plane representation (ax + by + c <= 0):\n";

            for (size_t i = 0; i < P.planes.size(); ++i)
            {
                const auto& pl = P.planes[i];
                std::cout << "   [" << std::setw(2) << i << "]  "
                        << std::fixed << std::setprecision(6)
                        << pl.x() << " * x  +  "
                        << pl.y() << " * y  +  "
                        << pl.z() << " <= 0\n";
            }
        }


        std::cout << "\n=============================================\n";
    }


    void generateAndStorePolytope(
    std::vector<int>& cluster_x_idx,
    std::vector<int>& cluster_y_idx)
    {
        // -------------------------------
        // 0. Safety
        // -------------------------------
        if (cluster_x_idx.empty() || cluster_y_idx.empty())
            return;

        // -------------------------------
        // 1. Seed (first grid point)
        // -------------------------------
        const int seed_gx = cluster_x_idx[0];
        const int seed_gy = cluster_y_idx[0];

        Polytope P;
        P.setSeed(Eigen::Vector2d(seed_gx, seed_gy));

        // -------------------------------
        // 2. Polygon grid clustering
        // -------------------------------
        polygonGeneration(cluster_x_idx, cluster_y_idx);

        if (cluster_x_idx.empty()) {
            std::cout << "[polyGen] Empty cluster, skipping polytope\n";
            return;
        }

        // -------------------------------
        // 3. Store cluster grid indices
        // -------------------------------
        P.addClusterIndices(cluster_x_idx, cluster_y_idx);

        // -------------------------------
        // 4. Compute center (mean of grids)
        // -------------------------------
        Eigen::Vector2d center(0.0, 0.0);
        for (size_t i = 0; i < cluster_x_idx.size(); ++i) {
            center.x() += cluster_x_idx[i];
            center.y() += cluster_y_idx[i];
        }
        center /= static_cast<double>(cluster_x_idx.size());
        P.setCenter(center);

        // -------------------------------
        // 5. Build point cloud for QuickHull
        // -------------------------------
        std::vector<sfc_quickhull::Vector2<double>> pts;
        pts.reserve(cluster_x_idx.size());

        for (size_t i = 0; i < cluster_x_idx.size(); ++i) {
            pts.emplace_back(
                static_cast<double>(cluster_x_idx[i]),
                static_cast<double>(cluster_y_idx[i])
            );
        }

        // -------------------------------
        // 6. Compute convex hull
        // -------------------------------
        sfc_quickhull::QuickHulld qh;
        auto hull = qh.getConvexHull(pts, true, false);

        if (hull.vertices().size() < 3) {
            std::cout << "[polyGen] Degenerate hull, skipping\n";
            return;
        }

        // -------------------------------
        // 7. Store convex hull vertices
        // -------------------------------
        std::vector<line_eq::LineEquations::Point> hull_pts;
        hull_pts.reserve(hull.vertices().size());

        for (size_t i = 0; i < hull.vertices().size(); ++i) {
            const auto& v = hull.vertices()[i];

            Eigen::Vector2d hv(v.x, v.y);
            P.addConvexVertex(hv);

            hull_pts.emplace_back(
                static_cast<int>(std::round(v.x)),
                static_cast<int>(std::round(v.y))
            );
        }

                cout<<"from polyhedron = "<<" "<<P.seed.x()<<" "<<P.seed.y()<<endl;

        line_eq::LineEquations::Point seed_pt(
            static_cast<int>(std::round(P.seed.x())),
            static_cast<int>(std::round(P.seed.y()))
        );


        // -------------------------------
        // 8. Compute Ax ≤ b (inequality)
        // -------------------------------
        line_eq::LineEquations le(hull_pts,seed_pt);
        le.compute();

        P.setInequality(le.A(), le.b());

        // -------------------------------
        // 9. Convert Ax ≤ b → ax + by + c ≤ 0
        // -------------------------------
        const Eigen::MatrixXd& A = le.A();
        const Eigen::VectorXd& b = le.b();

        for (int i = 0; i < A.rows(); ++i) {
            const double a = A(i, 0);
            const double bb = A(i, 1);
            const double c = -b(i);   // IMPORTANT

            P.appendPlane(Eigen::Vector3d(a, bb, c));
        }

        // Right wall (Max X):  1x + 0y - limit_max_x <= 0
        P.appendPlane(Eigen::Vector3d(1.0, 0.0, -static_cast<double>(limit_max_x)));
        
        // Left wall (Min X):  -1x + 0y + limit_min_x <= 0
        P.appendPlane(Eigen::Vector3d(-1.0, 0.0, static_cast<double>(limit_min_x)));
        
        // Top wall (Max Y):    0x + 1y - limit_max_y <= 0
        P.appendPlane(Eigen::Vector3d(0.0, 1.0, -static_cast<double>(limit_max_y)));
        
        // Bottom wall (Min Y): 0x - 1y + limit_min_y <= 0
        P.appendPlane(Eigen::Vector3d(0.0, -1.0, static_cast<double>(limit_min_y)));


        // --------------------------------------------------------------------------
        // VALIDATION: Check for Internal Obstacles (Continuation)
        // --------------------------------------------------------------------------
        
        bool obstacle_inside = false;
        const double eps = 1e-3; // Tolerance for floating point math

        for (const auto& obs : _occupied_cells) {
            double ox = static_cast<double>(obs.x());
            double oy = static_cast<double>(obs.y());

            // 1. FAST CHECK: Bounding Box
            // If the obstacle is outside the rectangular limits we just set, 
            // it cannot be inside the polytope. Skip it.
            if (ox < limit_min_x || ox > limit_max_x || 
                oy < limit_min_y || oy > limit_max_y) {
                continue;
            }

            // 2. PRECISE CHECK: Hull Equations
            // The obstacle is inside the box, now check if it is inside the angled lines (Ax <= b).
            // A point is inside if it satisfies ALL plane equations.
            bool inside_hull = true;
            for (int i = 0; i < A.rows(); ++i) {
                // Evaluate: ax + by - b <= 0
                double val = A(i, 0) * ox + A(i, 1) * oy - b(i);
                
                // If val > 0, the point is outside this specific line -> Safe
                if (val > eps) { 
                    inside_hull = false; 
                    break; 
                }
            }

            // If it passed the Box check AND the Hull check, it is inside the Polytope.
            if (inside_hull) {
                obstacle_inside = true;
                break; // Found a collision, no need to check further
            }
        }

        // --------------------------------------------------------------------------
        // DECISION: Save or Discard
        // --------------------------------------------------------------------------
        
        if (obstacle_inside) {
            // Discard logic
            std::cout << "[CorridorBuilder] Polytope discarded due to internal obstacle." << std::endl;
        } else {
            // Valid! Store the polytope in your final list
            // Assuming you have a std::vector<Polytope> called 'polytopes' or similar:

            // -------------------------------
            // 10. Store polytope
            // -------------------------------
            polytope_corridor.push_back(P);

            std::cout << "[polyGen] Stored polytope | planes="
            << P.planes.size()
            << " | hull_vertices="
            << P.convex_vertices.size()
            << " | cells="
            << P.cluster_x_idx.size()
            << "\n";
        }





    }



    // ====== NEW COORDINATE CONVERSION FUNCTIONS ======
    void coord2Index(int& id_x, int& id_y, const double& pt_x, const double& pt_y) {
        id_x = std::min(std::max(int((pt_x - _map_lower_x) * _resolution), 0), _max_x_id - 1);
        id_y = std::min(std::max(int((pt_y - _map_lower_y) * _resolution), 0), _max_y_id - 1);
    }

    Vector2i coord2Index(Vector2d coord) {
        int id_x, id_y;
        coord2Index(id_x, id_y, coord(0), coord(1));
        return Vector2i(id_x, id_y);
    }

    Vector2d index2Coord(Vector2i index) {
        return Vector2d(index(0) * _resolution + 0.5 * _resolution + _map_lower_x,
                        index(1) * _resolution + 0.5 * _resolution + _map_lower_y);
    }

    Vector2d index2Coord(int id_x, int id_y) {
        return Vector2d(id_x * _resolution + 0.5 * _resolution + _map_lower_x,
                        id_y * _resolution + 0.5 * _resolution + _map_lower_y);
    }

    // ====== NEW POLYTOPE CONTAINMENT CHECKS ======
    bool isOutsidePolytope(const Vector2i& cur_coord, const Polytope& polytope) {
        for(const auto& plane : polytope.planes) {
            double res = cur_coord.cast<double>().dot(plane.head<2>()) + plane(2);

            if(res > 0.01) return true;
        }
        return false;
    }

    bool isInsidePolytope(const Vector2i& cur_coord, const Polytope& polytope) {
        return !isOutsidePolytope(cur_coord, polytope);
    }

    bool isOutsideLatestPolytope(const Vector2i& cur_coord) {
        if(polytope_corridor.empty()) return true;
        return isOutsidePolytope(cur_coord, polytope_corridor.back());
    }

    bool isInsideLastSecondPolytope(const Vector2i& cur_coord) {
        if(polytope_corridor.size() < 2) return false;
        return isInsidePolytope(cur_coord, polytope_corridor[polytope_corridor.size()-2]);
    }


    // ====== NEW POLYTOPE CONVERSION FUNCTION ======
    Polytope convertToPolytope(const vector<int>& poly_x, const vector<int>& poly_y) {
        Polytope polytope;
        
        // QuickHull
        vector<sfc_quickhull::Vector2<double>> points;
        for(size_t i = 0; i < poly_x.size(); i++) {
            points.emplace_back(poly_x[i], poly_y[i]);
        }

        sfc_quickhull::QuickHull<double> qh;
        auto hull = qh.getConvexHull(points.data(), points.size(), true, false);
        
        // LineEquations
        vector<line_eq::LineEquations::Point> hull_points;
        for(const auto& v : hull.vertices()) {
            hull_points.emplace_back(round(v.x), round(v.y));
        }

        std::pair<int,int> seed_ij{
            static_cast<int>(std::round(polytope.seed.x())),
            static_cast<int>(std::round(polytope.seed.y()))
        };

        line_eq::LineEquations le(hull_points, seed_ij);
        le.compute();
        MatrixXd A = le.A();
        VectorXd b = le.b();
        
        // Convert to halfspaces ax+by+c<=0
        Vector2d seed(poly_x[0], poly_y[0]);
        for(int i = 0; i < A.rows(); i++) {
            double lhs = A.row(i).dot(seed);
            if(lhs > b(i)) { 
                A.row(i) *= -1; 
                b(i) *= -1; 
            }
            polytope.appendPlane(Vector3d(A(i,0), A(i,1), -b(i)));
        }
        
        // Center
        Vector2d center = Vector2d::Zero();
        for(const auto& v : hull.vertices()) {
            center += Vector2d(v.x, v.y);
        }
        center /= hull.vertices().size();
        polytope.setCenter(center);
        polytope.setSeed(seed);

        return polytope;
    }



    void printInsideData() {
        // cout << "\nInside Data (" << _max_x_id << "x" << _max_y_id << " grid):" << endl;
        // cout << "Legend: 0=boundary/free, 1=invalid, 2=inside polygon" << endl;
        for(int i = 0; i < _max_x_id; i++) {
            for(int j = 0; j < _max_y_id; j++) {
                if(inside_data(i, j) != 2) {
                    inside_data(i, j) = map_data(i, j);
                }
                // cout << inside_data(i, j) << " ";
            }
            cout << endl;
        }
    }

    void printMapData() {
        // cout << "\nMap Data (" << _max_x_id << "x" << _max_y_id << " grid - 1=obstacle, 0=free):" << endl;
        // for(int i = 0; i < _max_x_id; i++) {20
        //     for(int j = 0; j < _max_y_id; j++) {
        //         cout << map_data(i, j) << " ";
        //     }
        //     cout << endl;
        // }
    }

    void printVertexIdx(const VectorXi& vertex_idx, const string& name) {
        // cout << "\n" << name << ":" << endl;
        // cout << "X indices: [" << vertex_idx(0) << ", " << vertex_idx(1) << ", " 
        //      << vertex_idx(2) << ", " << vertex_idx(3) << "]" << endl;
        // cout << "Y indices: [" << vertex_idx(4) << ", " << vertex_idx(5) << ", " 
        //      << vertex_idx(6) << ", " << vertex_idx(7) << "]" << endl;
        // cout << "Bounding box: X[" << vertex_idx(2) << "," << vertex_idx(0) 
        //      << "] Y[" << vertex_idx(4) << "," << vertex_idx(5) << "]" << endl;
    }


 
    void loadMapData(const vector<uint8_t>& data) {
        if(data.size() != (size_t)_grid_num) {
            cerr << "Error: Map data size mismatch!" << endl;
            return;
        }

        // 1. Clear previous obstacles
        _occupied_cells.clear(); 
        // Optional: reserve memory to prevent re-allocations if you know avg density
        // _occupied_cells.reserve(_grid_num / 10); 

        for(int j = 0; j < _max_y_id; j++) {
            for(int i = 0; i < _max_x_id; i++) {
                uint8_t val = data[j * _max_x_id + i];
                
                // Store to grid
                map_data(i, j) = val;

                // 2. Store Occupied Coordinate
                // Assuming val > 0 means obstacle (adjust if 0 is obstacle in your map)
                if(val > 0) {
                    _occupied_cells.emplace_back(i, j);
                }
            }
        }
    }

    MatrixXi getInsideData() {
        return inside_data;
    }

    void printPolygonResult(const vector<int>& cluster_x, const vector<int>& cluster_y) {
        // cout << "\n=== FINAL POLYGON RESULT ===" << endl;
        // cout << "Total vertices: " << cluster_x.size() << endl;
        
        // Create visualization grid
        MatrixXi viz_grid = map_data;
        
        // Mark polygon vertices with 'P' (3)
        for(size_t i = 0; i < cluster_x.size(); i++) {
            viz_grid(cluster_x[i], cluster_y[i]) = 3;
        }
        
        cout << "\nVisualization (0=free, 1=obstacle, 2=inflated_cube, 3=polygon_vertex):" << endl;
        for(int i = 0; i < _max_x_id; i++) {
             for(int j = 0; j < _max_y_id; j++) {
                 if(inside_data(i, j) == 2) {
                     cout << "P ";
                 } else if(viz_grid(i, j) == 3) {
                     cout << "P ";
                 } else {
                     cout << viz_grid(i, j) << " ";
                 }
             }
             cout << endl;
         }
        
        // cout << "\nPolygon vertices (first 20):" << endl;
        for(size_t i = 0; i < min(cluster_x.size(), size_t(20)); i++) {
            // cout << "  (" << cluster_x[i] << ", " << cluster_y[i] << ")";
            // if((i + 1) % 5 == 0) cout << endl;
        }
        if(cluster_x.size() > 20) {
            // cout << "\n  ... and " << (cluster_x.size() - 20) << " more vertices" << endl;
        }
    }

    // ============================
    // ACCESSORS
    // ============================
    
    const std::vector<Polytope>& getCorridor() const {
        return polytope_corridor;
    }

    std::vector<Polytope>& getCorridorMutable() {
        return polytope_corridor;
    }

    void clearCorridor() {
        polytope_corridor.clear();
    }

    void addPolytope(const Polytope& P) {
        polytope_corridor.push_back(P);
    }

    // ... your existing public variables and functions ..

private:
    int _max_x_id;
    int _max_y_id;
    int _grid_num;
    
    int itr_inflate_max;
    int itr_cluster_max;
    int inf_step;
    int _cluster_buffer_size;
    int _candidate_buffer_size;
    int _cluster_buffer_size_square;
    
    // Eigen data structures
    MatrixXi map_data;
    MatrixXi use_data;
    MatrixXi invalid_data;
    MatrixXi inside_data;
    VectorXi vertex_idx;
    VectorXi vertex_idx_lst;

    std::vector<Eigen::Vector2i> _occupied_cells;               // to store the indexes of obstacles
    
    // Cluster data structures
    MatrixXi active_xy_id;
    MatrixXi cluster_xy_id;
    MatrixXi candidate_cell;
    //Added
    // BD bounding box limits — constrain inflation to local region
    int limit_min_x ;
    int limit_max_x ;
    int limit_min_y ;
    int limit_max_y ;

    


    inline float signum_cpu(const int &x)
    {
        return x == 0 ? 0 : x < 0 ? -1.0f : 1.0f;
    }

    inline float mod_cpu(const float &value, const float &modulus)
    {
        return fmod(fmod(value, modulus) + modulus, modulus);
    }

    inline float intbound_cpu(float s, int ds)
    {
        if (ds == 0)
        {
            return numeric_limits<float>::max();
        }
        else if (ds < 0)
        {
            return intbound_cpu(-s, -ds);
        }
        else
        {
            s = mod_cpu(s, 1.0f);
            return (1 - s) / ds;
        }
    }
};

#endif // POLYHEDRON_GENERATOR_CLASS_HPP