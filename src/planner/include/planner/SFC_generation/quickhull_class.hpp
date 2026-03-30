#ifndef QUICKHULL_2D_HPP_
#define QUICKHULL_2D_HPP_
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <cassert>
#include <limits>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <cstddef>
#include <memory>
#include <iostream>

namespace sfc_quickhull {

using IndexType = std::size_t;

// ============================================================
//                    2D VECTOR CLASS
// ============================================================

template<typename T>
struct Vector2 {
    T x, y;

    Vector2() : x(0), y(0) {}
    Vector2(T X, T Y) : x(X), y(Y) {}

    bool operator==(const Vector2& v) const { return x == v.x && y == v.y; }
    bool operator!=(const Vector2& v) const { return !(*this == v); }
};

// For debugging
template<typename T>
std::ostream& operator<<(std::ostream& os, const Vector2<T>& v) {
    return (os << "(" << v.x << "," << v.y << ")");
}

// ============================================================
//                 VERTEX DATA SOURCE WRAPPER
// ============================================================

template<typename T>
class VertexDataSource {
    const Vector2<T>* m_ptr;
    size_t m_count;

public:
    VertexDataSource() : m_ptr(nullptr), m_count(0) {}

    VertexDataSource(const Vector2<T>* ptr, size_t count)
        : m_ptr(ptr), m_count(count) {}

    VertexDataSource(const std::vector<Vector2<T>>& vec)
        : m_ptr(vec.data()), m_count(vec.size()) {}

    size_t size() const { return m_count; }

    const Vector2<T>& operator[](size_t i) const { return m_ptr[i]; }

    const Vector2<T>* begin() const { return m_ptr; }
    const Vector2<T>* end()   const { return m_ptr + m_count; }
};

// ============================================================
//                RESULT CONVEX HULL STRUCTURE
// ============================================================

template<typename T>
class ConvexHull {
public:
    // Either references the original buffer, or our optimized one.
    VertexDataSource<T> m_vertices;

    // Either hull indices into original data or sequential {0..H-1}
    std::vector<IndexType> m_indices;

private:
    std::unique_ptr<std::vector<Vector2<T>>> m_internalBuffer;

public:
    ConvexHull() = default;

    const std::vector<IndexType>& indices() const { return m_indices; }
    const VertexDataSource<T>&   vertices() const { return m_vertices; }

    // Export OBJ (for visualization)
    void writeOBJ(const std::string& filename) const {
        std::ofstream f(filename);
        if(!f.is_open()) return;

        for(size_t i=0;i<m_vertices.size();i++)
            f << "v " << m_vertices[i].x << " " << m_vertices[i].y << " 0\n";

        for(size_t i=0;i<m_indices.size();i++){
            size_t j = (i+1)%m_indices.size();
            f << "l " << m_indices[i]+1 << " " << m_indices[j]+1 << "\n";
        }
    }

    template<typename U> friend class QuickHull;
};

// ============================================================
//                     QUICKHULL 2D (MONOTONE CHAIN)
// ============================================================

template<typename T>
class QuickHull {
public:
    static constexpr T EPSILON = (std::is_same<T,double>::value ? T(1e-12) : T(1e-6));

private:

    // Cross product (2D version)
    static inline T cross(const Vector2<T>& O, const Vector2<T>& A, const Vector2<T>& B) {
        return (A.x - O.x)*(B.y - O.y) - (A.y - O.y)*(B.x - O.x);
    }

    // Actual convex hull computation
    std::vector<IndexType> computeHull(const VertexDataSource<T>& pts) const {
        size_t n = pts.size();
        std::vector<IndexType> H;

        if(n == 0) return H;
        if(n == 1) { H.push_back(0); return H; }

        // 1) Sort lexicographically by (x,y)
        std::vector<IndexType> order(n);
        for(size_t i=0;i<n;i++) order[i]=i;

        std::sort(order.begin(), order.end(), [&](IndexType a, IndexType b){
            if (std::fabs(pts[a].x - pts[b].x) > EPSILON)
                return pts[a].x < pts[b].x;
            return pts[a].y < pts[b].y;
        });

        // 2) Remove duplicate points
        std::vector<IndexType> uniq;
        uniq.reserve(n);
        uniq.push_back(order[0]);

        for(size_t i=1;i<n;i++){
            const auto& p = pts[order[i]];
            const auto& q = pts[uniq.back()];
            if(std::fabs(p.x-q.x)>EPSILON || std::fabs(p.y-q.y)>EPSILON)
                uniq.push_back(order[i]);
        }

        if(uniq.size() == 1){
            return uniq;
        }

        // 3) Build lower hull
        std::vector<IndexType> lower;
        for(IndexType idx : uniq) {
            while(lower.size() >= 2 &&
                  cross(pts[lower[lower.size()-2]], pts[lower.back()], pts[idx]) <= 0)
                lower.pop_back();
            lower.push_back(idx);
        }

        // 4) Build upper hull
        std::vector<IndexType> upper;
        for(int i = (int)uniq.size()-1; i >= 0; i--) {
            IndexType idx = uniq[i];
            while(upper.size() >= 2 &&
                  cross(pts[upper[upper.size()-2]], pts[upper.back()], pts[idx]) <= 0)
                upper.pop_back();
            upper.push_back(idx);
        }

        // 5) Remove duplicates and combine
        lower.pop_back();
        upper.pop_back();

        H.reserve(lower.size() + upper.size());
        H.insert(H.end(), lower.begin(), lower.end());
        H.insert(H.end(), upper.begin(), upper.end());

        return H;   // CCW
    }

public:

    // ------------------------------------------------------------
    //    MAIN API #1: From std::vector<Vector2<T>>
    // ------------------------------------------------------------
    ConvexHull<T> getConvexHull(const std::vector<Vector2<T>>& cloud,
                                bool CCW = true,
                                bool useOriginalIndices = false)
    {
        VertexDataSource<T> pts(cloud);
        auto hullIdx = computeHull(pts);

        ConvexHull<T> out;

        if(useOriginalIndices){
            out.m_vertices = pts;
            out.m_indices  = hullIdx;
        }
        else {
            // compact vertex buffer
            out.m_internalBuffer = std::make_unique<std::vector<Vector2<T>>>();
            for(auto id : hullIdx)
                out.m_internalBuffer->push_back(pts[id]);

            out.m_vertices = VertexDataSource<T>(*out.m_internalBuffer);
            out.m_indices.resize(hullIdx.size());
            for(size_t i=0;i<hullIdx.size();i++)
                out.m_indices[i]=i;
        }

        if(!CCW)
            std::reverse(out.m_indices.begin(), out.m_indices.end());

        return out;
    }

    // ------------------------------------------------------------
    //      MAIN API #2: From raw Vector2<T> pointer
    // ------------------------------------------------------------
    ConvexHull<T> getConvexHull(const Vector2<T>* ptsArr,
                                size_t count,
                                bool CCW=true,
                                bool useOriginalIndices=false)
    {
        VertexDataSource<T> pts(ptsArr, count);
        auto hullIdx = computeHull(pts);

        ConvexHull<T> out;

        if(useOriginalIndices){
            out.m_vertices = pts;
            out.m_indices = hullIdx;
        }
        else {
            out.m_internalBuffer = std::make_unique<std::vector<Vector2<T>>>();
            for(auto id : hullIdx)
                out.m_internalBuffer->push_back(pts[id]);

            out.m_vertices = VertexDataSource<T>(*out.m_internalBuffer);

            out.m_indices.resize(hullIdx.size());
            for(size_t i=0;i<hullIdx.size();i++)
                out.m_indices[i]=i;
        }

        if(!CCW)
            std::reverse(out.m_indices.begin(), out.m_indices.end());

        return out;
    }

    // ------------------------------------------------------------
    //      MAIN API #3: From flat coordinate array [x0,y0,x1,y1,...]
    // ------------------------------------------------------------
    ConvexHull<T> getConvexHull(const T* flatXY,
                                size_t count,
                                bool CCW=true,
                                bool useOriginalIndices=false)
    {
        std::vector<Vector2<T>> pts;
        pts.reserve(count);

        for(size_t i=0;i<count;i++)
            pts.emplace_back(flatXY[2*i], flatXY[2*i+1]);

        return getConvexHull(pts, CCW, useOriginalIndices);
    }
};

using QuickHullf = QuickHull<float>;
using QuickHulld = QuickHull<double>;

} // namespace sfc_quickhull

#endif // QUICKHULL_2D_HPP_