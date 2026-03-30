#ifndef OBSTACLE_INFLATOR_HPP
#define OBSTACLE_INFLATOR_HPP

#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>

class ObstacleInflator {
public:
    /**
     * @brief Inflate obstacles in a 2D grid map based on drone radius
     *
     * @param map_flat               Row-major occupancy grid (0=free, 1=occupied)
     * @param grid_width             Grid width (cells)
     * @param grid_height            Grid height (cells)
     * @param resolution_m_per_cell  Cell resolution in meters
     * @param droneRadiusCm          Drone radius in centimeters
     *
     * @return Inflated occupancy grid (same layout, same size)
     */
    static std::vector<uint8_t> inflateObstacles(
        const std::vector<uint8_t>& map_flat,
        int grid_width,
        int grid_height,
        double resolution_m_per_cell,
        double droneRadiusCm
    ) {
        if (grid_width <= 0 || grid_height <= 0 ||
            map_flat.size() != static_cast<size_t>(grid_width * grid_height)) {
            return map_flat;
        }

        // Convert drone radius from cm → meters → grid cells
        double droneRadiusM = droneRadiusCm / 100.0;
        int radiusCells = static_cast<int>(
            std::ceil(droneRadiusM / resolution_m_per_cell)
        );
        int radiusSq = radiusCells * radiusCells;

        std::vector<uint8_t> safeMap(grid_width * grid_height, 0);

        // ---- 1) Copy original obstacles ----
        for (int y = 0; y < grid_height; ++y) {
            for (int x = 0; x < grid_width; ++x) {
                int idx = y * grid_width + x;
                if (map_flat[idx] == 1) {
                    safeMap[idx] = 1;
                }
            }
        }

        // ---- 2) Inflate obstacles ----
        for (int y = 0; y < grid_height; ++y) {
            for (int x = 0; x < grid_width; ++x) {
                int idx = y * grid_width + x;
                if (safeMap[idx] == 1) continue;

                int startY = std::max(0, y - radiusCells);
                int endY   = std::min(grid_height - 1, y + radiusCells);
                int startX = std::max(0, x - radiusCells);
                int endX   = std::min(grid_width - 1, x + radiusCells);

                bool tooClose = false;

                for (int yy = startY; yy <= endY && !tooClose; ++yy) {
                    for (int xx = startX; xx <= endX; ++xx) {
                        int obs_idx = yy * grid_width + xx;
                        if (map_flat[obs_idx] == 1) {
                            int dx = yy - y;
                            int dy = xx - x;
                            if (dx * dx + dy * dy <= radiusSq) {
                                tooClose = true;
                                break;
                            }
                        }
                    }
                }

                if (tooClose) {
                    safeMap[idx] = 1;
                }
            }
        }

        return safeMap;
    }
};

#endif // OBSTACLE_INFLATOR_HPP