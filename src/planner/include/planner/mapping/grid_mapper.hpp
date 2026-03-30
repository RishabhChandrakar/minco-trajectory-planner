#ifndef OCCUPANCY_MAPPER_HPP
#define OCCUPANCY_MAPPER_HPP

// Necessary includes for ROS messages
#include <rclcpp/rclcpp.hpp>
#include <nav_msgs/msg/occupancy_grid.hpp>
#include <geometry_msgs/msg/pose.hpp>

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

class OccupancyMapper {
public:
    /**
     * @brief Constructor
     * @param width_m Map width in meters (e.g., 20.0)
     * @param height_m Map height in meters (e.g., 20.0)
     * @param resolution Map resolution (e.g., 0.1)
     * @param sensor_max_range LiDAR Max range to clamp INF values (e.g., 20.0)
     */
    OccupancyMapper(float width_m, float height_m, float resolution, float sensor_max_range) 
        : width_m_(width_m), 
          height_m_(height_m), 
          resolution_(resolution), 
          max_range_(sensor_max_range) 
    {
        // Grid Dimensions
        width_cells_ = static_cast<int>(width_m_ / resolution_);
        height_cells_ = static_cast<int>(height_m_ / resolution_);
        
        // Center Origin
        origin_x_cell_ = width_cells_ / 2;
        origin_y_cell_ = height_cells_ / 2;

        // Initialize with 0.0 log odds
        log_odds_grid_.resize(width_cells_ * height_cells_, 0.0f);
    }

    /**
     * @brief Updates the map based on scan data
     */
    void update_map(float robot_x, float robot_y, float robot_theta, 
                    const std::vector<float>& ranges, 
                    float angle_min, float angle_increment) {
        
        int robot_gx, robot_gy;
        world_to_grid(robot_x, robot_y, robot_gx, robot_gy);

        float current_angle = angle_min;

        for (float r : ranges) {
            float range_m;
            
            // Logic: Handle INF by capping at max_range
            if (std::isinf(r) || std::isnan(r)) {
                range_m = max_range_;
            } else {
                range_m = r;
            }

            // Logic: Calculate Hit Point
            // NOTE: Using 'minus sin' to strictly match your Python logic

            // since we are doing it in ENU frame 

            float global_angle = robot_theta + current_angle;
            float hit_x = robot_x + (range_m * std::cos(global_angle));
            float hit_y = robot_y + (range_m * std::sin(global_angle));

            int hit_gx, hit_gy;
            world_to_grid(hit_x, hit_y, hit_gx, hit_gy);

            // Bresenham Raycasting
            int x0 = robot_gx; int y0 = robot_gy;
            int x1 = hit_gx;   int y1 = hit_gy;
            int dx = std::abs(x1 - x0);
            int dy = -std::abs(y1 - y0);
            int sx = (x0 < x1) ? 1 : -1;
            int sy = (y0 < y1) ? 1 : -1;
            int err = dx + dy;

            while (true) {
                bool is_last_point = (x0 == x1 && y0 == y1);
                int idx = y0 * width_cells_ + x0;

                if (is_last_point) {
                    // Update Occupied
                    if (is_in_grid_bounds(x0, y0)) {
                        log_odds_grid_[idx] += LOG_ODDS_OCC;
                        clamp_value(idx);
                    }
                    break;
                } else {
                    // Update Free
                    if (is_in_grid_bounds(x0, y0)) {
                        log_odds_grid_[idx] += LOG_ODDS_FREE;
                        clamp_value(idx);
                    }
                }

                int e2 = 2 * err;
                if (e2 >= dy) { err += dy; x0 += sx; }
                if (e2 <= dx) { err += dx; y0 += sy; }
            }
            current_angle += angle_increment;
        }
    }

    /**
     * @brief Generates the final ROS OccupancyGrid message
     * @param frame_id The frame to publish in (e.g., "map" or "odom")
     * @param stamp The timestamp for the message header
     * @return Fully populated nav_msgs::msg::OccupancyGrid
     */
    nav_msgs::msg::OccupancyGrid get_grid_message(const std::string& frame_id, const rclcpp::Time& stamp) const {
        nav_msgs::msg::OccupancyGrid msg;

        // 1. Fill Header
        msg.header.frame_id = frame_id;
        msg.header.stamp = stamp;

        // 2. Fill MetaData
        msg.info.resolution = resolution_;
        msg.info.width = width_cells_;
        msg.info.height = height_cells_;
        
        // Origin Calculation (Bottom-Left corner relative to center)
        msg.info.origin.position.x = -width_m_ / 2.0f;
        msg.info.origin.position.y = -height_m_ / 2.0f;
        msg.info.origin.position.z = 0.0f;
        msg.info.origin.orientation.w = 1.0f; // Identity quaternion

        // 3. Fill Data
        // Resize payload to match grid size
        msg.data.resize(log_odds_grid_.size());

        for (size_t i = 0; i < log_odds_grid_.size(); ++i) {
            float log_val = log_odds_grid_[i];
            
            // Convert Log Odds to Probability (0.0 to 1.0)
            // Formula matches Python: 1 - (1 / (1 + exp(val)))
            float probability = 1.0f - (1.0f / (1.0f + std::exp(log_val)));

            // Apply Python thresholds
            if (probability < 0.2f) {
                msg.data[i] = 0;   // Free
            } else if (probability > 0.8f) {
                msg.data[i] = 100; // Occupied (ROS standard for wall)
            } else {
                msg.data[i] = -1;  // Unknown
            }
        }

        return msg;
    }

private:
    float width_m_;
    float height_m_;
    float resolution_;
    float max_range_;
    
    int width_cells_;
    int height_cells_;
    int origin_x_cell_;
    int origin_y_cell_;

    // Python Logic Parameters
    const float LOG_ODDS_OCC = 2.2f;
    const float LOG_ODDS_FREE = -0.7f;
    const float LOG_ODDS_MAX = 5.0f;
    const float LOG_ODDS_MIN = -5.0f;

    std::vector<float> log_odds_grid_;

    void world_to_grid(float wx, float wy, int& gx, int& gy) {
        gx = static_cast<int>((wx / resolution_) + origin_x_cell_);
        gy = static_cast<int>((wy / resolution_) + origin_y_cell_);
    }

    bool is_in_grid_bounds(int gx, int gy) {
        return (gx >= 0 && gx < width_cells_ && gy >= 0 && gy < height_cells_);
    }

    void clamp_value(int idx) {
        if (log_odds_grid_[idx] > LOG_ODDS_MAX) log_odds_grid_[idx] = LOG_ODDS_MAX;
        else if (log_odds_grid_[idx] < LOG_ODDS_MIN) log_odds_grid_[idx] = LOG_ODDS_MIN;
    }
};

#endif // OCCUPANCY_MAPPER_HPP