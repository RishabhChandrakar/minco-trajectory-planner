#ifndef VEHICLE_UTILS_HPP
#define VEHICLE_UTILS_HPP

#include <string>

namespace vehicle_utils
{

// ---------- 3D Vector ----------
struct Vector3
{
    double x{0.0};
    double y{0.0};
    double z{0.0};
};

// ---------- Euler Orientation ----------
struct Euler
{
    double roll{0.0};
    double pitch{0.0};
    double yaw{0.0};
};

// ---------- Quaternion ----------
struct Quaternion
{
    double x{0.0};
    double y{0.0};
    double z{0.0};
    double w{1.0};
};

// ---------- Angular Velocity ----------
struct AngularVelocity
{
    double roll_rate{0.0};
    double pitch_rate{0.0};
    double yaw_rate{0.0};
};

// ---------- Vehicle Odometry ----------
struct VehicleOdometry
{
    Vector3 position;
    Vector3 velocity;

    Quaternion quaternion;
    Euler orientation;

    AngularVelocity angular_velocity;
};

// -------- WAYPOINT ------
struct Waypoint {
        double x, y, z;
    };

} // namespace vehicle_utils

#endif
