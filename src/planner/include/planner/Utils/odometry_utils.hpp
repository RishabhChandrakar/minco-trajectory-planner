#ifndef ODOMETRY_UTILS_HPP
#define ODOMETRY_UTILS_HPP

#include <cmath>
#include "vehicle_utils.hpp"

namespace odom_utils
{using namespace vehicle_utils;

inline void quaternionToEuler(
    const vehicle_utils::Quaternion &q,
    Euler &euler)
{
    double qx = q.x;
    double qy = q.y;
    double qz = q.z;
    double qw = q.w;

    // Roll
    double sinr_cosp = 2.0 * (qw * qx + qy * qz);
    double cosr_cosp = 1.0 - 2.0 * (qx * qx + qy * qy);
    euler.roll = std::atan2(sinr_cosp, cosr_cosp);

    // Pitch
    double sinp = 2.0 * (qw * qy - qz * qx);
    if (std::abs(sinp) >= 1)
        euler.pitch = std::copysign(M_PI / 2, sinp);
    else
        euler.pitch = std::asin(sinp);

    // Yaw
    double siny_cosp = 2.0 * (qw * qz + qx * qy);
    double cosy_cosp = 1.0 - 2.0 * (qy * qy + qz * qz);
    euler.yaw = std::atan2(siny_cosp, cosy_cosp);
}

} // namespace odom_utils

#endif
