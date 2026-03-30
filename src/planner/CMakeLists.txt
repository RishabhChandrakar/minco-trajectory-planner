cmake_minimum_required(VERSION 3.8)
project(planner)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Core dependencies
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(nav_msgs REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(std_msgs REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(mavros_msgs REQUIRED)

# Include directories (package + shared Utils)
include_directories(
  include
)

if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()

##################################################
# Integrated drone mapping node with drone tracking (MAVROS)
##################################################

add_executable(integrated_drone_mapping_node_with_drone_tracking
  src/integrated_drone_mapping_node_with_drone_tracking.cpp
)

ament_target_dependencies(integrated_drone_mapping_node_with_drone_tracking
  rclcpp
  sensor_msgs
  nav_msgs
  geometry_msgs
  mavros_msgs
  std_msgs
)

target_link_libraries(integrated_drone_mapping_node_with_drone_tracking
  Eigen3::Eigen
)

target_compile_options(integrated_drone_mapping_node_with_drone_tracking PRIVATE
  -O3 -march=native -DNDEBUG
)

##################################################
# Install
##################################################

install(TARGETS
  integrated_drone_mapping_node_with_drone_tracking
  DESTINATION lib/${PROJECT_NAME}
)

install(DIRECTORY
  launch
  config
  DESTINATION share/${PROJECT_NAME}
)

ament_package()
