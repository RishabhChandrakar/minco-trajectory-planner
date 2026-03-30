import os
from ament_index_python.packages import get_package_share_directory
from launch import LaunchDescription
from launch_ros.actions import Node

def generate_launch_description():
    # Get the package directory
    pkg_dir = get_package_share_directory('planner')

    # Define the path to the parameter file
    params_file = os.path.join(pkg_dir, 'config', 'mapping_tracking_params.yaml')

    return LaunchDescription([
        Node(
            package='planner',
            executable='integrated_drone_mapping_node_with_drone_tracking',
            name='integrated_mapping_sfc_gcopter',
            output='screen',
            parameters=[params_file]
        )
    ])
