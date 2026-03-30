#!/usr/bin/env python3

import rclpy
from rclpy.node import Node
from geometry_msgs.msg import TwistStamped

import matplotlib.pyplot as plt
import os

class PlottingNode(Node):
    def __init__(self, namespace="uav1"):
        super().__init__('tracking_error_plotter')

        self.times = []
        
        # Velocity Errors (linear)
        self.v_err_x = []
        self.v_err_y = []
        self.v_err_z = []
        
        # Position Errors (angular)
        self.p_err_x = []
        self.p_err_y = []
        self.p_err_z = []

        self.start_time = None
        self.msg_count = 0

        topic_name = f'/{namespace}/tracking_error'
        self.get_logger().info(f"Subscribing to {topic_name}. Press Ctrl+C to stop recording and generate the plot.")

        self.subscription = self.create_subscription(
            TwistStamped,
            topic_name,
            self.error_callback,
            100)

    def error_callback(self, msg):
        current_time = msg.header.stamp.sec + msg.header.stamp.nanosec * 1e-9
        
        if self.start_time is None:
            self.start_time = current_time

        t = current_time - self.start_time

        self.times.append(t)
        
        # linear contains velocity errors
        self.v_err_x.append(msg.twist.linear.x)
        self.v_err_y.append(msg.twist.linear.y)
        self.v_err_z.append(msg.twist.linear.z)

        # angular contains position errors
        self.p_err_x.append(msg.twist.angular.x)
        self.p_err_y.append(msg.twist.angular.y)
        self.p_err_z.append(msg.twist.angular.z)

        self.msg_count += 1
        if self.msg_count % 100 == 0:
            self.get_logger().info(f"Recorded {self.msg_count} messages...")

def main(args=None):
    rclpy.init(args=args)
    node = PlottingNode()

    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        node.get_logger().info("KeyboardInterrupt received. Stopping recording and generating plot.")

    if len(node.times) > 0:
        # Generate the plot
        import matplotlib
        import datetime
        matplotlib.use('Agg')  # headless generation
        
        fig, (ax_pos, ax_vel) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        fig.suptitle("Drone Tracking Error Tuning (KP / KV)", fontsize=16)

        # Position Error Plot
        ax_pos.set_title("Position Error (m)")
        ax_pos.set_ylabel("Error (m)")
        ax_pos.grid(True)
        ax_pos.plot(node.times, node.p_err_x, label="X err", color='r', alpha=0.8)
        ax_pos.plot(node.times, node.p_err_y, label="Y err", color='g', alpha=0.8)
        ax_pos.plot(node.times, node.p_err_z, label="Z err", color='b', alpha=0.8)
        ax_pos.legend(loc="upper right")

        # Velocity Error Plot
        ax_vel.set_title("Velocity Error (m/s)")
        ax_vel.set_xlabel("Time (s)")
        ax_vel.set_ylabel("Error (m/s)")
        ax_vel.grid(True)
        ax_vel.plot(node.times, node.v_err_x, label="Vx err", color='r', alpha=0.8)
        ax_vel.plot(node.times, node.v_err_y, label="Vy err", color='g', alpha=0.8)
        ax_vel.plot(node.times, node.v_err_z, label="Vz err", color='b', alpha=0.8)
        ax_vel.legend(loc="upper right")

        plt.tight_layout()

        # Save to PNG
        timestamp_str = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"tracking_error_plot_{timestamp_str}.png"
        filepath = os.path.join(os.getcwd(), filename)
        
        fig.savefig(filepath, dpi=300)
        node.get_logger().info(f"Successfully saved tracking error plot to: {filepath}")
    else:
        node.get_logger().warn("No data recorded. The plot was not generated.")

    node.destroy_node()
    rclpy.shutdown()

if __name__ == '__main__':
    main()
