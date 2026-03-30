import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# =========================
# Parameters (match C++)
# =========================
GRID_W = 50.0      # meters
GRID_H = 50.0
RES = 0.05         # resolution

# =========================
# Load Map
# =========================
import numpy as np

def load_map_csv(filename):
    rows = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # remove trailing comma safely
            if line.endswith(","):
                line = line[:-1]
            row = [int(x) for x in line.split(",") if x != ""]
            rows.append(row)
    return np.array(rows)

map_data = load_map_csv("map.csv")
map_data = map_data.astype(int)

grid_h, grid_w = map_data.shape

# World coordinate extents
x_min = -GRID_W / 2
x_max = GRID_W / 2
y_min = -GRID_H / 2
y_max = GRID_H / 2

# =========================
# Load Executed Trajectory
# =========================
executed = pd.read_csv("executed_path.csv")
exec_x = executed["x"].values
exec_y = executed["y"].values

# =========================
# Load Planned Trajectories (supports multiple replans)
# =========================
planned_trajs = []
current_traj = []

with open("planned_traj.csv", "r") as f:
    for line in f:
        line = line.strip()

        if not line or "REPLAN_ID" in line:
            if current_traj:
                planned_trajs.append(np.array(current_traj))
                current_traj = []
            continue

        if line.startswith("t,") or line.startswith("x,"):
            continue

        parts = line.split(",")
        if len(parts) >= 4:
            # If format contains time column
            if len(parts) >= 10:
                x = float(parts[1])
                y = float(parts[2])
            else:
                x = float(parts[0])
                y = float(parts[1])

            current_traj.append([x, y])

    if current_traj:
        planned_trajs.append(np.array(current_traj))

# =========================
# Plot
# =========================
# =========================
plt.figure(figsize=(10, 10))

# Plot occupancy map
plt.imshow(
    map_data,
    origin="lower",
    extent=[x_min, x_max, y_min, y_max],
    cmap="gray_r",
)

# Plot planned trajectories
for i, traj in enumerate(planned_trajs):
    if len(traj) > 0:
        plt.plot(
            traj[:, 0],
            traj[:, 1],
            linestyle="--",
            linewidth=2,
            label=f"Planned Replan {i+1}",
        )

# Plot executed trajectory
if len(exec_x) > 0:
    plt.plot(
        exec_x,
        exec_y,
        linewidth=3,
        label="Executed Path",
    )

plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.title("Map + Planned Trajectories + Executed Path")
plt.axis("equal")
plt.grid(True)
plt.legend()

# Save to PNG (no GUI window)
plt.savefig("map_with_trajectories.png", dpi=300, bbox_inches="tight")
plt.close()