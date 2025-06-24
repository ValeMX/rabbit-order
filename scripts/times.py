import pandas as pd
import numpy as np
import os
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt

# Ask user for input file
input_file = input("Enter the CSV file name (e.g., out.csv): ").strip()

# Load the CSV
df = pd.read_csv(input_file)
output_dir = os.path.dirname(os.path.abspath(input_file))

# Create output directories if they do not exist
for subdir in ["general_edges", "general_nodes", "specific_edges", "specific_nodes"]:
    dir_path = os.path.join(output_dir, subdir)
    os.makedirs(dir_path, exist_ok=True)

# Ensure numeric fields are of the correct type
df["size"] = df["size"].astype(int)
df["totalTime"] = df["totalTime"].astype(float)

# ----------- PER-SIZE DETAILED TIME PLOTS -----------
time_columns = [
    "partitionTime",
    "distributionTime",
    "collectingRemoteTime",
    "detectionTime",
    "exchangingCommunitiesTime",
    "collectingMissingTime",
    "coarsenTime",
]

# Ensure all time columns exist and are float
for col in time_columns:
    if col in df.columns:
        df[col] = df[col].astype(float)
    else:
        df[col] = np.nan  # Fill missing columns with NaN

# Compute mean of all time_columns grouped only by totalEdges and size
df_time_mean = df.groupby(["totalEdges", "size"], as_index=False)[time_columns].mean()

# ----------- PLOT TIME BREAKDOWN FOR EACH SIZE BY EDGES -----------
for size in sorted(df_time_mean["size"].unique()):
    plt.figure(figsize=(10, 7))
    df_size = df_time_mean[df_time_mean["size"] == size].copy()
    df_size = df_size.sort_values("totalEdges")
    x = df_size["totalEdges"].values
    for col in time_columns:
        if df_size[col].isnull().all():
            continue
        # Convert seconds to milliseconds for plotting
        plt.plot(x, df_size[col] * 1000, marker="o", linestyle="-", label=col)
    plt.xlabel("Total Edges")
    plt.ylabel("Time (ms)")
    plt.xscale("log", base=10)
    plt.yscale("log", base=10)
    plt.grid(True, which="both", ls="--")
    plt.legend(fontsize="small", loc="best")
    plt.tight_layout()
    plt.subplots_adjust(left=0.1)
    plt.gca().yaxis.set_major_formatter(ScalarFormatter())
    plt.gca().ticklabel_format(style="plain", axis="y")
    plt.savefig(
        os.path.join(
            output_dir, "general_edges", f"time_breakdown_edges_size_{size}.png"
        )
    )
    plt.close()

df_time_mean = df.groupby(["totalNodes", "size"], as_index=False)[time_columns].mean()

# ----------- PLOT TIME BREAKDOWN FOR EACH SIZE BY NODES -----------
for size in sorted(df_time_mean["size"].unique()):
    plt.figure(figsize=(10, 7))
    df_size = df_time_mean[df_time_mean["size"] == size].copy()
    df_size = df_size.sort_values("totalNodes")
    x = df_size["totalNodes"].values
    for col in time_columns:
        if df_size[col].isnull().all():
            continue
        # Convert seconds to milliseconds for plotting
        plt.plot(x, df_size[col] * 1000, marker="o", linestyle="-", label=col)
    plt.xlabel("Total Nodes")
    plt.ylabel("Time (ms)")
    plt.xscale("log", base=10)
    plt.yscale("log", base=10)
    plt.grid(True, which="both", ls="--")
    plt.legend(fontsize="small", loc="best")
    plt.tight_layout()
    plt.subplots_adjust(left=0.1)
    plt.gca().yaxis.set_major_formatter(ScalarFormatter())
    plt.gca().ticklabel_format(style="plain", axis="y")
    plt.savefig(
        os.path.join(
            output_dir, "general_nodes", f"time_breakdown_nodes_size_{size}.png"
        )
    )
    plt.close()

# ----------- PLOT SPECIFIC TIME BREAKDOWN FOR EACH SIZE BY EDGES -----------
time_columns = [
    "partitionTime",
]

df_time_mean = df.groupby(["totalEdges", "size"], as_index=False)[time_columns].mean()

for size in sorted(df_time_mean["size"].unique()):
    plt.figure(figsize=(10, 7))
    df_size = df_time_mean[df_time_mean["size"] == size].copy()
    df_size = df_size.sort_values("totalEdges")
    x = df_size["totalEdges"].values
    for col in time_columns:
        if df_size[col].isnull().all():
            continue
        # Convert seconds to milliseconds for plotting
        plt.plot(x, df_size[col] * 1000, marker="o", linestyle="-", label=col)
    plt.xlabel("Total Edges")
    plt.ylabel("Time (ms)")
    plt.xscale("log", base=10)
    plt.yscale("log", base=10)
    plt.grid(True, which="both", ls="--")
    plt.legend(fontsize="small", loc="best")
    plt.tight_layout()
    plt.subplots_adjust(left=0.1)
    plt.gca().yaxis.set_major_formatter(ScalarFormatter())
    plt.gca().ticklabel_format(style="plain", axis="y")
    plt.savefig(
        os.path.join(
            output_dir, "specific_edges", f"time_breakdown_edges_size_{size}.png"
        )
    )
    plt.close()

# ----------- PLOT TIME BREAKDOWN FOR EACH SIZE BY NODES -----------
time_columns = [
    "distributionTime",
    "collectingRemoteTime",
    "detectionTime",
    "exchangingCommunitiesTime",
    "collectingMissingTime",
    "coarsenTime",
]

df_time_mean = df.groupby(["totalNodes", "size"], as_index=False)[time_columns].mean()

for size in sorted(df_time_mean["size"].unique()):
    plt.figure(figsize=(10, 7))
    df_size = df_time_mean[df_time_mean["size"] == size].copy()
    df_size = df_size.sort_values("totalNodes")
    x = df_size["totalNodes"].values
    for col in time_columns:
        if df_size[col].isnull().all():
            continue
        # Convert seconds to milliseconds for plotting
        plt.plot(x, df_size[col] * 1000, marker="o", linestyle="-", label=col)
    plt.xlabel("Total Nodes")
    plt.ylabel("Time (ms)")
    plt.xscale("log", base=10)
    plt.yscale("log", base=10)
    plt.grid(True, which="both", ls="--")
    plt.legend(fontsize="small", loc="best")
    plt.tight_layout()
    plt.subplots_adjust(left=0.1)
    plt.gca().yaxis.set_major_formatter(ScalarFormatter())
    plt.gca().ticklabel_format(style="plain", axis="y")
    plt.savefig(
        os.path.join(
            output_dir, "specific_nodes", f"time_breakdown_nodes_size_{size}.png"
        )
    )
    plt.close()
