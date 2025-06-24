import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# Ask user for input file
input_file = input("Enter the CSV file name (e.g., out.csv): ").strip()

# Load the CSV
df = pd.read_csv(input_file)

# Ensure numeric fields are of the correct type
df["size"] = df["size"].astype(int)
df["totalTime"] = df["totalTime"].astype(float)

# Compute the mean for each fileName and size
df_mean = df.groupby(["fileName", "size"], as_index=False)[["totalTime"]].mean()

# Prepare plots
plt.figure(figsize=(8, 6))
for file_name in df_mean["fileName"].unique():
    df_file = df_mean[df_mean["fileName"] == file_name].copy()
    df_file = df_file.sort_values("size")

    # Find base time for size == 1
    base_total = df_file[df_file["size"] == 1]["totalTime"].values
    if len(base_total) == 0:
        continue  # Skip if size 1 is missing
    base_total = base_total[0]

    # Compute strong scaling (speedup)
    df_file["speedup_total"] = base_total / df_file["totalTime"]

    x = df_file["size"].values
    label = file_name.rsplit(".", 1)[0]
    plt.plot(x, df_file["speedup_total"], marker="o", linestyle="-", label=label)

# Add ideal speedup line
xticks = sorted(df_mean["size"].unique())
plt.plot(xticks, xticks, color="k", linestyle="--", label="Ideal speedup")

plt.xscale("log", base=2)
plt.xlabel("Number of processes (size)")
plt.ylabel("Strong scaling (speedup)")
plt.grid(True, which="both", ls="--")
plt.xticks(xticks, labels=[str(x) for x in xticks])
plt.legend()
plt.tight_layout()
output_dir = os.path.dirname(os.path.abspath(input_file))
plt.savefig(os.path.join(output_dir, "strong_scaling_all.png"))
plt.close()

# ----------- EFFICIENCY PLOT -----------
plt.figure(figsize=(8, 6))
y_max = 0
for file_name in df_mean["fileName"].unique():
    df_file = df_mean[df_mean["fileName"] == file_name].copy()
    df_file = df_file.sort_values("size")

    base_total = df_file[df_file["size"] == 1]["totalTime"].values
    if len(base_total) == 0:
        continue
    base_total = base_total[0]

    df_file["speedup_total"] = base_total / df_file["totalTime"]
    df_file["eff_total"] = (df_file["speedup_total"] / df_file["size"]) * 100

    x = df_file["size"].values
    label = file_name.rsplit(".", 1)[0]
    plt.plot(x, df_file["eff_total"], marker="o", linestyle="-", label=label)
    y_max = max(y_max, df_file["eff_total"].max())

plt.xscale("log", base=2)
plt.xlabel("Number of processes (size)")
plt.ylabel("Efficiency (%)")
plt.grid(True, which="both", ls="--")
plt.xticks(xticks, labels=[str(x) for x in xticks])
y_limit = np.ceil(y_max / 10) * 10
plt.ylim(0, max(110, y_limit + 10))
plt.legend()
plt.tight_layout()
output_dir = os.path.dirname(os.path.abspath(input_file))
plt.savefig(os.path.join(output_dir, "efficiency_all.png"))
plt.close()
