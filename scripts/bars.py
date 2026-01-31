import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# Ask user for input file
input_file = input("Enter the CSV file name (e.g., out.csv): ").strip()

# Load the CSV
df = pd.read_csv(input_file)
output_dir = os.path.dirname(os.path.abspath(input_file))
dir_path = os.path.join(output_dir, "bars")
os.makedirs(dir_path, exist_ok=True)

# Ensure numeric fields are of the correct type
df["size"] = df["size"].astype(int)
df["totalTime"] = df["totalTime"].astype(float)

# Compute the mean for each fileName and size
df_mean = df.groupby(["fileName", "size"], as_index=False)[["totalTime"]].mean()

time_columns = [
    "totalTime",
    "partitionTime",
    "distributionTime",
    "collectingRemoteTime",
    "detectionTime",
    "exchangingCommunitiesTime",
    "collectingMissingTime",
    "coarsenTime",
]

# Compute the mean for each time column as well
df_time_means = df.groupby(["fileName", "size"], as_index=False)[time_columns].mean()

# Calculate percentages over totalTime for each time column except totalTime itself
for col in time_columns:
    if col != "totalTime":
        df_time_means[col] = (df_time_means[col] / df_time_means["totalTime"]) * 100

# Get unique sizes and filenames
sizes = df_time_means["size"].unique()
file_names = df_time_means["fileName"].unique()

for size in sizes:
    df_size = df_time_means[df_time_means["size"] == size]
    # Prepare data for stacked bar (exclude totalTime, which is always 100%)
    bar_data = [df_size[col].values for col in time_columns if col != "totalTime"]
    ind = np.arange(len(df_size))
    bottom = np.zeros(len(df_size))
    plt.figure(figsize=(10, 6))
    labels = [col for col in time_columns if col != "totalTime"]
    # Calculate sum of known percentages for each bar
    known_sum = np.sum(bar_data, axis=0)
    # Calculate "Other Computation" as the remaining percentage to 100%
    other_computation = 100 - known_sum
    # Add "Other Computation" to bar_data and labels
    bar_data.append(other_computation)
    labels.append("other computation")
    for i, col in enumerate(labels):
        plt.bar(ind, bar_data[i], bottom=bottom, label=col)
        bottom += bar_data[i]
    plt.xticks(ind, [fn.replace('.txt', '') for fn in df_size["fileName"]], rotation=45, ha="right")
    plt.ylabel("Time (%)")
    # plt.title(f"Stacked Bar Chart of Times (Percentage) for size={size}")
    plt.legend(loc="center left")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "bars", f"stacked_bar_size_{size}.png"))
