import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Carica il CSV
df = pd.read_csv("out.csv")

# Assicurati che i campi numerici siano del tipo corretto
df["size"] = df["size"].astype(int)
df["totalTime"] = df["totalTime"].astype(float)
df["partitionTime"] = df["partitionTime"].astype(float)
df["detectionTime"] = df["detectionTime"].astype(float)

# Calcola la media per ogni fileName e size
df_mean = df.groupby(["fileName", "size"], as_index=False)[
    ["totalTime", "partitionTime", "detectionTime"]
].mean()

# Trova tutti i fileName unici
for file_name in df_mean["fileName"].unique():
    df_file = df_mean[df_mean["fileName"] == file_name].copy()
    df_file = df_file.sort_values("size")

    # Trova i tempi base per size == 1
    base_total = df_file[df_file["size"] == 1]["totalTime"].values
    base_partition = df_file[df_file["size"] == 1]["partitionTime"].values
    base_detection = df_file[df_file["size"] == 1]["detectionTime"].values

    if len(base_total) == 0:
        continue  # Salta se manca size 1

    base_total = base_total[0]
    base_partition = base_partition[0]
    base_detection = base_detection[0]

    # Calcola strong scaling (speedup)
    df_file["speedup_total"] = base_total / df_file["totalTime"]
    df_file["speedup_partition"] = base_partition / df_file["partitionTime"]
    df_file["speedup_detection"] = base_detection / df_file["detectionTime"]

    # Calcola efficiency in %
    df_file["eff_total"] = (df_file["speedup_total"] / df_file["size"]) * 100
    df_file["eff_partition"] = (df_file["speedup_partition"] / df_file["size"]) * 100
    df_file["eff_detection"] = (df_file["speedup_detection"] / df_file["size"]) * 100

    # Ordina per size e prendi i dati
    df_file = df_file.sort_values("size")
    x = df_file["size"].values

    # ----------- GRAFICO SPEEDUP -----------
    plt.figure()
    plt.plot(x, df_file["speedup_total"], marker="o", linestyle="-", label="Total Time")
    plt.plot(
        x,
        df_file["speedup_partition"],
        marker="s",
        linestyle="--",
        label="Partition Time",
    )
    plt.plot(
        x,
        df_file["speedup_detection"],
        marker="^",
        linestyle=":",
        label="Detection Time",
    )

    plt.xscale("log", base=2)
    plt.xlabel("Numero di processi (size)")
    plt.ylabel("Strong scaling (speedup)")
    plt.title(f"Strong Scaling per {file_name}")
    plt.grid(True, which="both", ls="--")
    plt.xticks(x)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"strong_scaling_{file_name}.png")
    plt.close()

    # ----------- GRAFICO EFFICIENCY -----------
    plt.figure()
    plt.plot(x, df_file["eff_total"], marker="o", linestyle="-", label="Total Time")
    plt.plot(
        x, df_file["eff_partition"], marker="s", linestyle="--", label="Partition Time"
    )
    plt.plot(
        x, df_file["eff_detection"], marker="^", linestyle=":", label="Detection Time"
    )

    plt.xscale("log", base=2)
    plt.xlabel("Numero di processi (size)")
    plt.ylabel("Efficiency (%)")
    plt.title(f"Efficiency per {file_name}")
    plt.grid(True, which="both", ls="--")
    plt.xticks(x)

    # Limite asse Y dinamico
    y_max = max(
        df_file["eff_total"].max(),
        df_file["eff_partition"].max(),
        df_file["eff_detection"].max(),
    )
    y_limit = np.ceil(y_max / 10) * 10
    plt.ylim(0, max(110, y_limit + 10))

    plt.legend()
    plt.tight_layout()
    plt.savefig(f"efficiency_{file_name}.png")
    plt.close()
