"""
Figure 1C - Shannon diversity for 16S

INPUT:
    datasets_after_yamas/16S/
    ├── PRJNA669650/
    │   └── otu_PRJNA669650.csv
    └── PRJNA1254708/
        └── otu_PRJNA1254708.csv

OUTPUT:
    1C_plots/1C_16S.png

Notes:
    - Shannon is calculated from the original OTU table.
    - No taxonomy grouping is performed.
    - No MIPMLP preprocessed table is used.
    - If paired-end samples exist, _2 is removed and _1 suffix is cleaned.
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


DATA_DIR = Path("datasets_after_yamas/16S")
OUT_DIR = Path("1C_plots")
OUT_DIR.mkdir(exist_ok=True)

DATASETS = {
    "Pregnant": ("PRJNA1254708", "Pregnant"),
    "Control": ("PRJNA669650", "Control"),
}


def load_otu_table(folder):
    otu_path = next(folder.glob("otu_*.csv"))

    otu = pd.read_csv(otu_path, index_col=0)

    # OTU table: rows = OTUs, columns = samples
    # Shannon needs rows = samples, columns = OTUs
    df = otu.T
    df.index.name = "SampleID"

    # paired-end: keep only _1, remove _2
    df = df[~df.index.astype(str).str.endswith("_2")]
    df.index = df.index.astype(str).str.replace(r"_1$", "", regex=True)

    # safety: remove possible taxonomy row if it exists
    df = df[~df.index.astype(str).str.lower().str.contains("tax")]

    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)

    return df


def calculate_shannon(df):
    shannon = []

    for _, row in df.iterrows():
        values = row.values.astype(float)
        values = values[values > 0]

        if len(values) == 0:
            shannon.append(0)
        else:
            p = values / values.sum()
            shannon.append(-np.sum(p * np.log(p)))

    return np.array(shannon)


all_data = []

for name, (folder_name, group) in DATASETS.items():
    folder = DATA_DIR / folder_name

    df = load_otu_table(folder)

    print(f"\n{name}")
    print(f"Shape: {df.shape}")
    print(df.sum(axis=1).describe())

    all_data.append(pd.DataFrame({
        "Shannon": calculate_shannon(df),
        "Dataset": name,
        "Group": group
    }))

all_data = pd.concat(all_data, ignore_index=True)


# =====================
# Plot
# =====================
fig, ax = plt.subplots(figsize=(8, 6))

dataset_names = list(DATASETS.keys())
positions = np.arange(len(dataset_names))

colors = {
    "Pregnant": "lightcoral",
    "Control": "lightblue"
}

point_colors = {
    "Pregnant": "firebrick",
    "Control": "navy"
}

box_data = [
    all_data.loc[all_data["Dataset"] == name, "Shannon"]
    for name in dataset_names
]

bp = ax.boxplot(box_data, positions=positions, patch_artist=True)

for patch, name in zip(bp["boxes"], dataset_names):
    group = DATASETS[name][1]
    patch.set_facecolor(colors[group])
    patch.set_alpha(0.5)

np.random.seed(42)

for i, name in enumerate(dataset_names):
    group = DATASETS[name][1]
    subset = all_data[all_data["Dataset"] == name]

    ax.scatter(
        np.random.normal(i, 0.04, size=len(subset)),
        subset["Shannon"],
        color=point_colors[group],
        alpha=0.7,
        s=18
    )

ax.set_xticks(positions)
ax.set_xticklabels(dataset_names)
ax.set_ylabel("Shannon Index", fontsize=13)
ax.set_title("Shannon diversity - 16S")
ax.grid(True, alpha=0.3)

plt.tight_layout()

out_path = OUT_DIR / "1C_16S.png"
plt.savefig(out_path, dpi=300, bbox_inches="tight")
plt.show()

print(f"\nSaved: {out_path}")