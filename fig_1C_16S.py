"""
Figure 1C - Shannon diversity for 16S

INPUT:
    Original YAMAS OTU tables:
        datasets_after_yamas/16S/PRJNA669650/otu_PRJNA669650.csv
        datasets_after_yamas/16S/PRJNA1254708/otu_PRJNA1254708.csv

    Extra preprocessed Omri table:
        datasets_after_MIPMLP/omri/omri_stool_Pregnant_16S.csv

OUTPUT:
    1C_plots/1C_16S.png

Notes:
    - Shannon is calculated from non-normalized abundance/count tables.
    - For the original datasets: rows = OTUs, columns = samples, so we transpose.
    - For Omri: table is already after preprocess, rows = samples, columns = taxa.
    - No phylum grouping is used for Shannon.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# =====================
# Paths
# =====================
BASE_DIR = Path("/home/aharonox/Yoram_Omri_Preg")

DATA_DIR = BASE_DIR / "datasets_after_yamas/16S"
OUT_DIR = BASE_DIR / "1C_plots"
OUT_DIR.mkdir(exist_ok=True)

OMRI_PATH = BASE_DIR / "datasets_after_MIPMLP/omri/omri_stool_Pregnant_16S.csv"


# =====================
# Datasets
# =====================
DATASETS = {
    "Pregnant": ("PRJNA1254708", "Pregnant"),
    "Control": ("PRJNA669650", "Control"),
}

EXTRA_TABLES = {
    "Omri_Pregnant": (OMRI_PATH, "Pregnant"),
}


# =====================
# Load functions
# =====================
def load_otu_table(folder):
    """
    Load original YAMAS OTU table.

    Input format:
        rows = OTUs
        columns = samples

    Output format:
        rows = samples
        columns = OTUs
    """
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


def load_preprocessed_table(path):
    """
    Load Omri preprocessed table.

    Input format:
        rows = samples
        columns = taxa

    Output format:
        rows = samples
        columns = taxa
    """
    df = pd.read_csv(path, index_col=0)

    # clean names, just in case
    df.columns = df.columns.astype(str).str.strip()
    df.columns = df.columns.str.replace(r"\s*;\s*", ";", regex=True)

    # remove non-informative columns if they exist
    df = df.drop(columns=["Bacteria", "Unassigned"], errors="ignore")

    # safety: remove possible taxonomy row if it exists
    df = df[~df.index.astype(str).str.lower().str.contains("tax")]

    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)

    return df


def calculate_shannon(df):
    """
    Calculate Shannon diversity per sample.
    """
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


# =====================
# Calculate Shannon
# =====================
all_data = []

# Original YAMAS datasets
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


# Extra Omri table
for name, (path, group) in EXTRA_TABLES.items():
    if not path.exists():
        print(f"WARNING: file not found: {path}")
        continue

    df = load_preprocessed_table(path)

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
fig, ax = plt.subplots(figsize=(9, 6))

dataset_names = all_data["Dataset"].drop_duplicates().tolist()
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
    group = all_data.loc[all_data["Dataset"] == name, "Group"].iloc[0]
    patch.set_facecolor(colors[group])
    patch.set_alpha(0.5)

np.random.seed(42)

for i, name in enumerate(dataset_names):
    group = all_data.loc[all_data["Dataset"] == name, "Group"].iloc[0]
    subset = all_data[all_data["Dataset"] == name]

    ax.scatter(
        np.random.normal(i, 0.04, size=len(subset)),
        subset["Shannon"],
        color=point_colors[group],
        alpha=0.7,
        s=18
    )

ax.set_xticks(positions)
ax.set_xticklabels(dataset_names, rotation=20, ha="right")
ax.set_ylabel("Shannon Index", fontsize=13)
ax.set_title("Shannon diversity - 16S")
ax.grid(True, alpha=0.3)

plt.tight_layout()

out_path = OUT_DIR / "1C_16S.png"
plt.savefig(out_path, dpi=300, bbox_inches="tight")
plt.show()

print(f"\nSaved: {out_path}")