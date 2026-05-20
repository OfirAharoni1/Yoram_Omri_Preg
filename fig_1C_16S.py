import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


def load_clean_16s(path):
    df = pd.read_csv(path, index_col=0, sep=None, engine="python")

    # remove taxonomy row smartly
    tax_mask = df.index.astype(str).str.lower().str.contains("tax")
    if tax_mask.any():
        df = df.loc[~tax_mask]
    else:
        df = df.iloc[:-1]

    # paired-end: keep only _1, drop _2
    idx = df.index.astype(str)
    paired_mask = idx.str.endswith("_1") | idx.str.endswith("_2")

    if paired_mask.any():
        df = df.loc[~idx.str.endswith("_2")]
        df.index = df.index.astype(str).str.replace(r"_1$", "", regex=True)

    return df.apply(pd.to_numeric, errors="coerce").fillna(0)


def calculate_shannon(df):
    counts = df.to_numpy(dtype=float)
    shannon = []

    for row in counts:
        row = row[row > 0]
        if len(row) == 0:
            shannon.append(0)
        else:
            p = row / row.sum()
            shannon.append(-np.sum(p * np.log(p)))

    return np.array(shannon)


datasets = {
    "Preg_1": ("datasets_after_yamas/16S/PRJNA1254708_Pregnant_16S_for_MIPMLP.csv", "Pregnant"),
    "Ctrl_1": ("datasets_after_yamas/16S/PRJNA669650_Control_16S_for_MIPMLP.csv", "Control"),
}


all_data = []

for name, (path, group) in datasets.items():
    df_numeric = load_clean_16s(path)

    print(f"\n{name}")
    print(df_numeric.sum(axis=1).describe())

    all_data.append(pd.DataFrame({
        "Shannon": calculate_shannon(df_numeric),
        "Dataset": name,
        "Group": group
    }))


all_data = pd.concat(all_data, ignore_index=True)

fig, ax = plt.subplots(figsize=(8, 6))

dataset_names = list(datasets.keys())
positions = np.arange(len(dataset_names))

colors = {"Pregnant": "lightcoral", "Control": "lightblue"}
point_colors = {"Pregnant": "firebrick", "Control": "navy"}

box_data = [
    all_data.loc[all_data["Dataset"] == name, "Shannon"]
    for name in dataset_names
]

bp = ax.boxplot(box_data, positions=positions, patch_artist=True)

for patch, name in zip(bp["boxes"], dataset_names):
    group = datasets[name][1]
    patch.set_facecolor(colors[group])
    patch.set_alpha(0.5)

np.random.seed(42)

for i, name in enumerate(dataset_names):
    subset = all_data[all_data["Dataset"] == name]

    ax.scatter(
        np.random.normal(i, 0.04, size=len(subset)),
        subset["Shannon"],
        color=point_colors[datasets[name][1]],
        alpha=0.7,
        s=18
    )

ax.set_xticks(positions)
ax.set_xticklabels(dataset_names, rotation=45)
ax.set_ylabel("Shannon Index", fontsize=13)
ax.grid(True, alpha=0.3)

plt.tight_layout()

os.makedirs("Yoram_Omri_Preg/1C_plots", exist_ok=True)
plt.savefig("Yoram_Omri_Preg/1C_plots/1C_16S.png", dpi=300)

plt.show()