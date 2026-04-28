import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


# ======================
# calculate Shannon
# ======================
def calculate_shannon(df):
    shannon_vals = []

    for i in range(df.shape[0]):
        counts = df.iloc[i].values.astype(float)
        counts = counts[counts > 0]

        if len(counts) == 0:
            shannon_vals.append(0)
            continue

        p = counts / counts.sum()
        shannon = -np.sum(p * np.log(p))
        shannon_vals.append(shannon)

    return np.array(shannon_vals)


# ======================
# insert datasets
# ======================
datasets = {
    "Preg_1": ("preg1.csv", "Pregnant"),
    "Preg_2": ("preg2.csv", "Pregnant"),
    "Ctrl_1": ("control1.csv", "Control"),
    "Ctrl_2": ("control2.csv", "Control"),
}


all_data = []

for name, (path, group) in datasets.items():
    df = pd.read_csv(path, index_col=0)

    shannon_vals = calculate_shannon(df)

    temp_df = pd.DataFrame({
        "Shannon": shannon_vals,
        "Dataset": name,
        "Group": group
    })

    all_data.append(temp_df)

all_data = pd.concat(all_data, ignore_index=True)


# ======================
# ציור הפיגר
# ======================
fig, ax = plt.subplots(figsize=(10, 6))

dataset_names = list(datasets.keys())
positions = np.arange(len(dataset_names))

colors = {
    "Pregnant": "lightcoral",
    "Control": "lightblue"
}

point_colors = {
    "Pregnant": "firebrick",
    "Control": "navy"
}

box_data = []

for name in dataset_names:
    vals = all_data[all_data["Dataset"] == name]["Shannon"]
    box_data.append(vals)

bp = ax.boxplot(box_data, positions=positions, patch_artist=True)

# צבעים לפי קבוצה
for patch, name in zip(bp["boxes"], dataset_names):
    group = datasets[name][1]
    patch.set_facecolor(colors[group])
    patch.set_alpha(0.5)

# נקודות
np.random.seed(42)
for i, name in enumerate(dataset_names):
    subset = all_data[all_data["Dataset"] == name]

    x = np.random.normal(i, 0.04, size=len(subset))
    y = subset["Shannon"]

    ax.scatter(
        x, y,
        color=point_colors[datasets[name][1]],
        alpha=0.7,
        s=18
    )

# עיצוב
ax.set_xticks(positions)
ax.set_xticklabels(dataset_names, rotation=45)
ax.set_ylabel("Shannon Index", fontsize=13)

ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("shannon_by_dataset.png", dpi=300)
plt.show()