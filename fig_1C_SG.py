import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


# ======================
# calculate Shannon
# ======================
def calculate_shannon(df):
    shannon_vals = []

    # נשאיר רק עמודות מספריות
    df_numeric = df.select_dtypes(include=[np.number])

    for i in range(df_numeric.shape[0]):
        counts = df_numeric.iloc[i].values.astype(float)
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
    "Preg_1": ("Yoram_Omri_Preg/datasets_after_yamas/SG/ERP020710_Pregnant_SG_formatted.csv", "Pregnant"),
    "Preg_2": ("Yoram_Omri_Preg/datasets_after_yamas/SG/PRJNA1247940_Pregnant_SG_formatted.csv", "Pregnant"),
    "Ctrl_1": ("Yoram_Omri_Preg/datasets_after_yamas/SG/PRJEB37731_Control_SG_formatted.csv", "Control"),
    "Ctrl_2": ("Yoram_Omri_Preg/datasets_after_yamas/SG/PRJNA48479_Control_SG_formatted.csv", "Control"),
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
# create figure
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

# design
ax.set_xticks(positions)
ax.set_xticklabels(dataset_names, rotation=45)
ax.set_ylabel("Shannon Index", fontsize=13)

ax.grid(True, alpha=0.3)

plt.tight_layout()

os.makedirs("Yoram_Omri_Preg/1C_plots", exist_ok=True)

plt.savefig("Yoram_Omri_Preg/1C_plots/1C_SG.png", dpi=300)
plt.show()