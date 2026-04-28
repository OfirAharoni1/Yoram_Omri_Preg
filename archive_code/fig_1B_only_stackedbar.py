import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# ===== SETTINGS =====
input_folder = "/home/aharonox/Yoram_Omri_Preg/datasets_after_yamas"
output_folder = "/home/aharonox/Yoram_Omri_Preg/1B_plots"
os.makedirs(output_folder, exist_ok=True)

sample_id_col = "SampleID"
top_n = 12   # keep top phyla, merge the rest into "Other"

def extract_phylum(taxon_name):
    """
    Extract phylum level from taxonomy string.
    Example:
    k__Bacteria,p__Firmicutes,c__Clostridia -> Bacteria;Firmicutes
    """
    parts = str(taxon_name).split(",")
    kingdom = None
    phylum = None

    for part in parts:
        part = part.strip()
        if part.startswith("k__"):
            kingdom = part.replace("k__", "")
        elif part.startswith("p__"):
            phylum = part.replace("p__", "")

    if phylum:
        if kingdom:
            return f"{kingdom};{phylum}"
        return phylum

    # If no phylum exists, return None so we can skip it
    return None


csv_files = glob.glob(os.path.join(input_folder, "*.csv"))
print(f"Found {len(csv_files)} CSV files")

for file_path in csv_files:
    file_name = os.path.splitext(os.path.basename(file_path))[0]
    title_name = file_name.split("_")[0]

    print(f"Processing: {file_name}")

    # ===== LOAD =====
    df = pd.read_csv(file_path)

    # Rename first column if needed
    if df.columns[0] != sample_id_col:
        df = df.rename(columns={df.columns[0]: sample_id_col})

    df = df.set_index(sample_id_col)

    # Keep numeric columns only
    abundance = df.select_dtypes(include="number").copy()

    if abundance.shape[1] == 0:
        print(f"Skipping {file_name}: no numeric abundance columns found")
        continue

    # Remove zero-sum samples
    row_sums = abundance.sum(axis=1)
    abundance = abundance.loc[row_sums > 0]

    if abundance.empty:
        print(f"Skipping {file_name}: all samples have zero total abundance")
        continue

    # ===== COLLAPSE TO PHYLUM LEVEL =====
    phylum_map = {}
    for col in abundance.columns:
        phylum = extract_phylum(col)
        if phylum is not None:
            phylum_map[col] = phylum

    if len(phylum_map) == 0:
        print(f"Skipping {file_name}: could not extract phylum names")
        continue

    abundance_phylum = abundance[list(phylum_map.keys())].copy()
    abundance_phylum.columns = [phylum_map[c] for c in abundance_phylum.columns]

    # Merge duplicate phylum columns by summing
    abundance_phylum = abundance_phylum.T.groupby(level=0).sum().T
    # Convert to relative abundance (%)
    rel = abundance_phylum.div(abundance_phylum.sum(axis=1), axis=0) * 100

    # ===== KEEP TOP PHYLA =====
    phylum_means = rel.mean(axis=0).sort_values(ascending=False)
    top_phyla = phylum_means.head(top_n).index.tolist()

    rel_top = rel[top_phyla].copy()

    other_phyla = [c for c in rel.columns if c not in top_phyla]
    if other_phyla:
        rel_top["Other"] = rel[other_phyla].sum(axis=1)

    # Sort taxa by mean abundance (important for clean visualization)
    rel_top = rel_top[rel_top.mean(axis=0).sort_values(ascending=False).index]

    # Sort samples by the dominant phylum
    if "Bacteria;Firmicutes" in rel_top.columns:
        rel_top = rel_top.sort_values(by="Bacteria;Firmicutes", ascending=False)

    # Order by mean abundance
    final_order = rel_top.mean(axis=0).sort_values(ascending=False).index.tolist()
    rel_top = rel_top[final_order]

    # ===== COLORS =====
    colors = (
        list(plt.cm.tab20.colors) +
        list(plt.cm.Set3.colors) +
        list(plt.cm.Pastel1.colors)
    )
    color_map = {taxon: colors[i % len(colors)] for i, taxon in enumerate(rel_top.columns)}

    # ===== PLOT =====
    fig, axes = plt.subplots(
        1, 2,
        figsize=(16, 6),
        gridspec_kw={"width_ratios": [5, 1]}
    )

    ax_samples, ax_mean = axes

    # Left panel: samples
    x_positions = range(len(rel_top))
    bottom = [0] * len(rel_top)

    for taxon in rel_top.columns:
        values = rel_top[taxon].values
        ax_samples.bar(
            x_positions,
            values,
            bottom=bottom,
            width=0.9,
            color=color_map[taxon],
            edgecolor="none"
        )
        bottom = [b + v for b, v in zip(bottom, values)]

    ax_samples.set_ylabel("Relative abundance (%)", fontsize=12)
    ax_samples.set_ylim(0, 100)
    ax_samples.set_xticks([])
    ax_samples.set_title("Samples", fontsize=13)
    ax_samples.spines["top"].set_visible(False)
    ax_samples.spines["right"].set_visible(False)

    # Right panel: mean
    mean_values = rel_top.mean(axis=0)
    bottom = 0

    for taxon in rel_top.columns:
        value = mean_values[taxon]
        ax_mean.bar(
            [0],
            [value],
            bottom=bottom,
            width=0.9,
            color=color_map[taxon],
            edgecolor="none"
        )
        bottom += value

    ax_mean.set_ylim(0, 100)
    ax_mean.set_xticks([0])
    ax_mean.set_xticklabels(["Mean"], fontsize=11)
    ax_mean.set_title("Average", fontsize=13)
    ax_mean.spines["top"].set_visible(False)
    ax_mean.spines["right"].set_visible(False)
    ax_mean.spines["left"].set_visible(False)
    ax_mean.yaxis.set_visible(False)

    # Legend
    legend_handles = [
        Rectangle((0, 0), 1, 1, color=color_map[taxon], label=taxon)
        for taxon in rel_top.columns
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=4,
        frameon=False,
        fontsize=9
    )

    fig.suptitle(title_name, fontsize=16, fontweight="bold")

    plt.tight_layout(rect=[0, 0.12, 1, 1])

    out_path = os.path.join(output_folder, f"{file_name}_phylum_composition.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved to: {out_path}")

print("Done.")