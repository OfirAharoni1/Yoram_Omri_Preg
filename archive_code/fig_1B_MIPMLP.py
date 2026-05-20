import os
import glob
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

warnings.filterwarnings("ignore")

# ============================================================
# SETTINGS
# ============================================================
input_folder = "/home/aharonox/Yoram_Omri_Preg/datasets_after_MIPMLP"
output_folder = "/home/aharonox/Yoram_Omri_Preg/1B_plots/16S"
os.makedirs(output_folder, exist_ok=True)

top_n = 12


# ============================================================
# TAXONOMY CLEANING
# ============================================================
def to_phylum(col):
    col = str(col).strip()

    if col.lower() in ["unassigned", "nan", "none", ""]:
        return "Unassigned"

    # Handle comma format:
    # k:Bacteria,p:Firmicutes,c:Bacilli
    if "," in col and ":" in col:
        kingdom = None
        phylum = None

        for part in col.split(","):
            part = part.strip()

            if part.startswith("k:"):
                kingdom = part.split(":", 1)[1]

            elif part.startswith("p:"):
                phylum = part.split(":", 1)[1]

        if phylum:
            return f"{kingdom};{phylum}" if kingdom else phylum

        if kingdom:
            return kingdom

    # Handle semicolon format:
    # k:Bacteria;p:Firmicutes;c:Bacilli
    if ";" in col and ":" in col:
        kingdom = None
        phylum = None

        for part in col.split(";"):
            part = part.strip()

            if part.startswith("k:"):
                kingdom = part.split(":", 1)[1]

            elif part.startswith("p:"):
                phylum = part.split(":", 1)[1]

        if phylum:
            return f"{kingdom};{phylum}" if kingdom else phylum

        if kingdom:
            return kingdom

    # Already collapsed:
    # Bacteria;Firmicutes
    if ";" in col:
        parts = [x.strip() for x in col.split(";") if x.strip()]
        if len(parts) >= 2:
            return f"{parts[0]};{parts[1]}"
        if len(parts) == 1:
            return parts[0]

    return col
    
# ============================================================
# PREPARE MIPMLP TABLE
# ============================================================
def prepare_mipmlp_table(df):
    """
    Input:
    - First column = SampleID
    - Other columns = taxa after MIPMLP

    Output:
    - Relative abundance table (%) collapsed to phylum level
    """

    id_col = df.columns[0]
    df = df.rename(columns={id_col: "SampleID"})
    df["SampleID"] = df["SampleID"].astype(str)

    # Remove paired-end read 2 samples
    df = df[~df["SampleID"].str.endswith("_2")]

    counts = df.set_index("SampleID")
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0)

    # Remove empty samples
    counts = counts.loc[counts.sum(axis=1) > 0]

    if counts.empty:
        return None

    # Collapse taxa columns to phylum level
    counts.columns = [to_phylum(c) for c in counts.columns]
    counts = counts.T.groupby(level=0).sum().T

    # Convert to relative abundance (%)
    rel = counts.div(counts.sum(axis=1), axis=0) * 100

    print("Taxa after collapse:", rel.columns.tolist())

    return rel


# ============================================================
# LOAD FILES
# ============================================================
csv_files = sorted(glob.glob(os.path.join(input_folder, "*.csv")))
print(f"Found {len(csv_files)} CSV files")


# ============================================================
# BUILD GLOBAL COLOR MAP
# ============================================================
all_top_taxa = []

for file_path in csv_files:
    file_name = os.path.splitext(os.path.basename(file_path))[0]
    print(f"Scanning for colors: {file_name}")

    df = pd.read_csv(file_path)
    rel = prepare_mipmlp_table(df)

    if rel is None:
        print(f"Skipping {file_name} in color-map stage")
        continue

    taxa_means = rel.mean(axis=0).sort_values(ascending=False)
    top_taxa = taxa_means.head(top_n).index.tolist()
    all_top_taxa.extend(top_taxa)

global_taxa_order = list(dict.fromkeys(all_top_taxa))

colors = (
    list(plt.cm.tab20.colors) +
    list(plt.cm.Set3.colors) +
    list(plt.cm.Pastel1.colors)
)

global_color_map = {
    taxon: colors[i % len(colors)]
    for i, taxon in enumerate(global_taxa_order)
}

global_color_map["Other"] = "#BDBDBD"


# ============================================================
# MAIN LOOP
# ============================================================
all_records = []

for file_path in csv_files:
    file_name = os.path.splitext(os.path.basename(file_path))[0]
    title_name = file_name.replace("_MIPMLP", "").split("_")[0]

    print(f"\nProcessing: {file_name}")

    try:
        df = pd.read_csv(file_path)
        rel = prepare_mipmlp_table(df)

        if rel is None:
            print(f"Skipping {file_name}: could not prepare table")
            continue

        taxa_means = rel.mean(axis=0).sort_values(ascending=False)
        top_taxa = taxa_means.head(top_n).index.tolist()

        if "Control" in file_name:
            group_name = "Control"
        elif "Pregnant" in file_name:
            group_name = "Pregnant"
        else:
            group_name = "Unknown"

        # Collect records for combined violin plot
        for taxon in top_taxa:
            all_records.append(pd.DataFrame({
                "dataset": title_name,
                "group": group_name,
                "taxon": taxon,
                "value": rel[taxon].values
            }))

        # Prepare stacked bar data
        rel_top = rel[top_taxa].copy()

        other_taxa = [c for c in rel.columns if c not in top_taxa]
        if other_taxa:
            rel_top["Other"] = rel[other_taxa].sum(axis=1)

        rel_top = rel_top[rel_top.mean(axis=0).sort_values(ascending=False).index]

        if "Bacteria;Firmicutes" in rel_top.columns:
            rel_top = rel_top.sort_values(by="Bacteria;Firmicutes", ascending=False)

        final_order = rel_top.mean(axis=0).sort_values(ascending=False).index.tolist()
        rel_top = rel_top[final_order]

        # Make sure every taxon has a color
        for taxon in rel_top.columns:
            if taxon not in global_color_map:
                global_color_map[taxon] = colors[len(global_color_map) % len(colors)]

        color_map = {taxon: global_color_map[taxon] for taxon in rel_top.columns}

        # ========================================================
        # STACKED BAR PLOT
        # ========================================================
        fig, axes = plt.subplots(
            1, 2,
            figsize=(16, 6),
            gridspec_kw={"width_ratios": [5, 1]}
        )

        ax_samples, ax_mean = axes

        # Left panel: samples
        x_positions = range(len(rel_top))
        bottom = np.zeros(len(rel_top))

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

            bottom += values

        ax_samples.set_ylabel("Relative abundance (%)", fontsize=12)
        ax_samples.set_ylim(0, 100)
        ax_samples.set_xticks([])
        ax_samples.set_title("Samples", fontsize=13)
        ax_samples.spines["top"].set_visible(False)
        ax_samples.spines["right"].set_visible(False)

        # Right panel: average composition
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

        out_path = os.path.join(output_folder, f"{file_name}_composition.png")
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close()

        print(f"Saved stacked plot to: {out_path}")

    except Exception as e:
        print(f"Skipping {file_name}: {repr(e)}")


# ============================================================
# COMBINED SPLIT VIOLIN PLOT
# ============================================================
if len(all_records) == 0:
    print("No valid records collected for combined violin plot.")

else:
    violin_df = pd.concat(all_records, ignore_index=True)

    violin_df["group"] = violin_df["group"].astype(str).str.strip()
    violin_df = violin_df[violin_df["group"].isin(["Control", "Pregnant"])]

    if violin_df.empty:
        print("No Control/Pregnant rows found for violin plot.")

    else:
        taxon_order = (
            violin_df.groupby("taxon")["value"]
            .mean()
            .sort_values(ascending=False)
            .index
            .tolist()
        )

        fig, ax = plt.subplots(figsize=(max(14, len(taxon_order) * 0.8), 7))

        control_color = "#4C72B0"
        pregnant_color = "#DD8452"

        for i, taxon in enumerate(taxon_order):
            control_vals = violin_df[
                (violin_df["taxon"] == taxon) &
                (violin_df["group"] == "Control")
            ]["value"].dropna().values

            pregnant_vals = violin_df[
                (violin_df["taxon"] == taxon) &
                (violin_df["group"] == "Pregnant")
            ]["value"].dropna().values

            center = i

            # Control = left half
            if len(control_vals) > 1:
                parts = ax.violinplot(
                    [control_vals],
                    positions=[center],
                    widths=0.8,
                    showmeans=False,
                    showmedians=True,
                    showextrema=False
                )

                for body in parts["bodies"]:
                    body.set_facecolor(control_color)
                    body.set_edgecolor("black")
                    body.set_alpha(0.7)

                    verts = body.get_paths()[0].vertices
                    verts[:, 0] = np.minimum(verts[:, 0], center)

                if "cmedians" in parts:
                    parts["cmedians"].set_color("black")
                    parts["cmedians"].set_linewidth(1.2)

            # Pregnant = right half
            if len(pregnant_vals) > 1:
                parts = ax.violinplot(
                    [pregnant_vals],
                    positions=[center],
                    widths=0.8,
                    showmeans=False,
                    showmedians=True,
                    showextrema=False
                )

                for body in parts["bodies"]:
                    body.set_facecolor(pregnant_color)
                    body.set_edgecolor("black")
                    body.set_alpha(0.7)

                    verts = body.get_paths()[0].vertices
                    verts[:, 0] = np.maximum(verts[:, 0], center)

                if "cmedians" in parts:
                    parts["cmedians"].set_color("black")
                    parts["cmedians"].set_linewidth(1.2)

        ax.set_xticks(range(len(taxon_order)))
        ax.set_xticklabels(taxon_order, rotation=45, ha="right")
        ax.set_ylabel("Relative abundance (%)")
        ax.set_title("Top phyla across 16S datasets: Control vs Pregnant")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        ax.legend(handles=[
            Patch(facecolor=control_color, edgecolor="black", label="Control"),
            Patch(facecolor=pregnant_color, edgecolor="black", label="Pregnant")
        ], frameon=False)

        plt.tight_layout()

        violin_out_path = os.path.join(output_folder, "violin_plot_16S_MIPMLP.png")
        plt.savefig(violin_out_path, dpi=300, bbox_inches="tight")
        plt.close()

        print(f"Saved combined split violin plot to: {violin_out_path}")


print("\nDone.")