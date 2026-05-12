import os
import glob
import traceback
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

warnings.filterwarnings("ignore")

# ===== SETTINGS =====
input_folder = "/home/aharonox/Yoram_Omri_Preg/datasets_after_yamas/16S"
output_folder = "/home/aharonox/Yoram_Omri_Preg/1B_plots/16S"
os.makedirs(output_folder, exist_ok=True)

top_n = 12


# ===== TAXONOMY FUNCTIONS =====
def clean_taxonomy_string(x):
    if pd.isna(x):
        return None

    x = str(x).strip()

    # QIIME/SILVA format:
    # k__Bacteria;p__Firmicutes;c__Clostridia
    if ";" in x and "__" in x:
        parts = []
        for part in x.split(";"):
            part = part.strip()
            if "__" in part:
                part = part.split("__", 1)[1]
            if part:
                parts.append(part)
        return ";".join(parts)

    # SINTAX/vsearch format:
    # k:Bacteria,p:Firmicutes,c:Clostridia
    if "," in x:
        parts = []
        for part in x.split(","):
            part = part.strip()
            if ":" in part:
                part = part.split(":", 1)[1]
            if part:
                parts.append(part)
        return ";".join(parts)

    return x


def extract_phylum_16s(taxon_name):
    """
    Convert taxonomy to phylum level:
    k__Bacteria;p__Firmicutes;c__Clostridia -> Bacteria;Firmicutes
    Bacteria;Firmicutes;Clostridia -> Bacteria;Firmicutes
    """
    tax = clean_taxonomy_string(taxon_name)

    if tax is None:
        return None

    parts = [p.strip() for p in str(tax).split(";") if p.strip()]

    if len(parts) < 2:
        return None

    return f"{parts[0]};{parts[1]}"


# ===== PREPARE 16S TABLE =====
def prepare_relative_phylum_table_16s(df):
    """
    Expected input format:
    - first column = ID
    - one row where ID == taxonomy
    - other rows = samples
    - columns = OTUs/ASVs/features
    """
    id_col = df.columns[0]

    taxonomy_mask = df[id_col].astype(str).str.lower().eq("taxonomy")

    if not taxonomy_mask.any():
        return None

    tax_idx = df.index[taxonomy_mask][0]

    taxonomy_row = df.loc[tax_idx].drop(labels=id_col)
    samples = df.drop(index=tax_idx).copy()

    sample_ids = samples[id_col].astype(str)
    counts = samples.drop(columns=[id_col])

    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0)
    counts.index = sample_ids

    phylum_map = {}

    for col in counts.columns:
        phylum = extract_phylum_16s(taxonomy_row[col])
        if phylum is not None:
            phylum_map[col] = phylum

    if len(phylum_map) == 0:
        return None

    abundance_phylum = counts[list(phylum_map.keys())].copy()
    abundance_phylum.columns = [phylum_map[c] for c in abundance_phylum.columns]

    # Merge duplicate phylum columns by summing
    abundance_phylum = abundance_phylum.T.groupby(level=0).sum().T

    # Remove zero-sum samples
    row_sums = abundance_phylum.sum(axis=1)
    abundance_phylum = abundance_phylum.loc[row_sums > 0]

    if abundance_phylum.empty:
        return None

    # Convert to relative abundance (%)
    rel = abundance_phylum.div(abundance_phylum.sum(axis=1), axis=0) * 100

    return rel


# ===== FILES =====
csv_files = sorted(glob.glob(os.path.join(input_folder, "*.csv")))
print(f"Found {len(csv_files)} CSV files")


# ===== BUILD GLOBAL COLOR MAP =====
all_top_taxa = []

for file_path in csv_files:
    file_name = os.path.splitext(os.path.basename(file_path))[0]
    print(f"Scanning for colors: {file_name}")

    try:
        df = pd.read_csv(file_path)
        rel = prepare_relative_phylum_table_16s(df)

        if rel is None:
            print(f"Skipping {file_name} in color-map stage")
            continue

        phylum_means = rel.mean(axis=0).sort_values(ascending=False)
        top_phyla = phylum_means.head(top_n).index.tolist()
        all_top_taxa.extend(top_phyla)

    except Exception as e:
        print(f"Skipping {file_name} in color-map stage: {repr(e)}")


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


# ===== MAIN LOOP =====
all_records = []

for file_path in csv_files:
    file_name = os.path.splitext(os.path.basename(file_path))[0]
    title_name = file_name.split("_")[0]

    print(f"\nProcessing: {file_name}")

    try:
        df = pd.read_csv(file_path)
        rel = prepare_relative_phylum_table_16s(df)

        if rel is None:
            print(f"Skipping {file_name}: could not prepare relative phylum table")
            continue

        # ===== KEEP TOP PHYLA FOR THIS DATASET =====
        phylum_means = rel.mean(axis=0).sort_values(ascending=False)
        top_phyla = phylum_means.head(top_n).index.tolist()

        # ===== COLLECT RECORDS FOR COMBINED VIOLIN =====
        if "Control" in file_name:
            group_name = "Control"
        elif "Pregnant" in file_name:
            group_name = "Pregnant"
        else:
            group_name = "Unknown"

        for taxon in top_phyla:
            tmp = pd.DataFrame({
                "dataset": title_name,
                "group": group_name,
                "taxon": taxon,
                "value": rel[taxon].values
            })
            all_records.append(tmp)

        # ===== STACKED PLOT FOR THIS DATASET =====
        rel_top = rel[top_phyla].copy()

        other_phyla = [c for c in rel.columns if c not in top_phyla]

        if other_phyla:
            rel_top["Other"] = rel[other_phyla].sum(axis=1)

        # Sort taxa by mean abundance
        rel_top = rel_top[rel_top.mean(axis=0).sort_values(ascending=False).index]

        # Sort samples same as SG code
        if "Bacteria;Firmicutes" in rel_top.columns:
            rel_top = rel_top.sort_values(by="Bacteria;Firmicutes", ascending=False)

        # Final taxon order
        final_order = rel_top.mean(axis=0).sort_values(ascending=False).index.tolist()
        rel_top = rel_top[final_order]

        # Colors
        for taxon in rel_top.columns:
            if taxon not in global_color_map:
                global_color_map[taxon] = colors[len(global_color_map) % len(colors)]

        color_map = {taxon: global_color_map[taxon] for taxon in rel_top.columns}

        # ===== STACKED BAR PLOT =====
        fig, axes = plt.subplots(
            1, 2,
            figsize=(16, 6),
            gridspec_kw={"width_ratios": [5, 1]}
        )

        ax_samples, ax_mean = axes

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

        # Mean panel
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

        print(f"Saved stacked plot to: {out_path}")

    except Exception as e:
        print(f"Skipping {file_name}: {repr(e)}")
        traceback.print_exc()


# ===== COMBINED SPLIT VIOLIN PLOT =====
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
                    med = parts["cmedians"]
                    med.set_color("black")
                    med.set_linewidth(1.2)

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
                    med = parts["cmedians"]
                    med.set_color("black")
                    med.set_linewidth(1.2)

        ax.set_xticks(range(len(taxon_order)))
        ax.set_xticklabels(taxon_order, rotation=45, ha="right")
        ax.set_ylabel("Relative abundance (%)")
        ax.set_title("Top phyla across 16S datasets: Control vs Pregnant")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        legend_handles = [
            Patch(facecolor=control_color, edgecolor="black", label="Control"),
            Patch(facecolor=pregnant_color, edgecolor="black", label="Pregnant")
        ]

        ax.legend(handles=legend_handles, frameon=False)

        plt.tight_layout()

        violin_out_path = os.path.join(output_folder, "violin_plot_16S.png")
        plt.savefig(violin_out_path, dpi=300, bbox_inches="tight")
        plt.close()

        print(f"Saved combined split violin plot to: {violin_out_path}")


print("\nDone.")