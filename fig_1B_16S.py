import os, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

# ===== SETTINGS =====
input_folder = "/home/aharonox/Yoram_Omri_Preg/datasets_after_yamas/16S"
output_folder = "/home/aharonox/Yoram_Omri_Preg/1B_plots/16S"
os.makedirs(output_folder, exist_ok=True)

top_n = 12


# ===== TAXONOMY =====
def extract_phylum(x):
    if pd.isna(x):
        return None

    x = str(x).replace(",", ";").replace('"', "").strip()
    parts = []

    for p in x.split(";"):
        p = p.strip()
        if "__" in p:
            p = p.split("__")[-1]
        elif ":" in p:
            p = p.split(":")[-1]
        if p:
            parts.append(p)

    return f"{parts[0]};{parts[1]}" if len(parts) >= 2 else None


# ===== PREPARE 16S TABLE =====
def prepare_16s(df):
    id_col = df.columns[0]

    tax_mask = df[id_col].astype(str).str.lower().eq("taxonomy")
    if not tax_mask.any():
        return None

    tax_row = df.loc[df.index[tax_mask][0]].drop(id_col)
    samples = df.loc[~tax_mask].copy()
    samples[id_col] = samples[id_col].astype(str)

    phylum_map = {c: extract_phylum(tax_row[c]) for c in df.columns if c != id_col}
    phylum_map = {c: p for c, p in phylum_map.items() if p is not None}

    if not phylum_map:
        return None

    # If paired (_1/_2), choose the side with lower Proteobacteria
    has_1 = samples[id_col].str.endswith("_1").any()
    has_2 = samples[id_col].str.endswith("_2").any()

    if has_1 and has_2:
        scores = {}

        for suffix in ["_1", "_2"]:
            sub = samples[samples[id_col].str.endswith(suffix)].copy()
            counts = sub[list(phylum_map.keys())].apply(pd.to_numeric, errors="coerce").fillna(0)
            counts.columns = [phylum_map[c] for c in counts.columns]
            counts = counts.T.groupby(level=0).sum().T
            rel = counts.div(counts.sum(axis=1), axis=0) * 100
            scores[suffix] = rel.mean().get("Bacteria;Proteobacteria", 100)

        keep_suffix = min(scores, key=scores.get)
        print(f"Paired detected -> keeping {keep_suffix} | Proteobacteria scores: {scores}")

        samples = samples[samples[id_col].str.endswith(keep_suffix)].copy()
        samples[id_col] = samples[id_col].str.replace(r"_[12]$", "", regex=True)

    # Build abundance table
    counts = samples[[id_col] + list(phylum_map.keys())].copy()
    counts = counts.set_index(id_col)
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0)

    counts.columns = [phylum_map[c] for c in counts.columns]
    counts = counts.T.groupby(level=0).sum().T
    counts = counts.loc[counts.sum(axis=1) > 0]

    if counts.empty:
        return None

    return counts.div(counts.sum(axis=1), axis=0) * 100


# ===== FILES =====
csv_files = sorted(glob.glob(os.path.join(input_folder, "*.csv")))
print(f"Found {len(csv_files)} CSV files")


# ===== GLOBAL COLOR MAP =====
all_top_taxa = []

for f in csv_files:
    rel = prepare_16s(pd.read_csv(f, low_memory=False))
    if rel is not None:
        all_top_taxa += rel.mean().sort_values(ascending=False).head(top_n).index.tolist()

global_taxa_order = list(dict.fromkeys(all_top_taxa))

colors = (
    list(plt.cm.tab20.colors) +
    list(plt.cm.Set3.colors) +
    list(plt.cm.Pastel1.colors)
)

global_color_map = {t: colors[i % len(colors)] for i, t in enumerate(global_taxa_order)}
global_color_map["Other"] = "#BDBDBD"


# ===== MAIN LOOP =====
all_records = []

for f in csv_files:
    file_name = os.path.splitext(os.path.basename(f))[0]
    print(f"\nProcessing: {file_name}")

    rel = prepare_16s(pd.read_csv(f, low_memory=False))
    if rel is None:
        print(f"Skipping {file_name}")
        continue

    top_taxa = rel.mean().sort_values(ascending=False).head(top_n).index.tolist()

    group = "Control" if "Control" in file_name else "Pregnant" if "Pregnant" in file_name else "Unknown"

    for taxon in top_taxa:
        all_records.append(pd.DataFrame({
            "dataset": file_name,
            "group": group,
            "taxon": taxon,
            "value": rel[taxon].values
        }))

    rel_top = rel[top_taxa].copy()

    other_taxa = [c for c in rel.columns if c not in top_taxa]
    if other_taxa:
        rel_top["Other"] = rel[other_taxa].sum(axis=1)

    rel_top = rel_top[rel_top.mean(axis=0).sort_values(ascending=False).index]

    if "Bacteria;Firmicutes" in rel_top.columns:
        rel_top = rel_top.sort_values(by="Bacteria;Firmicutes", ascending=False)

    final_order = rel_top.mean(axis=0).sort_values(ascending=False).index.tolist()
    rel_top = rel_top[final_order]

    color_map = {t: global_color_map[t] for t in rel_top.columns}

    # ===== STACKED BAR PLOT =====
    fig, axes = plt.subplots(
        1, 2,
        figsize=(16, 6),
        gridspec_kw={"width_ratios": [5, 1]}
    )

    ax_samples, ax_mean = axes

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

    mean_values = rel_top.mean(axis=0)
    bottom = 0

    for taxon in rel_top.columns:
        ax_mean.bar(
            [0],
            [mean_values[taxon]],
            bottom=bottom,
            width=0.9,
            color=color_map[taxon],
            edgecolor="none"
        )
        bottom += mean_values[taxon]

    ax_mean.set_ylim(0, 100)
    ax_mean.set_xticks([0])
    ax_mean.set_xticklabels(["Mean"], fontsize=11)
    ax_mean.set_title("Average", fontsize=13)
    ax_mean.spines["top"].set_visible(False)
    ax_mean.spines["right"].set_visible(False)
    ax_mean.spines["left"].set_visible(False)
    ax_mean.yaxis.set_visible(False)

    handles = [
        Rectangle((0, 0), 1, 1, color=color_map[t], label=t)
        for t in rel_top.columns
    ]

    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=4,
        frameon=False,
        fontsize=9
    )

    fig.suptitle(file_name, fontsize=16, fontweight="bold")
    plt.tight_layout(rect=[0, 0.12, 1, 1])

    out_path = os.path.join(output_folder, f"{file_name}_phylum_composition.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved stacked plot to: {out_path}")


# ===== COMBINED SPLIT VIOLIN PLOT =====
if all_records:
    violin_df = pd.concat(all_records, ignore_index=True)
    violin_df = violin_df[violin_df["group"].isin(["Control", "Pregnant"])]

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
        for group, side, color in [
            ("Control", "left", control_color),
            ("Pregnant", "right", pregnant_color),
        ]:
            vals = violin_df[
                (violin_df["taxon"] == taxon) &
                (violin_df["group"] == group)
            ]["value"].dropna().values

            if len(vals) > 1:
                parts = ax.violinplot(
                    [vals],
                    positions=[i],
                    widths=0.8,
                    showmeans=False,
                    showmedians=True,
                    showextrema=False
                )

                for body in parts["bodies"]:
                    body.set_facecolor(color)
                    body.set_edgecolor("black")
                    body.set_alpha(0.7)

                    verts = body.get_paths()[0].vertices
                    if side == "left":
                        verts[:, 0] = np.minimum(verts[:, 0], i)
                    else:
                        verts[:, 0] = np.maximum(verts[:, 0], i)

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

    violin_out = os.path.join(output_folder, "violin_plot_16S.png")
    plt.savefig(violin_out, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved combined split violin plot to: {violin_out}")

print("Done.")