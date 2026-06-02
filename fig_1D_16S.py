"""
Figure 1D - 16S GIMIC analysis

INPUT:
    Original YAMAS 16S datasets:
        datasets_after_yamas/16S/PRJNA669650/
            otu_PRJNA669650.csv
            taxonomy_PRJNA669650.csv

        datasets_after_yamas/16S/PRJNA1254708/
            otu_PRJNA1254708.csv
            taxonomy_PRJNA1254708.csv

    Extra preprocessed 16S tables:
        datasets_after_MIPMLP/omri/omri_stool_Pregnant_16S.csv

OUTPUT:
    1D_16S/Control_VS_Pregnant/for_preprocess.csv
    1D_16S/Control_VS_Pregnant/tag.csv
    1D_16S/figure_plots/fig_1D_16S_Control_vs_Pregnant.png

Notes:
    - The analysis is always between two groups only:
        Control vs Pregnant

    - To add another original YAMAS dataset:
        add its folder name to DATASETS["Control"] or DATASETS["Pregnant"]

    - To add another ready preprocessed table:
        add its full path to EXTRA_TABLES["Control"] or EXTRA_TABLES["Pregnant"]

    - Omri is added to Pregnant as another table in the same group.

    - Tables inside each group are combined using pd.concat(axis=0).
      This means concat is done by samples, not by taxa:
          rows = samples are added one under another
          columns = taxa stay as features

    - Before concat, the code keeps only taxa that are shared across all tables:
          mutual_cols = mutual_cols.intersection(df.columns)

      Therefore, if many projects are added, the number of shared taxa may decrease.
      This is intentional, because GIMIC needs all samples to have the same columns.
"""

import re
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

from updated_gimic_package.SAMBA_metric import *


# =====================
# Config
# =====================
BASE_DIR = Path("/home/aharonox/Yoram_Omri_Preg")

DATA_DIR = BASE_DIR / "datasets_after_yamas/16S"
OUT_ROOT = BASE_DIR / "1D_16S"
COMPARISON = "Control_VS_Pregnant"

CUTOFF = 0.8
TOP_N = 25

# Original YAMAS folders
DATASETS = {
    "Control": ["PRJNA669650"],
    "Pregnant": ["PRJNA1254708"],
}

# Ready preprocessed tables
# Important: for GIMIC we use the full non-normalized table, not phylum.
EXTRA_TABLES = {
    "Control": [],
    "Pregnant": [
        BASE_DIR / "datasets_after_MIPMLP/omri/omri_stool_Pregnant_16S.csv"
    ],
}


# =====================
# Taxonomy cleaning
# =====================
def normalize_taxonomy_name(tax):
    tax = str(tax).strip()

    if tax.lower() in ["nan", "none", "", "0", "unassigned"]:
        return None

    parts = re.split(r"[;,]\s*", tax)

    level_order = ["k", "p", "c", "o", "f", "g", "s"]
    level_map = {
        "d": "k",
        "k": "k",
        "p": "p",
        "c": "c",
        "o": "o",
        "f": "f",
        "g": "g",
        "s": "s",
    }

    tax_dict = {}

    for p in parts:
        p = p.strip()

        if not p or p.lower() in ["nan", "none", "unassigned"]:
            continue

        # format like D_0__Bacteria
        p = re.sub(r"^D_0__", "k__", p)
        p = re.sub(r"^D_1__", "p__", p)
        p = re.sub(r"^D_2__", "c__", p)
        p = re.sub(r"^D_3__", "o__", p)
        p = re.sub(r"^D_4__", "f__", p)
        p = re.sub(r"^D_5__", "g__", p)
        p = re.sub(r"^D_6__", "s__", p)

        # format like k:Bacteria / p:Firmicutes
        m = re.match(r"^([dkpcfgos]):(.+)$", p)
        if m:
            code = level_map.get(m.group(1))
            value = m.group(2).strip()
            p = f"{code}__{value}"

        # format like d__Bacteria
        if p.startswith("d__"):
            p = "k__" + p.replace("d__", "", 1)

        # format like k__Bacteria
        m = re.match(r"^([kpcofgs])__(.*)$", p)
        if m:
            code = m.group(1)
            value = m.group(2).strip()

            if value == "":
                value = "unclassified"

            value = value.strip("[]")
            tax_dict[code] = value

    # If taxonomy has no prefixes, e.g.
    # Bacteria;Firmicutes;Bacilli;...
    # assign levels by order.
    if not tax_dict:
        plain_parts = [
            p.strip().strip("[]")
            for p in parts
            if p.strip() and p.strip().lower() not in ["nan", "none", "unassigned"]
        ]

        for code, value in zip(level_order, plain_parts):
            tax_dict[code] = value

    if not tax_dict:
        return None

    # Fill missing levels so GIMIC always gets 7 levels
    last_value = "unclassified"
    keep = []

    for code in level_order:
        value = tax_dict.get(code)

        if value is None or value == "":
            value = f"{last_value}_unclassified"
        else:
            last_value = value

        keep.append(f"{code}__{value}")

    return ";".join(keep)


# =====================
# Load original YAMAS OTU + taxonomy
# =====================
def load_otu_tax_for_gimic(folder):
    otu_path = next(folder.glob("otu_*.csv"))
    tax_path = next(folder.glob("taxonomy_*.csv"))

    otu = pd.read_csv(otu_path, index_col=0)
    tax = pd.read_csv(tax_path, index_col=0)

    # safer than tax.loc[otu.index] in case some OTUs are missing
    tax = tax.reindex(otu.index)

    new_cols = [normalize_taxonomy_name(x) for x in tax["Taxon"]]
    keep = [x is not None for x in new_cols]

    otu = otu.loc[keep]
    new_cols = [x for x in new_cols if x is not None]

    # OTU table: rows = OTUs, columns = samples
    # GIMIC input: rows = samples, columns = taxonomy
    df = otu.T
    df.columns = new_cols
    df.index.name = "SampleID"

    # paired-end: keep _1, drop _2
    df = df[~df.index.astype(str).str.endswith("_2")]
    df.index = df.index.astype(str).str.replace(r"_1$", "", regex=True)

    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)

    # Collapse identical full taxonomy names
    df = df.T.groupby(level=0).sum().T

    return df


# =====================
# Load ready preprocessed table
# =====================
def load_preprocessed_for_gimic(path):
    """
    Input format:
        rows = samples
        columns = taxonomy/taxa

    Output:
        rows = samples
        columns = cleaned full taxonomy names
    """
    df = pd.read_csv(path, index_col=0)

    df.columns = df.columns.astype(str).str.strip()
    df.columns = df.columns.str.replace(r"\s*;\s*", ";", regex=True)

    # remove non-informative columns if they exist
    df = df.drop(columns=["Bacteria", "Unassigned"], errors="ignore")

    # remove possible taxonomy row if it exists
    df = df[~df.index.astype(str).str.lower().str.contains("tax")]

    new_cols = [normalize_taxonomy_name(c) for c in df.columns]
    keep = [c is not None for c in new_cols]

    df = df.loc[:, keep]
    df.columns = [c for c in new_cols if c is not None]

    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)

    # Collapse identical full taxonomy names
    df = df.T.groupby(level=0).sum().T

    return df


# =====================
# Normalize for GIMIC
# =====================
def normalize_for_gimic(df):
    df = np.log10(df + 1)

    std = df.std().replace(0, np.nan)
    df = (df - df.mean()) / std

    return df.replace([np.inf, -np.inf], 0).fillna(0)


# =====================
# Load datasets
# =====================
tables = {
    "Control": [],
    "Pregnant": [],
}

# Load original YAMAS folders
for group, folder_names in DATASETS.items():
    for folder_name in folder_names:
        folder = DATA_DIR / folder_name

        if not folder.exists():
            print(f"Skipping missing folder: {folder}")
            continue

        df = load_otu_tax_for_gimic(folder)
        tables[group].append(df)

        print(f"Loaded {folder_name}: {group}, shape={df.shape}")


# Load extra preprocessed tables
for group, paths in EXTRA_TABLES.items():
    for path in paths:
        if not path.exists():
            print(f"Skipping missing file: {path}")
            continue

        df = load_preprocessed_for_gimic(path)
        tables[group].append(df)

        print(f"Loaded {path.name}: {group}, shape={df.shape}")


if not tables["Control"] or not tables["Pregnant"]:
    raise ValueError("Need at least one Control and one Pregnant 16S dataset.")


# =====================
# Shared taxa only
# =====================
all_tables = tables["Control"] + tables["Pregnant"]

mutual_cols = all_tables[0].columns
for df in all_tables[1:]:
    mutual_cols = mutual_cols.intersection(df.columns)

mutual_cols = sorted(mutual_cols)

print("Shared taxa:", len(mutual_cols))

if len(mutual_cols) == 0:
    raise ValueError("No shared taxa between Control and Pregnant.")


# Combine all projects inside each group
control = pd.concat([df[mutual_cols] for df in tables["Control"]], axis=0)
preg = pd.concat([df[mutual_cols] for df in tables["Pregnant"]], axis=0)

print("Combined Control shape:", control.shape)
print("Combined Pregnant shape:", preg.shape)


# =====================
# Normalize + tag
# =====================
control = normalize_for_gimic(control)
preg = normalize_for_gimic(preg)

control.index = [f"Control_{i}" for i in range(len(control))]
preg.index = [f"Pregnant_{i}" for i in range(len(preg))]

processed = pd.concat([control, preg], axis=0)

tag = pd.DataFrame(index=processed.index)
tag["Tag"] = [0 if x.startswith("Control") else 1 for x in processed.index]


# =====================
# Save GIMIC input
# =====================
comparison_dir = OUT_ROOT / COMPARISON
comparison_dir.mkdir(parents=True, exist_ok=True)

processed.to_csv(comparison_dir / "for_preprocess.csv")
tag.to_csv(comparison_dir / "tag.csv")

print("Saved:")
print(comparison_dir / "for_preprocess.csv")
print(comparison_dir / "tag.csv")


# =====================
# Run GIMIC-like comparison
# =====================
array_of_imgs, bact_names, ordered_df = micro2matrix(
    processed,
    str(comparison_dir / "2D_images"),
    save=False
)

tag0 = tag[tag["Tag"] == 0]
tag1 = tag[tag["Tag"] == 1]

imgs_names = list(ordered_df.index)

img0_index = [imgs_names.index(x) for x in tag0.index if x in imgs_names]
img1_index = [imgs_names.index(x) for x in tag1.index if x in imgs_names]

imgs0 = array_of_imgs[img0_index]
imgs1 = array_of_imgs[img1_index]

diff_df = apply_class_analysis(
    imgs1,
    imgs0,
    CUTOFF=CUTOFF,
    bact_names=bact_names
)

diff_df = diff_df.replace([np.inf, -np.inf], np.nan).fillna(0)
diff_df = diff_df.loc[diff_df.abs().sum(axis=1) > 0]


# =====================
# Filter + top taxa
# =====================
bad_patterns = r"unclassified|unknown|uncultured|metagenome|bacterium$"

diff_df = diff_df[
    ~diff_df.index.astype(str).str.contains(
        bad_patterns,
        case=False,
        regex=True,
        na=False
    )
]

if diff_df.empty:
    raise ValueError("No taxa left after filtering unclassified/unknown labels.")

diff_df["score"] = diff_df.abs().max(axis=1)
diff_df = diff_df.sort_values("score", ascending=False).head(TOP_N)
diff_df = diff_df.drop(columns="score")

if diff_df.empty:
    raise ValueError("No differential taxa found after GIMIC analysis.")


# =====================
# Plot
# =====================
levels = list(diff_df.columns)
taxa = list(diff_df.index)

plot_long = diff_df.reset_index().melt(
    id_vars="index",
    var_name="Level",
    value_name="Diff"
).rename(columns={"index": "Taxon"})

plot_long["x"] = plot_long["Level"].map({lvl: i for i, lvl in enumerate(levels)})
plot_long["y"] = plot_long["Taxon"].map({tax: i for i, tax in enumerate(taxa)})

max_abs = plot_long["Diff"].abs().max()

if max_abs == 0 or pd.isna(max_abs):
    plot_long["size"] = 80
    color_lim = 1
else:
    plot_long["size"] = 40 + 1200 * (plot_long["Diff"].abs() / max_abs)
    color_lim = max_abs

plt.figure(figsize=(10, max(7, len(taxa) * 0.35)))

sc = plt.scatter(
    plot_long["x"],
    plot_long["y"],
    s=plot_long["size"],
    c=plot_long["Diff"],
    cmap="bwr_r",
    vmin=-color_lim,
    vmax=color_lim,
    edgecolors="black",
    linewidths=1.2,
    alpha=0.85
)

plt.xticks(range(len(levels)), levels, rotation=45, ha="right")
plt.yticks(range(len(taxa)), taxa)

plt.grid(True, alpha=0.4)

cbar = plt.colorbar(sc)
cbar.set_label("GIMIC score: Pregnant - Control")

plt.title("16S GIMIC: Control vs Pregnant")
plt.xlabel("Taxonomic level")
plt.ylabel("Taxon")

plot_dir = OUT_ROOT / "figure_plots"
plot_dir.mkdir(parents=True, exist_ok=True)

out_plot = plot_dir / "fig_1D_16S_Control_vs_Pregnant.png"

plt.tight_layout()
plt.savefig(out_plot, dpi=300, bbox_inches="tight")
plt.close()

print("Done.")
print("Shared taxa:", len(mutual_cols))
print("Taxa shown:", len(taxa))
print("Color limit:", color_lim)
print("Plot saved in:", out_plot)