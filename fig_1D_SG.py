import re
import numpy as np
import pandas as pd
from pathlib import Path

from updated_gimic_package.SAMBA_metric import *

# =====================
# Config
# =====================
CUTOFF = 0.8

DATA_DIR = Path("datasets_after_yamas/SG")
OUT_ROOT = Path("1D_SG")
COMPARISON = "Control_VS_Pregnant"


def group_from_name(path):
    name = path.name.lower()
    if "pregnant" in name or "preg" in name:
        return "Pregnant"
    if "control" in name or "ctrl" in name:
        return "Control"
    return None


def trim_taxonomy_to_species(col):
    col = str(col).strip()

    # remove metadata columns
    if col in ["SampleID", "Group"]:
        return None

    # split full taxonomy path
    parts = re.split(r"[;,]\s*", col)

    keep = []
    allowed = ("k__", "p__", "c__", "o__", "f__", "g__", "s__")

    for p in parts:
        p = p.strip()

        if p.startswith("t__"):
            break

        if p.startswith(allowed):
            keep.append(p)

        if p.startswith("s__"):
            break

    if not keep:
        return None

    return ";".join(keep)


def load_table(path):
    df = pd.read_csv(path)

    if "SampleID" in df.columns:
        df = df.set_index("SampleID")

    df = df.drop(columns=["Group"], errors="ignore")

    new_cols = [trim_taxonomy_to_species(c) for c in df.columns]
    keep_mask = [c is not None for c in new_cols]

    df = df.loc[:, keep_mask]
    df.columns = [c for c in new_cols if c is not None]

    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)

    # if several t__ columns collapse into the same species, sum them
    df = df.groupby(df.columns, axis=1).sum()

    return df


def normalize_for_gimic(df):
    df = np.log10(df + 1)
    df = (df - df.mean()) / df.std()
    return df.replace([np.inf, -np.inf], 0).fillna(0)


# =====================
# Load datasets
# =====================
tables = {"Control": [], "Pregnant": []}

for path in sorted(DATA_DIR.glob("*SG_formatted.csv")):
    group = group_from_name(path)

    if group is None:
        print(f"Skipping: {path.name}")
        continue

    df = load_table(path)
    tables[group].append(df)

    print(f"Loaded {path.name}: {group}, shape after trimming={df.shape}")

if not tables["Control"] or not tables["Pregnant"]:
    raise ValueError("Need at least one Control and one Pregnant SG file.")


# =====================
# Shared taxa only
# =====================
all_tables = tables["Control"] + tables["Pregnant"]
mutual_cols = all_tables[0].columns

for df in all_tables[1:]:
    mutual_cols = mutual_cols.intersection(df.columns)

control = pd.concat([df[mutual_cols] for df in tables["Control"]], axis=0)
preg = pd.concat([df[mutual_cols] for df in tables["Pregnant"]], axis=0)


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
# Run GIMIC-like single comparison plot
# =====================
import matplotlib.pyplot as plt

array_of_imgs, bact_names, ordered_df = micro2matrix(
    processed,
    str(OUT_ROOT / COMPARISON / "2D_images"),
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

# remove empty / bad rows
diff_df = diff_df.replace([np.inf, -np.inf], np.nan).fillna(0)
diff_df = diff_df.loc[(diff_df.abs().sum(axis=1) > 0)]

# Optional visualization filter:
# remove non-informative GTDB-style clades
#bad_patterns = r"^(CFGB|GGB|SGB)\d+$|unclassified|unknown|uncultured|bacterium$"

#diff_df = diff_df[
#    ~diff_df.index.astype(str).str.contains(
#        bad_patterns,
#        case=False,
#        regex=True,
#        na=False
#    )
#]

# keep only top taxa, otherwise the plot is unreadable
TOP_N = 25
diff_df["score"] = diff_df.abs().max(axis=1)
diff_df = diff_df.sort_values("score", ascending=False).head(TOP_N)
diff_df = diff_df.drop(columns="score")

# plot
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
plot_long["size"] = 40 + 1200 * (plot_long["Diff"].abs() / max_abs)

plt.figure(figsize=(10, max(7, TOP_N * 0.35)))

sc = plt.scatter(
    plot_long["x"],
    plot_long["y"],
    s=plot_long["size"],
    c=plot_long["Diff"],
    cmap="bwr_r",
    vmin=-0.2,
    vmax=0.2,
    edgecolors="black",
    linewidths=1.2,
    alpha=0.85
)

plt.xticks(range(len(levels)), levels, rotation=45, ha="right")
plt.yticks(range(len(taxa)), taxa)

plt.grid(True, alpha=0.4)

cbar = plt.colorbar(sc)
cbar.set_label("Sum Value")

plt.title("SG GIMIC: Control vs Pregnant")
plt.xlabel("Taxonomic level")
plt.ylabel("Taxon")

Path("1D_SG/figure_plots").mkdir(exist_ok=True)
plt.tight_layout()
plt.savefig("1D_SG/figure_plots/fig_1D_SG_Control_vs_Pregnant.png", dpi=300, bbox_inches="tight")
plt.close()

print("Done.")
print("Shared taxa:", len(mutual_cols))
print("Plot saved in: 1D_SG/figure_plots/fig_1D_SG_Control_vs_Pregnant.png")