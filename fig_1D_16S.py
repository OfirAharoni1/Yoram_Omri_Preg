import re
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

from updated_gimic_package.SAMBA_metric import *

# =====================
# Config
# =====================
CUTOFF = 0.8

DATA_DIR = Path("datasets_after_yamas/16S")
OUT_ROOT = Path("1D_16S")
COMPARISON = "Control_VS_Pregnant"
TOP_N = 25


def group_from_name(path):
    name = path.name.lower()
    if "pregnant" in name or "preg" in name:
        return "Pregnant"
    if "control" in name or "ctrl" in name:
        return "Control"
    return None


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

        m = re.match(r"^([kpcofgs])__(.*)$", p)
        if m:
            code = m.group(1)
            value = m.group(2).strip()

            if value == "":
                value = "unclassified"

            tax_dict[code] = value

    if not tax_dict:
        return None

    # fill missing levels so GIMIC always gets 7 levels
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


def load_clean_16s_for_gimic(path):
    df = pd.read_csv(path, index_col=0, sep=None, engine="python")

    # find taxonomy row
    tax_mask = df.index.astype(str).str.lower().str.contains("tax")

    if tax_mask.any():
        tax_row = df.loc[tax_mask].iloc[0]
        df = df.loc[~tax_mask]
    else:
        tax_row = df.iloc[-1]
        df = df.iloc[:-1]

    # paired-end: keep only _1, drop _2
    idx = df.index.astype(str)
    if (idx.str.endswith("_1") | idx.str.endswith("_2")).any():
        df = df.loc[~idx.str.endswith("_2")]
        df.index = df.index.astype(str).str.replace(r"_1$", "", regex=True)

    # use taxonomy row as column names
    new_cols = [normalize_taxonomy_name(x) for x in tax_row.values]
    keep_mask = [x is not None for x in new_cols]

    df = df.loc[:, keep_mask]
    df.columns = [x for x in new_cols if x is not None]

    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)

    # collapse features with same taxonomy
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

for path in sorted(DATA_DIR.glob("*16S*.csv")):
    group = group_from_name(path)

    if group is None:
        print(f"Skipping: {path.name}")
        continue

    df = load_clean_16s_for_gimic(path)
    tables[group].append(df)

    print(f"Loaded {path.name}: {group}, shape={df.shape}")

if not tables["Control"] or not tables["Pregnant"]:
    raise ValueError("Need at least one Control and one Pregnant 16S file.")


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

diff_df = diff_df.replace([np.inf, -np.inf], np.nan).fillna(0)
diff_df = diff_df.loc[(diff_df.abs().sum(axis=1) > 0)]

# Optional visualization filter
# bad_patterns = r"unclassified|unknown|uncultured|metagenome|bacterium$"
# diff_df = diff_df[
#     ~diff_df.index.astype(str).str.contains(
#         bad_patterns,
#         case=False,
#         regex=True,
#         na=False
#     )
# ]

diff_df["score"] = diff_df.abs().max(axis=1)
diff_df = diff_df.sort_values("score", ascending=False).head(TOP_N)
diff_df = diff_df.drop(columns="score")


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
cbar.set_label("GIMIC score: Pregnant - Control")

plt.title("16S GIMIC: Control vs Pregnant")
plt.xlabel("Taxonomic level")
plt.ylabel("Taxon")

plot_dir = OUT_ROOT / "figure_plots"
plot_dir.mkdir(parents=True, exist_ok=True)

plt.tight_layout()
plt.savefig(plot_dir / "fig_1D_16S_Control_vs_Pregnant.png", dpi=300, bbox_inches="tight")
plt.close()

print("Done.")
print("Shared taxa:", len(mutual_cols))
print("Plot saved in:", plot_dir / "fig_1D_16S_Control_vs_Pregnant.png")