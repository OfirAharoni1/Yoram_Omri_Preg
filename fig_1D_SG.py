import os
import numpy as np
import pandas as pd
from pathlib import Path

from updated_gimic_package.SAMBA_metric import *

# =====================
# Config
# =====================
CLASS = False
CUTOFF = 0.8

DATA_DIR = Path("datasets_after_yamas/SG")
OUT_ROOT = Path("SG_for_gimic")
COMPARISON = "Control_VS_Pregnant"


def group_from_name(path):
    name = path.name.lower()
    if "pregnant" in name or "preg" in name:
        return "Pregnant"
    if "control" in name or "ctrl" in name:
        return "Control"
    return None


def load_table(path):
    df = pd.read_csv(path)

    if "SampleID" in df.columns:
        df = df.set_index("SampleID")

    df = df.drop(columns=["Group"], errors="ignore")
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)

    return df


def normalize_for_gimic(df):
    df = np.log10(df + 1)
    df = (df - df.mean()) / df.std()
    return df.replace([np.inf, -np.inf], 0).fillna(0)


# =====================
# Load all datasets
# =====================
tables = {"Control": [], "Pregnant": []}

for path in sorted(DATA_DIR.glob("*SG_formatted.csv")):
    group = group_from_name(path)

    if group is None:
        print(f"Skipping, unknown group: {path.name}")
        continue

    df = load_table(path)
    tables[group].append(df)

    print(f"Loaded {path.name}: {group}, shape={df.shape}")

if not tables["Control"] or not tables["Pregnant"]:
    raise ValueError("Need at least one Control and one Pregnant SG file.")


# =====================
# Keep only shared bacteria
# =====================
all_tables = tables["Control"] + tables["Pregnant"]
mutual_cols = all_tables[0].columns

for df in all_tables[1:]:
    mutual_cols = mutual_cols.intersection(df.columns)

tables["Control"] = [df[mutual_cols] for df in tables["Control"]]
tables["Pregnant"] = [df[mutual_cols] for df in tables["Pregnant"]]


# =====================
# Merge by group
# =====================
control = pd.concat(tables["Control"], axis=0)
preg = pd.concat(tables["Pregnant"], axis=0)

control = normalize_for_gimic(control)
preg = normalize_for_gimic(preg)

control.index = [f"Control_{i}" for i in range(len(control))]
preg.index = [f"Pregnant_{i}" for i in range(len(preg))]

processed = pd.concat([control, preg], axis=0)

tag = pd.DataFrame(index=processed.index)
tag["Tag"] = [0 if x.startswith("Control") else 1 for x in processed.index]


# =====================
# Save GIMIC inputs
# =====================
comparison_dir = OUT_ROOT / COMPARISON
comparison_dir.mkdir(parents=True, exist_ok=True)

processed.to_csv(comparison_dir / "for_preprocess.csv")
tag.to_csv(comparison_dir / "tag.csv")


# =====================
# Run GIMIC
# =====================
image_folder = OUT_ROOT / "2D_images_SG"

array_of_imgs, bact_names, ordered_df = micro2matrix(
    processed,
    str(image_folder),
    save=False
)

DM = build_SAMBA_distance_matrix(
    str(image_folder),
    imgs=array_of_imgs,
    ordered_df=ordered_df,
    bact_names=bact_names,
    class_=CLASS
)

apply_meta_analysis(
    str(OUT_ROOT),
    [COMPARISON],
    CUTOFF,
    "SG"
)

print("Done.")
print("Shared taxa:", len(mutual_cols))
print("Saved in:", OUT_ROOT)