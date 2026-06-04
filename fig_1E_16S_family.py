"""
Figure 1E - 16S miMic analysis

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
    1E_16S/Control_VS_Pregnant/
        processed_for_mimic.csv
        tag.csv
        miMic / SAMBA outputs and plots from apply_mimic

Notes:
    - The analysis is between two groups only:
        Control vs Pregnant

    - Original YAMAS datasets are loaded from raw OTU + taxonomy files.
    - Omri is loaded from the ready preprocessed full 16S table.
    - Do NOT use phylum-normalized tables for miMic.
    - Tables inside each group are combined using pd.concat(axis=0):
        rows = samples are added one under another
        columns = taxa stay as features

    - Before concat, the code keeps only taxa shared across all tables.
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd

from updated_miMic_test.miMic_test import apply_mimic


# =====================
# Config
# =====================

BASE_DIR = Path("/home/aharonox/Yoram_Omri_Preg")

DATA_DIR = BASE_DIR / "datasets_after_yamas/16S"
OUT_ROOT = BASE_DIR / "1E_16S"
COMPARISON = "Control_VS_Pregnant_family"

# Original YAMAS folders
DATASETS = {
    "Control": ["PRJNA669650"],
    "Pregnant": ["PRJNA1254708"],
}

TAXONOMY_LEVEL = 5  # 6=genus, 5=family, 7=species/full

# Ready preprocessed full 16S tables
# Important: use full non-phylum table, not phylum_normalized.
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

    # If taxonomy has no prefixes:
    # Bacteria;Firmicutes;Bacilli;...
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

    # Fill missing levels so all columns have comparable full taxonomy names
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

def load_otu_tax_for_mimic(folder):
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
    # miMic input: rows = samples, columns = taxonomy
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

def load_preprocessed_for_mimic(path):
    """
    Input:
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
# Normalize for miMic
# =====================

def normalize_for_mimic(df):
    """
    Convert counts to relative abundance per sample.

    We do this here because the original YAMAS tables are count tables,
    and different datasets can have different sequencing depths.

    This is not phylum normalization.
    It keeps the full taxa columns.
    """
    row_sums = df.sum(axis=1).replace(0, np.nan)
    df = df.div(row_sums, axis=0)

    return df.replace([np.inf, -np.inf], 0).fillna(0)


# =====================
# Collapse taxonomy level
# =====================

def collapse_to_taxonomy_level(df, level):
    """
    Collapse taxonomy columns to a chosen level before taking shared taxa.

    level:
        7 = species/full taxonomy
        6 = genus
        5 = family
        4 = order
        3 = class
        2 = phylum

    The output still keeps 7 taxonomy fields so miMic can work normally.
    Lower levels are filled as *_unclassified after the chosen level.
    """
    level = int(level)

    if level >= 7:
        return df.T.groupby(level=0).sum().T

    if level < 2:
        raise ValueError("TAXONOMY_LEVEL should be between 2 and 7")

    level_codes = ["k", "p", "c", "o", "f", "g", "s"]
    new_cols = []

    for col in df.columns:
        parts = str(col).split(";")
        values = {}

        for p in parts:
            p = p.strip()
            m = re.match(r"^([kpcofgs])__(.*)$", p)
            if m:
                code = m.group(1)
                value = m.group(2).strip()
                values[code] = value if value else "unclassified"

        last_value = "unclassified"
        keep = []

        for i, code in enumerate(level_codes):
            if i < level:
                value = values.get(code)
                if value is None or value == "":
                    value = f"{last_value}_unclassified"
                else:
                    last_value = value
            else:
                value = f"{last_value}_unclassified"

            keep.append(f"{code}__{value}")

        new_cols.append(";".join(keep))

    out = df.copy()
    out.columns = new_cols

    # Collapse identical taxa after truncating level
    out = out.T.groupby(level=0).sum().T

    return out


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

        df = load_otu_tax_for_mimic(folder)
        tables[group].append(df)

        print(f"Loaded {folder_name}: {group}, shape={df.shape}")


# Load extra preprocessed tables
for group, paths in EXTRA_TABLES.items():
    for path in paths:
        if not path.exists():
            print(f"Skipping missing file: {path}")
            continue

        df = load_preprocessed_for_mimic(path)
        tables[group].append(df)

        print(f"Loaded {path.name}: {group}, shape={df.shape}")


if not tables["Control"] or not tables["Pregnant"]:
    raise ValueError("Need at least one Control and one Pregnant 16S dataset.")


# Collapse taxonomy before intersecting shared taxa
if TAXONOMY_LEVEL < 7:
    for group in tables:
        tables[group] = [
            collapse_to_taxonomy_level(df, TAXONOMY_LEVEL)
            for df in tables[group]
        ]

    print(f"\nCollapsed taxonomy to level {TAXONOMY_LEVEL} before shared taxa.")


# =====================
# Shared taxa only
# =====================

all_tables = tables["Control"] + tables["Pregnant"]

mutual_cols = all_tables[0].columns
for df in all_tables[1:]:
    mutual_cols = mutual_cols.intersection(df.columns)

mutual_cols = sorted(mutual_cols)

print("\nShared taxa:", len(mutual_cols))

if len(mutual_cols) == 0:
    raise ValueError("No shared taxa between Control and Pregnant.")


# =====================
# Combine groups
# =====================

control = pd.concat([df[mutual_cols] for df in tables["Control"]], axis=0)
pregnant = pd.concat([df[mutual_cols] for df in tables["Pregnant"]], axis=0)

print("Combined Control shape:", control.shape)
print("Combined Pregnant shape:", pregnant.shape)


# =====================
# Normalize + tag
# =====================

control = normalize_for_mimic(control)
pregnant = normalize_for_mimic(pregnant)

control.index = [f"Control_{i}" for i in range(len(control))]
pregnant.index = [f"Pregnant_{i}" for i in range(len(pregnant))]

processed = pd.concat([control, pregnant], axis=0)

# remove taxa that are zero in all samples
processed = processed.loc[:, processed.sum(axis=0) > 0]

tag = pd.DataFrame(index=processed.index)
tag["Tag"] = [0 if x.startswith("Control") else 1 for x in processed.index]

# make sure tag and processed have the same order
tag = tag.loc[processed.index]


# =====================
# Save miMic input
# =====================

comparison_dir = OUT_ROOT / COMPARISON
comparison_dir.mkdir(parents=True, exist_ok=True)

processed_path = comparison_dir / "processed_for_mimic.csv"
tag_path = comparison_dir / "tag.csv"

processed.to_csv(processed_path)
tag.to_csv(tag_path)

print("\nSaved:")
print(processed_path)
print(tag_path)

print("\nFinal processed shape:", processed.shape)
print("Final tag shape:", tag.shape)
print("Control samples:", (tag["Tag"] == 0).sum())
print("Pregnant samples:", (tag["Tag"] == 1).sum())
print("Final taxa:", processed.shape[1])


# =====================
# Run miMic
# =====================

taxonomy_selected, samba_output = apply_mimic(
    str(comparison_dir),
    tag,
    eval="man",
    threshold_p=0.05,
    processed=processed,
    apply_samba=True,
    save=True,
)

if taxonomy_selected is not None:
    apply_mimic(
        str(comparison_dir),
        tag,
        mode="plot",
        tax=taxonomy_selected,
        eval="man",
        sis="fdr_bh",
        samba_output=samba_output,
        save=True,
        threshold_p=0.05,
        THRESHOLD_edge=0.5,
    )
else:
    print("No taxonomy selected by miMic, skipping plot.")

print("\nDone.")