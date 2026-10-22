"""Figure 1E: miMic Control vs Pregnant, using five prepared abundance tables.

Run in mipmlp_env: python fig_1E_16S.py
Taxonomy paths and database suffixes are matched exactly; no invented species.
Abundances are normalized before matching, over each sample's original taxa.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd

# Save figures without opening blocking GUI windows (also works on the server).
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("MPLBACKEND", "Agg")
from updated_miMic_test.miMic_test import apply_mimic

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "datasets_after_MIPMLP"
OUT_DIR = BASE_DIR / "1E_16S/Control_VS_Pregnant_five_projects"
PROJECTS = {
    "Control": ["PRJNA669650", "PRJNA388263", "PRJNA545251"],
    "Pregnant": ["PRJEB58050", "PRJNA1254708"],
}


def load_table(project, group):
    """Load the prepared, nonnegative abundance table; preserve taxonomic ranks."""
    df = pd.read_csv(DATA_DIR / f"{project}_{group}_16S.csv", index_col=0)
    df = df.loc[~df.index.astype(str).str.lower().isin(["taxonomy", "taxon", "taxa"])]
    print(f"{project}: dropping {df.isna().all(axis=1).sum()} entirely missing rows")
    df = df.dropna(how="all")
    df = df.apply(pd.to_numeric, errors="raise")
    if not np.isfinite(df.to_numpy()).all() or (df < 0).any().any():
        raise ValueError(f"{project}: expected finite, nonnegative abundances.")
    # Input headers are semicolon-separated rank names, without prefixes.
    df.columns = [";".join(f"{rank}__{name.strip() or 'unclassified'}"
                          for rank, name in zip("kpcofgs", col.strip().split(";")))
                  for col in df.columns]
    df = df.T.groupby(level=0, sort=False).sum().T
    totals = df.sum(axis=1)
    print(f"{project} ({group}): {len(df)} samples; {(totals == 0).sum()} empty")
    df = df.loc[totals > 0].div(totals[totals > 0], axis=0)
    df.index = [f"{group}_{project}_{sample}" for sample in df.index]
    return df


def main():
    tables, tags = [], []
    for label, (group, projects) in enumerate(PROJECTS.items()):
        for project in projects:
            table = load_table(project, group)
            tables.append(table)
            tags.append(pd.Series(label, index=table.index, name="Tag"))

    shared = sorted(set.intersection(*(set(table.columns) for table in tables)))
    shared = [name for name in shared if ";" in name]  # exclude kingdom-only entries
    # SAMBA expects one column per tree leaf, not both a parent and its children.
    shared = [name for name in shared
              if not any(other.startswith(name + ";") for other in shared)]
    if len(shared) < 2:
        raise ValueError("Too few shared taxa. Check taxonomy compatibility.")
    processed = pd.concat([table.loc[:, shared] for table in tables])
    processed = processed.loc[:, (processed > 0).any(axis=0)]
    nonempty = processed.sum(axis=1) > 0
    print(f"Excluding {(~nonempty).sum()} samples without matched abundance")
    processed = processed.loc[nonempty]
    tag = pd.concat(tags).loc[processed.index].to_frame()
    if not processed.index.is_unique or tag.Tag.nunique() != 2:
        raise ValueError("Need unique sample IDs and samples in both groups.")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    processed.to_csv(OUT_DIR / "processed_for_mimic.csv")
    tag.to_csv(OUT_DIR / "tag.csv")
    print(f"miMic input: {processed.shape}; groups: {tag.Tag.value_counts().to_dict()}")
    np.random.seed(0)
    tax, samba_output = apply_mimic(
        str(OUT_DIR), tag, eval="man", threshold_p=0.05,
        processed=processed, apply_samba=True, save=True,
    )
    if tax not in (None, "nosignificant"):
        apply_mimic(
            str(OUT_DIR), tag, mode="plot", tax=tax, eval="man", sis="fdr_bh",
            samba_output=samba_output, save=True, threshold_p=0.05, THRESHOLD_edge=0.5,
        )
        print(f"Figures saved to {OUT_DIR / 'plots'}")
    else:
        print("No significant taxa; no cladogram generated.")


if __name__ == "__main__":
    main()
