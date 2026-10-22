"""Figure 1E: miMic Control vs Pregnant, using four shotgun abundance tables.

Run in mipmlp_env: python fig_1E_SG.py
Taxonomy paths and database suffixes are matched exactly; no invented species.
Abundances are normalized before matching, over each sample's species (not the repeated ancestor abundances).
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
DATA_DIR = BASE_DIR / "datasets_after_yamas/SG"
OUT_DIR = BASE_DIR / "1E_SG/Control_VS_Pregnant_four_projects"
PROJECTS = {
    "Control": ["PRJEB37731", "PRJNA48479"],
    "Pregnant": ["ERP020710", "PRJNA1247940"],
}


def load_table(project, group):
    """Load the prepared, nonnegative abundance table; preserve taxonomic ranks."""
    df = pd.read_csv(DATA_DIR / f"{project}_{group}_SG_formatted.csv", index_col=0)
    df = df.loc[~df.index.astype(str).str.lower().isin(["taxonomy", "taxon", "taxa"])]
    print(f"{project}: dropping {df.isna().all(axis=1).sum()} entirely missing rows")
    df = df.dropna(how="all")
    df = df.apply(pd.to_numeric, errors="raise")
    if not np.isfinite(df.to_numpy()).all() or (df < 0).any().any():
        raise ValueError(f"{project}: expected finite, nonnegative abundances.")
    # Shotgun tables contain cumulative abundance at every rank. Use species
    # only, so parent abundances are not counted again during normalization.
    df.columns = df.columns.str.replace(",", ";", regex=False)
    species = [col for col in df if len(col.split(";")) == 7
               and all(part.startswith(rank + "__") and len(part) > 3
                       for rank, part in zip("kpcofgs", col.split(";")))]
    df = df.loc[:, species]
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
