from pathlib import Path
import re

import pandas as pd
import MIPMLP

INPUT_ROOT = Path("datasets_after_yamas/16S")
OUTPUT_FOLDER = Path("datasets_after_MIPMLP")
OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)


RANKS = ("k", "p", "c", "o", "f", "g", "s")
D_RANKS = {f"D_{i}": rank for i, rank in enumerate(RANKS)}
MISSING_TAXON_VALUES = {"", "nan", "none", "unassigned", "unknown"}


def normalize_taxonomy_name(taxon):
    """Return one canonical taxonomy string for the formats used by YAMAS.

    Rank prefixes are converted to ``k__/p__/...`` and numeric confidence
    annotations such as ``Bacteria(0.71)`` are removed. Empty ranked levels
    are retained so that a missing level cannot shift all subsequent ranks.
    """
    value = str(taxon).strip()
    if value.lower() in MISSING_TAXON_VALUES:
        return "Unassigned"

    parts = re.split(r"\s*[;,|]\s*", value)
    ranked = {}
    plain = []

    for part in parts:
        part = part.strip()
        if not part or part.lower() in MISSING_TAXON_VALUES:
            continue

        match = re.match(
            r"^(?:(D_[0-6])__|([dkpcofgs])(?:__|:))\s*(.*)$",
            part,
            flags=re.IGNORECASE,
        )

        if match:
            rank = D_RANKS.get(
                (match.group(1) or "").upper(),
                (match.group(2) or "").lower(),
            )
            name = match.group(3).strip()
        else:
            rank = None
            name = part

        # Some classifiers append a per-rank confidence score to the name.
        # Restrict removal to a trailing number in [0, 1], so meaningful
        # parenthesized text in a taxon name is preserved.
        name = re.sub(
            r"\s*\((?:0(?:\.\d+)?|1(?:\.0+)?)\)\s*$",
            "",
            name,
        ).strip()

        if rank:
            ranked[rank] = name
        elif name:
            plain.append(name)

    if ranked:
        deepest_rank = max(RANKS.index(rank) for rank in ranked)
        return ";".join(
            f"{rank}__{ranked.get(rank, '')}"
            for rank in RANKS[:deepest_rank + 1]
        )

    return ";".join(plain) if plain else "Unassigned"


for folder in INPUT_ROOT.iterdir():
    if not folder.is_dir():
        continue

    name = folder.name
    otu_path = next(folder.glob("otu_*.csv"), None)
    tax_path = next(folder.glob("taxonomy_*.csv"), None)

    if otu_path is None or tax_path is None:
        print(f"Skipping {name}")
        continue

    print(f"Processing {name}")

    otu = pd.read_csv(otu_path, index_col=0)
    tax = pd.read_csv(tax_path, index_col=0)

    tax = tax.loc[otu.index]

    # Standardize ranks and remove classifier confidence annotations before
    # MIPMLP groups features. This makes columns comparable across projects.
    tax["Taxon"] = tax["Taxon"].map(normalize_taxonomy_name)

    otu["taxonomy"] = tax["Taxon"]

    df = otu.T
    df.index.name = "ID"

    df = df[~df.index.astype(str).str.endswith("_2")]
    df.index = df.index.astype(str).str.replace("_1$", "", regex=True)

    df = df.reset_index()

    processed_norm = MIPMLP.preprocess(
        df,
        taxnomy_group="mean",
        normalization="none",
        rare_bacteria_threshold=None,
        taxonomy_level=2
    )

    processed_norm.index.name = "SampleID"

    processed_norm = processed_norm.div(processed_norm.sum(axis=1), axis=0) * 100

    mean_abundance = processed_norm.mean(axis=0)
    high_abundance_taxa = mean_abundance[mean_abundance >= 0].index
    processed_norm = processed_norm[high_abundance_taxa].copy()

    processed_norm = processed_norm.div(processed_norm.sum(axis=1), axis=0) * 100

    out_path = OUTPUT_FOLDER / f"{name}_GROUP_16S_formatted.csv"
    processed_norm.to_csv(out_path)

    # save one without normalization
    processed = MIPMLP.preprocess(
        df,
        normalization="none",
        rare_bacteria_threshold=None
    )

    processed.index.name = "SampleID"

    processed = processed.div(processed.sum(axis=1), axis=0) * 100

    mean_abundance = processed.mean(axis=0)
    high_abundance_taxa = mean_abundance[mean_abundance >= 0].index
    processed_ = processed[high_abundance_taxa].copy()
    processed = processed.div(processed.sum(axis=1), axis=0) * 100

    out_path = OUTPUT_FOLDER / f"{name}_GROUP_16S.csv"
    processed.to_csv(out_path)

    print(f"Saved: {out_path}")

print("done")
