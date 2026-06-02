from pathlib import Path
import pandas as pd
import MIPMLP

INPUT_ROOT = Path("datasets_after_yamas/16S")
OUTPUT_FOLDER = Path("datasets_after_MIPMLP")
OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)

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

    # Standardize taxonomy format
    tax["Taxon"] = tax["Taxon"].astype(str)
    tax["Taxon"] = tax["Taxon"].str.replace(",", ";", regex=False)
    tax["Taxon"] = tax["Taxon"].str.replace(r"(^|;)\s*[a-zA-Z]:", r"\1", regex=True)

    otu["taxonomy"] = tax["Taxon"]

    df = otu.T
    df.index.name = "ID"

    df = df[~df.index.astype(str).str.endswith("_2")]
    df.index = df.index.astype(str).str.replace("_1$", "", regex=True)

    df = df.reset_index()

    processed = MIPMLP.preprocess(
        df,
        taxnomy_group="mean",
        normalization="none",
        taxonomy_level=2
    )

    processed.index.name = "SampleID"

    processed = processed.div(processed.sum(axis=1), axis=0) * 100

    mean_abundance = processed.mean(axis=0)
    high_abundance_taxa = mean_abundance[mean_abundance >= 0].index
    processed = processed[high_abundance_taxa].copy()

    processed = processed.div(processed.sum(axis=1), axis=0) * 100

    out_path = OUTPUT_FOLDER / f"{name}_GROUP_16S_formatted.csv"
    processed.to_csv(out_path)

    print(f"Saved: {out_path}")

print("done")