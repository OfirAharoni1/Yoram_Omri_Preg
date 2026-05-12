import pandas as pd
from pathlib import Path

base_path = Path("Yoram_Omri_Preg/datasets_after_yamas/16S/PRJNA669650")

otu_path = base_path / "otu_PRJNA669650.csv"
taxonomy_path = base_path / "taxonomy_PRJNA669650.csv"
preprocess_path = "Yoram_Omri_Preg/datasets_after_yamas/16S/PRJNA669650_Control_16S_for_MIPMLP.csv"

otu = pd.read_csv(otu_path, index_col=0)
tax = pd.read_csv(taxonomy_path, index_col=0)
tax = tax.loc[otu.index]

# Add taxonomy to OTU and transpose
otu["taxonomy"] = tax["Taxon"]

otu_t = otu.T
otu_t.index.name = 'ID'

print(tax["Taxon"].isna().sum())
print(tax[tax["Taxon"].isna()].head())

otu_t.to_csv(preprocess_path)