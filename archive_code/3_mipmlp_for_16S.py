import glob, os, pandas as pd, MIPMLP

INPUT_FOLDER = "datasets_after_yamas/16S"
OUTPUT_FOLDER = "datasets_after_MIPMLP"
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

for f in glob.glob(f"{INPUT_FOLDER}/*.csv"):
    df = pd.read_csv(f)

    # Remove rows where the "ID" column ends with "_2"
    df = df[~df["ID"].astype(str).str.endswith("_2")]

    processed = MIPMLP.preprocess(
        df,
        taxnomy_group="mean",
        normalization="none",
        taxonomy_level=2
    )
    processed.index.name = "SampleID"

    out = os.path.join(OUTPUT_FOLDER, os.path.basename(f).replace(".csv", "_MIPMLP.csv"))
    processed.to_csv(out)

print("done")