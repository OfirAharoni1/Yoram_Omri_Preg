'''
Per-project microbe-metadata correlation analysis.

Goal:
For each dataset/project separately, create a bar plot showing the number of
significant correlations (FDR-adjusted p < 0.05) between metadata features and
microbial taxa.

Input format:
1. Microbiome table:
   - rows = samples
   - columns = microbial taxa / bacteria
   - first column = sample IDs / index

2. Metadata table:
   - rows = samples
   - columns = metadata variables
   - first column = sample IDs / index

Output:
For each project, the script saves:
- all_microbe_metadata_correlations.csv
- significant_correlation_counts.csv
- number_of_significant_correlations.png
- number_of_significant_correlations.pdf

It also saves one combined table across all projects.
'''

import os
import re
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr


# ============================================================
# INPUTS
# ============================================================
# Add every dataset/project here.
# Each project gets its own independent correlation analysis and its own plot.
PROJECTS = [
    {
        "project_name": "omri_stool_Pregnant_16S",
        "microbiome_file": "/home/aharonox/Yoram_Omri_Preg/datasets_after_MIPMLP/omri/omri_stool_Pregnant_16S.csv",
        "metadata_file": "/home/aharonox/Yoram_Omri_Preg/datasets_after_MIPMLP/omri/omri_metadata_Pregnant_expanded_to_microbiome_samples.csv",
    },
    {
         "project_name": "PRJNA669650_Control_16S",
         "microbiome_file": "/home/aharonox/Yoram_Omri_Preg/datasets_after_MIPMLP/PRJNA669650_Control_16S.csv",
         "metadata_file": "/home/aharonox/Yoram_Omri_Preg/datasets_after_MIPMLP/PRJNA669650_Control_16S_metadata.csv",
    },
    {
         "project_name": "PRJNA1254708_Pregnant_16S",
         "microbiome_file": "/home/aharonox/Yoram_Omri_Preg/datasets_after_MIPMLP/PRJNA1254708_Pregnant_16S.csv",
         "metadata_file": "", #doesnt have metadata
    },
    {
         "project_name": "ERP020710_Pregnant_SG",
         "microbiome_file": "/home/aharonox/Yoram_Omri_Preg/datasets_after_yamas/SG/ERP020710_Pregnant_SG_formatted.csv",
         "metadata_file": "/home/aharonox/Yoram_Omri_Preg/datasets_after_yamas/SG/ERP020710_metadata_filtered.csv",
    },
    {
         "project_name": "PRJEB37731_Control_SG",
         "microbiome_file": "/home/aharonox/Yoram_Omri_Preg/datasets_after_yamas/SG/PRJEB37731_Control_SG_formatted.csv",
         "metadata_file": "/home/aharonox/Yoram_Omri_Preg/datasets_after_yamas/SG/PRJEB37731_metadata_filtered.csv",
    },
    {
         "project_name": "PRJNA48479_Control_SG",
         "microbiome_file": "/home/aharonox/Yoram_Omri_Preg/datasets_after_yamas/SG/PRJNA48479_Control_SG_formatted.csv",
         "metadata_file": "",  #doesnt have metadata
    },
    {
         "project_name": "PRJNA1247940_Pregnant_SG",
         "microbiome_file": "/home/aharonox/Yoram_Omri_Preg/datasets_after_yamas/SG/PRJNA1247940_Pregnant_SG_formatted.csv",
         "metadata_file": "/home/aharonox/Yoram_Omri_Preg/datasets_after_yamas/SG/PRJNA1247940_metadata_filtered.csv",
    },
    
]

OUTPUT_DIR = "fig_2A"

# Analysis settings
NORMALIZE_MICROBIOME_ROWS = True
P_VALUE_THRESHOLD = 0.05
MIN_SAMPLES_PER_TEST = 10
MIN_NONZERO_SAMPLES_PER_TAXON = 3
MAX_CATEGORIES_FOR_METADATA = 10
TOP_N_METADATA_TO_PLOT = 20

# What should the bar plot count?
# "metadata_feature" = numeric columns as-is, categorical columns as dummy variables such as sex_Female.
# "metadata_column"  = collapse categorical dummy variables back to the original metadata column name.
PLOT_BY = "metadata_feature"

# Metadata columns to exclude from correlation analysis
EXCLUDE_METADATA_COLUMNS = {
    "sample_id", "SampleID", "sample", "Sample", "sample_name", "SampleName",
    "subject_id", "SubjectID", "participant_id", "ParticipantID",
    "run", "Run", "project", "Project", "BioProject", "bioproject",
    "accession", "Accession", "ACC", "sra", "SRA",
}


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def safe_name(x):
    """Create a safe name for files/folders."""
    x = str(x)
    x = re.sub(r"[^A-Za-z0-9_.-]+", "_", x)
    x = x.strip("_")
    return x if x else "project"


def clean_label(x):
    """Create a cleaner label for plotting."""
    x = str(x)
    x = x.replace("_", " ")
    x = x.replace("|", " | ")
    return x


def normalize_sample_index(df):
    """Make sample IDs comparable between microbiome and metadata tables."""
    df = df.copy()
    df.index = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip()
    return df


def read_csv_flexible(path, index_col=0):
    """
    Read a CSV file while handling common encoding differences.

    Some metadata files are exported from Excel/Windows and are not UTF-8.
    This function tries several common encodings before failing.
    """
    encodings_to_try = ["utf-8", "utf-8-sig", "cp1255", "cp1252", "latin1"]
    last_error = None

    for enc in encodings_to_try:
        try:
            df = pd.read_csv(path, index_col=index_col, encoding=enc)
            if enc != "utf-8":
                print(f"Loaded with encoding {enc}: {path}")
            return df
        except UnicodeDecodeError as e:
            last_error = e

    raise UnicodeDecodeError(
        last_error.encoding,
        last_error.object,
        last_error.start,
        last_error.end,
        f"Could not decode {path} with any of these encodings: {encodings_to_try}"
    )


def fdr_bh(p_values):
    """
    Apply Benjamini-Hochberg FDR correction.

    Returns adjusted p-values in the same order as the input.
    """
    p_values = np.asarray(p_values, dtype=float)
    q_values = np.full_like(p_values, np.nan, dtype=float)

    valid = ~np.isnan(p_values)
    p = p_values[valid]

    if len(p) == 0:
        return q_values

    n = len(p)
    order = np.argsort(p)
    ranked_p = p[order]

    q = ranked_p * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)

    q_unsorted = np.empty_like(q)
    q_unsorted[order] = q

    q_values[valid] = q_unsorted
    return q_values


def load_microbiome_table(path):
    """
    Load a microbiome table.

    Expected format:
    rows = samples
    columns = taxa / bacteria
    first column = sample IDs
    """
    df = read_csv_flexible(path, index_col=0)

    # Remove a taxonomy row if it exists by mistake.
    df = df[~df.index.astype(str).str.lower().isin(["taxonomy", "taxa", "taxon"])]

    # Make sure all abundance values are numeric.
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(axis=1, how="all")
    df = df.dropna(axis=0, how="all")
    df = df.fillna(0)

    # Remove taxa with total abundance 0.
    df = df.loc[:, df.sum(axis=0) > 0]

    # Normalize each sample to relative abundance.
    if NORMALIZE_MICROBIOME_ROWS:
        row_sums = df.sum(axis=1)
        df = df.div(row_sums.replace(0, np.nan), axis=0).fillna(0)

    df = normalize_sample_index(df)
    return df


def load_metadata(path):
    """
    Load a metadata table.

    Expected format:
    rows = samples
    columns = metadata variables
    first column = sample IDs
    """
    meta = read_csv_flexible(path, index_col=0)
    meta = normalize_sample_index(meta)
    return meta


def encode_metadata_for_correlations(meta):
    """
    Prepare metadata variables for correlation analysis.

    Numeric columns are kept as numeric variables.
    Categorical columns are converted into one-hot encoded variables.

    Returns:
    - encoded metadata dataframe
    - mapping table from encoded metadata feature to original metadata column
    """
    encoded_parts = []
    mapping_rows = []

    for col in meta.columns:
        if col in EXCLUDE_METADATA_COLUMNS:
            continue

        s = meta[col]

        # Skip columns with too few non-missing values.
        if s.notna().sum() < MIN_SAMPLES_PER_TEST:
            continue

        # Skip constant columns.
        if s.nunique(dropna=True) < 2:
            continue

        numeric_s = pd.to_numeric(s, errors="coerce")
        numeric_non_na = numeric_s.notna().sum()

        # Numeric metadata variable, including binary numeric variables.
        if numeric_non_na >= MIN_SAMPLES_PER_TEST:
            encoded_parts.append(pd.DataFrame({col: numeric_s}, index=meta.index))
            mapping_rows.append({
                "metadata_feature": col,
                "metadata_column": col,
                "metadata_type": "numeric",
            })
            continue

        # Categorical metadata variable.
        s_cat = s.astype("category")
        n_categories = s_cat.nunique(dropna=True)

        if n_categories < 2:
            continue

        # Avoid including free-text columns with too many categories.
        if n_categories > MAX_CATEGORIES_FOR_METADATA:
            print(f"Skipping metadata column '{col}' because it has {n_categories} categories")
            continue

        dummies = pd.get_dummies(s_cat, prefix=col, dummy_na=False)
        dummies = dummies.astype(float)

        # Drop dummy variables that are too rare or constant.
        keep_cols = []
        for dummy_col in dummies.columns:
            if dummies[dummy_col].nunique(dropna=True) < 2:
                continue
            if dummies[dummy_col].sum() < MIN_SAMPLES_PER_TEST:
                continue
            keep_cols.append(dummy_col)
            mapping_rows.append({
                "metadata_feature": dummy_col,
                "metadata_column": col,
                "metadata_type": "categorical_dummy",
            })

        if keep_cols:
            encoded_parts.append(dummies[keep_cols])

    if not encoded_parts:
        return pd.DataFrame(index=meta.index), pd.DataFrame(
            columns=["metadata_feature", "metadata_column", "metadata_type"]
        )

    encoded = pd.concat(encoded_parts, axis=1)
    mapping = pd.DataFrame(mapping_rows).drop_duplicates()
    return encoded, mapping


def compute_correlations_for_project(microbiome, metadata, project_name):
    """
    Compute Spearman correlations between each taxon and each metadata feature
    for one project/dataset.
    """
    common_samples = microbiome.index.intersection(metadata.index)

    microbiome = microbiome.loc[common_samples].copy()
    metadata = metadata.loc[common_samples].copy()

    print("\n" + "=" * 80)
    print(project_name)
    print(f"Common samples: {len(common_samples)}")
    print(f"Microbiome shape: {microbiome.shape}")
    print(f"Metadata shape: {metadata.shape}")

    if len(common_samples) < MIN_SAMPLES_PER_TEST:
        print(f"Skipping {project_name}: not enough common samples")
        return pd.DataFrame()

    encoded_meta, metadata_mapping = encode_metadata_for_correlations(metadata)

    print(f"Usable metadata features after encoding: {encoded_meta.shape[1]}")

    if encoded_meta.shape[1] == 0:
        print(f"No usable metadata columns for {project_name}")
        return pd.DataFrame()

    feature_to_original_col = dict(
        zip(metadata_mapping["metadata_feature"], metadata_mapping["metadata_column"])
    )
    feature_to_type = dict(
        zip(metadata_mapping["metadata_feature"], metadata_mapping["metadata_type"])
    )

    results = []

    for metadata_feature in encoded_meta.columns:
        meta_values = encoded_meta[metadata_feature]

        # Skip constant metadata features.
        if meta_values.nunique(dropna=True) < 2:
            continue

        for taxon in microbiome.columns:
            taxon_values = microbiome[taxon]

            # Skip taxa that appear in too few samples.
            if (taxon_values > 0).sum() < MIN_NONZERO_SAMPLES_PER_TAXON:
                continue

            valid = meta_values.notna() & taxon_values.notna()

            if valid.sum() < MIN_SAMPLES_PER_TEST:
                continue

            x = meta_values[valid]
            y = taxon_values[valid]

            if x.nunique() < 2 or y.nunique() < 2:
                continue

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                corr, pval = spearmanr(x, y)

            if np.isnan(pval):
                continue

            results.append({
                "project": project_name,
                "metadata_feature": metadata_feature,
                "metadata_column": feature_to_original_col.get(metadata_feature, metadata_feature),
                "metadata_type": feature_to_type.get(metadata_feature, "unknown"),
                "bacteria_name": taxon,
                "correlation_coefficient": corr,
                "p_value": pval,
                "n_samples": int(valid.sum()),
            })

    res = pd.DataFrame(results)

    if len(res) == 0:
        print(f"No correlations calculated for {project_name}")
        return res

    # Apply FDR correction within the current project only.
    res["p_value_fdr"] = fdr_bh(res["p_value"].values)

    n_sig = (res["p_value_fdr"] < P_VALUE_THRESHOLD).sum()
    print(f"Total correlations calculated: {len(res)}")
    print(f"Significant correlations, FDR < {P_VALUE_THRESHOLD}: {n_sig}")

    return res


def make_significant_counts_table(corr_df, project_name, output_dir):
    """
    Count the number of significant correlations per metadata feature/column.
    """
    if corr_df.empty:
        return pd.DataFrame()

    if PLOT_BY not in {"metadata_feature", "metadata_column"}:
        raise ValueError("PLOT_BY must be either 'metadata_feature' or 'metadata_column'")

    sig = corr_df[corr_df["p_value_fdr"] < P_VALUE_THRESHOLD].copy()

    if sig.empty:
        counts = pd.DataFrame(columns=["project", PLOT_BY, "count"])
        counts_file = os.path.join(output_dir, "significant_correlation_counts.csv")
        counts.to_csv(counts_file, index=False)
        print(f"Saved empty counts table: {counts_file}")
        return counts

    counts = (
        sig
        .groupby(["project", PLOT_BY])
        .size()
        .reset_index(name="count")
        .sort_values("count", ascending=False)
    )

    counts_file = os.path.join(output_dir, "significant_correlation_counts.csv")
    counts.to_csv(counts_file, index=False)
    print(f"Saved counts table: {counts_file}")

    return counts


def plot_number_of_significant_correlations_for_project(corr_df, project_name, output_dir):
    """
    Create one horizontal bar plot for one project.

    The plot shows how many significant correlations were found between each
    metadata feature and all microbial taxa in this project.
    """
    os.makedirs(output_dir, exist_ok=True)

    counts = make_significant_counts_table(corr_df, project_name, output_dir)

    if counts.empty:
        fig, ax = plt.subplots(figsize=(10, 4), dpi=150)
        ax.text(
            0.5,
            0.5,
            "No significant correlations found\nFDR-adjusted p < "
            f"{P_VALUE_THRESHOLD}",
            ha="center",
            va="center",
            fontsize=12,
            transform=ax.transAxes,
        )
        ax.set_title(project_name)
        ax.axis("off")

        out_png = os.path.join(output_dir, "number_of_significant_correlations.png")
        out_pdf = os.path.join(output_dir, "number_of_significant_correlations.pdf")
        plt.savefig(out_png, bbox_inches="tight")
        plt.savefig(out_pdf, bbox_inches="tight")
        plt.close(fig)

        print(f"No significant correlations found for {project_name}. Saved empty plot:")
        print(out_png)
        print(out_pdf)
        return

    counts = counts.sort_values("count", ascending=True).tail(TOP_N_METADATA_TO_PLOT)

    fig_height = max(6, 0.35 * len(counts))
    fig, ax = plt.subplots(figsize=(12, fig_height), dpi=150)

    y_pos = np.arange(len(counts))
    values = counts["count"].values
    labels = counts[PLOT_BY].values

    ax.barh(y_pos, values)

    # Add the count at the end of each bar.
    max_value = max(values) if len(values) else 0
    label_offset = max(0.5, max_value * 0.01)
    for i, value in enumerate(values):
        ax.text(value + label_offset, i, str(int(value)), va="center", fontsize=9)

    ax.set_yticks(y_pos)
    ax.set_yticklabels([clean_label(x) for x in labels], fontsize=9)

    ax.set_xlabel("Number of significant correlations")
    ax.set_ylabel("Metadata feature" if PLOT_BY == "metadata_feature" else "Metadata column")
    ax.set_title(
        f"{project_name}\n"
        f"Significant microbe-metadata correlations, FDR-adjusted p < {P_VALUE_THRESHOLD}"
    )

    ax.set_xlim(0, max_value + label_offset * 8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()

    out_png = os.path.join(output_dir, "number_of_significant_correlations.png")
    out_pdf = os.path.join(output_dir, "number_of_significant_correlations.pdf")

    plt.savefig(out_png, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved plot: {out_png}")
    print(f"Saved plot: {out_pdf}")


def run_one_project(project):
    """Run the full analysis for one project."""
    project_name = project["project_name"]
    microbiome_file = str(project.get("microbiome_file", "")).strip()
    metadata_file = str(project.get("metadata_file", "")).strip()

    project_output_dir = os.path.join(OUTPUT_DIR, safe_name(project_name))
    os.makedirs(project_output_dir, exist_ok=True)

    if not microbiome_file:
        print(f"Skipping {project_name}: no microbiome file was provided")
        return pd.DataFrame()

    if not os.path.exists(microbiome_file):
        print(f"Skipping {project_name}: microbiome file not found:")
        print(microbiome_file)
        return pd.DataFrame()

    if not metadata_file:
        print(f"Skipping {project_name}: no metadata file was provided")
        return pd.DataFrame()

    if not os.path.exists(metadata_file):
        print(f"Skipping {project_name}: metadata file not found:")
        print(metadata_file)
        return pd.DataFrame()

    print(f"Loading project: {project_name}")
    print(f"Microbiome file: {microbiome_file}")
    print(f"Metadata file: {metadata_file}")

    microbiome = load_microbiome_table(microbiome_file)
    metadata = load_metadata(metadata_file)

    corr_df = compute_correlations_for_project(
        microbiome=microbiome,
        metadata=metadata,
        project_name=project_name,
    )

    if corr_df.empty:
        return corr_df

    corr_file = os.path.join(project_output_dir, "all_microbe_metadata_correlations.csv")
    corr_df.to_csv(corr_file, index=False)
    print(f"Saved full correlations table: {corr_file}")

    plot_number_of_significant_correlations_for_project(
        corr_df=corr_df,
        project_name=project_name,
        output_dir=project_output_dir,
    )

    return corr_df


# ============================================================
# MAIN
# ============================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    all_results = []

    for project in PROJECTS:
        corr_df = run_one_project(project)
        if not corr_df.empty:
            all_results.append(corr_df)

    if all_results:
        all_corr = pd.concat(all_results, ignore_index=True)
        all_corr_file = os.path.join(OUTPUT_DIR, "all_projects_microbe_metadata_correlations.csv")
        all_corr.to_csv(all_corr_file, index=False)
        print("\nSaved combined correlations table:")
        print(all_corr_file)

        if "p_value_fdr" in all_corr.columns:
            if PLOT_BY not in {"metadata_feature", "metadata_column"}:
                raise ValueError("PLOT_BY must be either 'metadata_feature' or 'metadata_column'")

            all_sig = all_corr[all_corr["p_value_fdr"] < P_VALUE_THRESHOLD].copy()
            all_counts = (
                all_sig
                .groupby(["project", PLOT_BY])
                .size()
                .reset_index(name="count")
                .sort_values(["project", "count"], ascending=[True, False])
            )
            all_counts_file = os.path.join(OUTPUT_DIR, "all_projects_significant_correlation_counts.csv")
            all_counts.to_csv(all_counts_file, index=False)
            print("\nSaved combined significant counts table:")
            print(all_counts_file)
    else:
        print("\nNo correlation results were created for any project.")


if __name__ == "__main__":
    main()