'''
THIS IS AN INITIAL DRAFT
NEED TO MAKE CHANGES
JUST A BASE FOR NOW
'''

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr


# ============================================================
# INPUTS
# ============================================================
# Rows = samples, columns = taxa / bacteria
PREG_MICROBIOME_FILE = "pregnancy_microbiome.csv"
CTRL_MICROBIOME_FILE = "control_microbiome.csv"

# Rows = samples, columns = metadata variables
PREG_METADATA_FILE = "pregnancy_metadata.csv"
CTRL_METADATA_FILE = "control_metadata.csv"

OUTPUT_DIR = "pregnancy_vs_control_significant_correlations"

# Analysis settings
NORMALIZE_MICROBIOME_ROWS = True
P_VALUE_THRESHOLD = 0.05
MIN_SAMPLES_PER_TEST = 10
MIN_NONZERO_SAMPLES_PER_TAXON = 3
MAX_CATEGORIES_FOR_METADATA = 10
TOP_N_METADATA_TO_PLOT = 20

# Metadata columns to exclude from correlation analysis
EXCLUDE_METADATA_COLUMNS = {
    "sample_id", "SampleID", "sample", "Sample",
    "subject_id", "SubjectID", "participant_id",
    "run", "Run", "project", "Project",
}


# ============================================================
# FUNCTIONS
# ============================================================

def clean_label(x):
    """Create a cleaner label for plotting."""
    x = str(x)
    x = x.replace("_", " ")
    return x


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
    """
    df = pd.read_csv(path, index_col=0)

    # Remove a taxonomy row if it exists by mistake
    df = df[~df.index.astype(str).str.lower().isin(["taxonomy", "taxa", "taxon"])]

    # Make sure all values are numeric
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(axis=1, how="all")
    df = df.dropna(axis=0, how="all")
    df = df.fillna(0)

    # Normalize each sample to relative abundance
    if NORMALIZE_MICROBIOME_ROWS:
        row_sums = df.sum(axis=1)
        df = df.div(row_sums.replace(0, np.nan), axis=0).fillna(0)

    df.index = df.index.astype(str)
    return df


def load_metadata(path):
    """
    Load a metadata table.

    Expected format:
    rows = samples
    columns = metadata variables
    """
    meta = pd.read_csv(path, index_col=0)
    meta.index = meta.index.astype(str)
    return meta


def encode_metadata_for_correlations(meta):
    """
    Prepare metadata variables for correlation analysis.

    - Numeric columns are kept as numeric variables.
    - Categorical columns are converted into one-hot encoded variables.
    """
    encoded_parts = []

    for col in meta.columns:
        if col in EXCLUDE_METADATA_COLUMNS:
            continue

        s = meta[col]

        # Skip columns with too few non-missing values
        if s.notna().sum() < MIN_SAMPLES_PER_TEST:
            continue

        numeric_s = pd.to_numeric(s, errors="coerce")
        numeric_non_na = numeric_s.notna().sum()

        # Numeric metadata variable
        if numeric_non_na >= MIN_SAMPLES_PER_TEST and s.nunique(dropna=True) > 2:
            encoded_parts.append(pd.DataFrame({col: numeric_s}, index=meta.index))
            continue

        # Binary numeric metadata variable
        if numeric_non_na >= MIN_SAMPLES_PER_TEST and s.nunique(dropna=True) == 2:
            encoded_parts.append(pd.DataFrame({col: numeric_s}, index=meta.index))
            continue

        # Categorical metadata variable
        s_cat = s.astype("category")
        n_categories = s_cat.nunique(dropna=True)

        if n_categories < 2:
            continue

        # Avoid including free-text columns with too many categories
        if n_categories > MAX_CATEGORIES_FOR_METADATA:
            print(f"Skipping metadata column '{col}' because it has {n_categories} categories")
            continue

        dummies = pd.get_dummies(s_cat, prefix=col, dummy_na=False)
        dummies = dummies.astype(float)
        encoded_parts.append(dummies)

    if not encoded_parts:
        return pd.DataFrame(index=meta.index)

    encoded = pd.concat(encoded_parts, axis=1)
    return encoded


def compute_correlations_for_group(microbiome, metadata, group_name):
    """
    Compute Spearman correlations between each taxon and each metadata variable
    for one group: Pregnancy or Control.
    """
    common_samples = microbiome.index.intersection(metadata.index)

    microbiome = microbiome.loc[common_samples].copy()
    metadata = metadata.loc[common_samples].copy()

    print("\n" + "=" * 80)
    print(group_name)
    print(f"Common samples: {len(common_samples)}")
    print(f"Microbiome shape: {microbiome.shape}")
    print(f"Metadata shape: {metadata.shape}")

    if len(common_samples) < MIN_SAMPLES_PER_TEST:
        raise ValueError(f"Not enough common samples for {group_name}")

    encoded_meta = encode_metadata_for_correlations(metadata)

    print(f"Usable metadata variables after encoding: {encoded_meta.shape[1]}")

    if encoded_meta.shape[1] == 0:
        print(f"No usable metadata columns for {group_name}")
        return pd.DataFrame()

    results = []

    for metadata_col in encoded_meta.columns:
        meta_values = encoded_meta[metadata_col]

        # Skip constant metadata variables
        if meta_values.nunique(dropna=True) < 2:
            continue

        for taxon in microbiome.columns:
            taxon_values = microbiome[taxon]

            # Skip taxa that appear in too few samples
            if (taxon_values > 0).sum() < MIN_NONZERO_SAMPLES_PER_TAXON:
                continue

            valid = meta_values.notna() & taxon_values.notna()

            if valid.sum() < MIN_SAMPLES_PER_TEST:
                continue

            x = meta_values[valid]
            y = taxon_values[valid]

            if x.nunique() < 2 or y.nunique() < 2:
                continue

            corr, pval = spearmanr(x, y)

            results.append({
                "group": group_name,
                "metadata_column": metadata_col,
                "bacteria_name": taxon,
                "correlation_coefficient": corr,
                "p_value": pval,
                "n_samples": int(valid.sum())
            })

    res = pd.DataFrame(results)

    if len(res) == 0:
        print(f"No correlations calculated for {group_name}")
        return res

    # Apply FDR correction within each group
    res["p_value_fdr"] = fdr_bh(res["p_value"].values)

    print(f"Total correlations calculated: {len(res)}")
    print(f"Significant correlations, FDR < {P_VALUE_THRESHOLD}: {(res['p_value_fdr'] < P_VALUE_THRESHOLD).sum()}")

    return res


def plot_number_of_significant_correlations(corr_df, output_dir):
    """
    Create a stacked horizontal bar plot.

    For each metadata variable, the plot shows how many significant correlations
    were found in the Pregnancy group and how many were found in the Control group.
    """
    os.makedirs(output_dir, exist_ok=True)

    sig = corr_df[corr_df["p_value_fdr"] < P_VALUE_THRESHOLD].copy()

    if len(sig) == 0:
        print("No significant correlations found.")
        return

    counts = (
        sig
        .groupby(["metadata_column", "group"])
        .size()
        .reset_index(name="count")
    )

    counts_file = os.path.join(output_dir, "significant_correlation_counts_by_group.csv")
    counts.to_csv(counts_file, index=False)

    print(f"Saved counts table: {counts_file}")

    pivot = counts.pivot_table(
        index="metadata_column",
        columns="group",
        values="count",
        fill_value=0
    )

    # Keep a consistent group order
    preferred_order = ["Pregnancy", "Control"]
    existing = list(pivot.columns)
    ordered_groups = [g for g in preferred_order if g in existing]
    ordered_groups += [g for g in existing if g not in ordered_groups]
    pivot = pivot[ordered_groups]

    # Keep only the top metadata variables according to the total number of correlations
    pivot["total"] = pivot.sum(axis=1)
    pivot = pivot.sort_values("total", ascending=True).tail(TOP_N_METADATA_TO_PLOT)

    totals = pivot["total"].copy()
    pivot = pivot.drop(columns=["total"])

    fig_height = max(6, 0.35 * len(pivot))
    fig, ax = plt.subplots(figsize=(14, fig_height), dpi=150)

    y_pos = np.arange(len(pivot))
    left = np.zeros(len(pivot))

    colors = {
        "Pregnancy": "#E76F51",
        "Control": "#457B9D",
    }

    for group in pivot.columns:
        values = pivot[group].values

        ax.barh(
            y_pos,
            values,
            left=left,
            label=group,
            color=colors.get(group, None)
        )

        left += values

    # Add the total number of significant correlations at the end of each bar
    for i, total in enumerate(totals.values):
        ax.text(total + 0.5, i, str(int(total)), va="center", fontsize=9)

    ax.set_yticks(y_pos)
    ax.set_yticklabels([clean_label(x) for x in pivot.index], fontsize=9)

    ax.set_xlabel("Number of Significant Correlations")
    ax.set_title("Significant microbe-metadata correlations: Pregnancy vs Control")

    ax.legend(title="Group", bbox_to_anchor=(1.02, 1), loc="upper left")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()

    out_png = os.path.join(output_dir, "pregnancy_vs_control_number_of_significant_correlations.png")
    out_pdf = os.path.join(output_dir, "pregnancy_vs_control_number_of_significant_correlations.pdf")

    plt.savefig(out_png, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.show()

    print(f"Saved plot: {out_png}")
    print(f"Saved plot: {out_pdf}")


# ============================================================
# MAIN
# ============================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    preg_micro = load_microbiome_table(PREG_MICROBIOME_FILE)
    ctrl_micro = load_microbiome_table(CTRL_MICROBIOME_FILE)

    preg_meta = load_metadata(PREG_METADATA_FILE)
    ctrl_meta = load_metadata(CTRL_METADATA_FILE)

    preg_corr = compute_correlations_for_group(
        preg_micro,
        preg_meta,
        group_name="Pregnancy"
    )

    ctrl_corr = compute_correlations_for_group(
        ctrl_micro,
        ctrl_meta,
        group_name="Control"
    )

    all_corr = pd.concat([preg_corr, ctrl_corr], ignore_index=True)

    all_corr_file = os.path.join(OUTPUT_DIR, "all_microbe_metadata_correlations_by_group.csv")
    all_corr.to_csv(all_corr_file, index=False)

    print("\nSaved full correlations table:")
    print(all_corr_file)

    plot_number_of_significant_correlations(all_corr, OUTPUT_DIR)


if __name__ == "__main__":
    main()