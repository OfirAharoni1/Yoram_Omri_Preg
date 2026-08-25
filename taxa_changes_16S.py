import pandas as pd
import numpy as np
import re
from scipy.stats import fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from pathlib import Path

# =========================
# Input files
# rows = samples, columns = taxa
# =========================

PREGNANT_FILES = [
    "/home/aharonox/Yoram_Omri_Preg/datasets_after_MIPMLP/omri/omri_stool_Pregnant_16S.csv",
    "/home/aharonox/Yoram_Omri_Preg/datasets_after_MIPMLP/PRJNA1254708_Pregnant_16S.csv",
]

CONTROL_FILES = [
    "/home/aharonox/Yoram_Omri_Preg/datasets_after_MIPMLP/PRJNA669650_Control_16S.csv",
    "/home/aharonox/Yoram_Omri_Preg/datasets_after_MIPMLP/PRJNA388263_Control_16S.csv",
    "/home/aharonox/Yoram_Omri_Preg/datasets_after_MIPMLP/PRJNA545251_Control_16S.csv",
]

OUTDIR = Path("/home/aharonox/Yoram_Omri_Preg/consistent_pregnancy_taxa_16S")
OUTDIR.mkdir(exist_ok=True)

# =========================
# Settings
# =========================

NORMALIZE_ROWS = True
PRESENT_THR = 0

# Genus level:
# Bacteria;Phylum;Class;Order;Family;Genus
TAXONOMY_LEVEL = 6

ALPHA = 0.05
PREV_HIGH = 0.20
PREV_LOW = 0.02
MIN_LOG2FC = 0.5
CONSISTENCY = 0.8


def clean_taxon_name(taxon):
    taxon = str(taxon).strip()
    taxon = taxon.replace("|", ";").replace(",", ";")

    prefixes = [
        "k__", "p__", "c__", "o__", "f__", "g__", "s__",
        "d__", "D_0__", "D_1__", "D_2__", "D_3__", "D_4__", "D_5__", "D_6__"
    ]

    parts = []
    for part in taxon.split(";"):
        part = part.strip()

        for prefix in prefixes:
            if part.startswith(prefix):
                part = part.replace(prefix, "", 1)

        # Defensive cleanup for outputs created before preprocess_16S.py
        # standardized per-rank confidence annotations.
        part = re.sub(
            r"\s*\((?:0(?:\.\d+)?|1(?:\.0+)?)\)\s*$",
            "",
            part,
        ).strip()

        if part and part.lower() not in {"nan", "none"}:
            parts.append(part)

    return ";".join(parts)


def keep_taxonomy_level(taxon, level=TAXONOMY_LEVEL):
    parts = [p.strip() for p in str(taxon).split(";") if p.strip()]

    # Keep only taxa that reach genus level
    if len(parts) < level:
        return None

    return ";".join(parts[:level])


def remove_uninformative_taxa(df):
    bad_exact = {
        "",
        "Unassigned",
        "unassigned",
        "Bacteria",
        "bacteria",
        "Unknown",
        "unknown",
    }

    bad_words = {
        "unassigned",
        "unknown",
        "uncultured",
        "unclassified",
        "metagenome",
    }

    keep_cols = []

    for col in df.columns:
        col_str = str(col).strip()
        parts = [p.strip().lower() for p in col_str.split(";") if p.strip()]
        last_level = parts[-1] if parts else ""

        if col_str in bad_exact:
            continue

        if last_level in {"", "bacteria"}:
            continue

        if any(bad in last_level for bad in bad_words):
            continue

        keep_cols.append(col)

    return df[keep_cols]


def load_project(path, group):
    df = pd.read_csv(path, index_col=0)

    # Remove taxonomy row if it exists
    df = df[~df.index.astype(str).str.lower().isin(["taxonomy", "taxa", "taxon"])]

    # Clean taxon names and keep only genus-level taxa
    new_cols = [
        keep_taxonomy_level(clean_taxon_name(c))
        for c in df.columns
    ]

    df.columns = new_cols

    # Drop taxa that did not reach genus level
    df = df.loc[:, df.columns.notna()]

    # Convert values to numeric
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)

    # Merge duplicated taxa created after cleaning
    df = df.T.groupby(level=0).sum().T

    # Remove general / uninformative taxa
    df = remove_uninformative_taxa(df)

    # Remove empty taxa
    df = df.loc[:, df.sum(axis=0) > 0]

    # Normalize to relative abundance
    if NORMALIZE_ROWS:
        df = df.div(df.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)

    project = Path(path).stem

    # Make sample IDs unique across projects
    df.index = [f"{project}__{sample}" for sample in df.index]

    meta = pd.DataFrame({
        "group": group,
        "project": project
    }, index=df.index)

    present_taxa = set(df.columns[(df > PRESENT_THR).any(axis=0)])

    print(group, project, "samples:", df.shape[0], "taxa:", df.shape[1])

    return df, meta, present_taxa


# =========================
# Load data
# =========================

preg_data = [load_project(f, "Pregnant") for f in PREGNANT_FILES]
ctrl_data = [load_project(f, "Control") for f in CONTROL_FILES]

preg_taxa_all_projects = set.intersection(*[x[2] for x in preg_data])
ctrl_taxa_all_projects = set.intersection(*[x[2] for x in ctrl_data])

taxa_to_test = sorted(preg_taxa_all_projects | ctrl_taxa_all_projects)

print("\nTaxa in all pregnancy projects:", len(preg_taxa_all_projects))
print("Taxa in all control projects:", len(ctrl_taxa_all_projects))
print("Total taxa to test:", len(taxa_to_test))

if len(taxa_to_test) == 0:
    print("\nNo taxa to test.")
    raise SystemExit


dfs = []
metas = []

for df, meta, _ in preg_data + ctrl_data:
    dfs.append(df.reindex(columns=taxa_to_test, fill_value=0))
    metas.append(meta)

data = pd.concat(dfs, axis=0)
meta = pd.concat(metas, axis=0)

# =========================
# Tests
# =========================

results = []

for taxon in taxa_to_test:
    x = data[taxon]

    preg = x[meta["group"] == "Pregnant"]
    ctrl = x[meta["group"] == "Control"]

    preg_present = (preg > PRESENT_THR).sum()
    ctrl_present = (ctrl > PRESENT_THR).sum()

    preg_prev = preg_present / len(preg)
    ctrl_prev = ctrl_present / len(ctrl)

    fisher_p = fisher_exact([
        [preg_present, len(preg) - preg_present],
        [ctrl_present, len(ctrl) - ctrl_present]
    ])[1]

    mw_p = mannwhitneyu(preg, ctrl, alternative="two-sided").pvalue

    preg_median = preg.median()
    ctrl_median = ctrl.median()

    log2fc = np.log2((preg_median + 1e-9) / (ctrl_median + 1e-9))

    preg_proj_medians = x[meta["group"] == "Pregnant"].groupby(
        meta.loc[meta["group"] == "Pregnant", "project"]
    ).median()

    ctrl_proj_medians = x[meta["group"] == "Control"].groupby(
        meta.loc[meta["group"] == "Control", "project"]
    ).median()

    pairwise_diffs = [
        p - c
        for p in preg_proj_medians
        for c in ctrl_proj_medians
    ]

    support_higher = np.mean([d > 0 for d in pairwise_diffs])
    support_lower = np.mean([d < 0 for d in pairwise_diffs])

    results.append({
        "taxon": taxon,

        "in_all_pregnancy_projects": taxon in preg_taxa_all_projects,
        "in_all_control_projects": taxon in ctrl_taxa_all_projects,

        "pregnancy_prevalence": preg_prev,
        "control_prevalence": ctrl_prev,
        "prevalence_diff_preg_minus_control": preg_prev - ctrl_prev,

        "pregnancy_median": preg_median,
        "control_median": ctrl_median,
        "log2FC_preg_vs_control": log2fc,

        "fisher_p": fisher_p,
        "mannwhitney_p": mw_p,

        "support_higher_in_pregnancy": support_higher,
        "support_lower_in_pregnancy": support_lower,
    })

res = pd.DataFrame(results)

if res.empty:
    print("\nNo results were created.")
    raise SystemExit

# =========================
# FDR correction
# =========================

res["fisher_q"] = multipletests(res["fisher_p"], method="fdr_bh")[1]
res["mannwhitney_q"] = multipletests(res["mannwhitney_p"], method="fdr_bh")[1]

# =========================
# Classification
# =========================

res["category"] = "not_candidate"

res.loc[
    (res["fisher_q"] < ALPHA) &
    (res["pregnancy_prevalence"] >= PREV_HIGH) &
    (res["control_prevalence"] <= PREV_LOW),
    "category"
] = "appears_in_pregnancy"

res.loc[
    (res["fisher_q"] < ALPHA) &
    (res["control_prevalence"] >= PREV_HIGH) &
    (res["pregnancy_prevalence"] <= PREV_LOW),
    "category"
] = "disappears_in_pregnancy"

res.loc[
    (res["category"] == "not_candidate") &
    (res["mannwhitney_q"] < ALPHA) &
    (res["log2FC_preg_vs_control"] >= MIN_LOG2FC) &
    (res["support_higher_in_pregnancy"] >= CONSISTENCY),
    "category"
] = "higher_in_pregnancy"

res.loc[
    (res["category"] == "not_candidate") &
    (res["mannwhitney_q"] < ALPHA) &
    (res["log2FC_preg_vs_control"] <= -MIN_LOG2FC) &
    (res["support_lower_in_pregnancy"] >= CONSISTENCY),
    "category"
] = "lower_in_pregnancy"

res = res.sort_values([
    "category",
    "fisher_q",
    "mannwhitney_q",
    "log2FC_preg_vs_control"
])

# =========================
# Save
# =========================

res.to_csv(OUTDIR / "all_tested_taxa_results.csv", index=False)
res[res["category"] != "not_candidate"].to_csv(OUTDIR / "candidate_taxa.csv", index=False)

print("\nDone")
print("Tested taxa:", len(res))
print(res["category"].value_counts())
print("Saved to:", OUTDIR)
