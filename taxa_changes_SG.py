import pandas as pd
import numpy as np
import re
from scipy.stats import fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from pathlib import Path

# =========================
# Input files
# rows = samples, columns = taxa
# Pregnancy cohorts: ERP020710 and PRJNA1247940.
# =========================

PREGNANT_FILES = [
    "/home/aharonox/Yoram_Omri_Preg/datasets_after_yamas/SG/ERP020710_Pregnant_SG_formatted.csv",
    "/home/aharonox/Yoram_Omri_Preg/datasets_after_yamas/SG/PRJNA1247940_Pregnant_SG_formatted.csv",
]

CONTROL_FILES = [
    "/home/aharonox/Yoram_Omri_Preg/datasets_after_yamas/SG/PRJEB37731_Control_SG_formatted.csv",
    "/home/aharonox/Yoram_Omri_Preg/datasets_after_yamas/SG/PRJNA48479_Control_SG_formatted.csv",
]

OUTDIR = Path("/home/aharonox/Yoram_Omri_Preg/consistent_pregnancy_taxa_SG")
OUTDIR.mkdir(exist_ok=True)

# =========================
# Settings
# =========================

NORMALIZE_ROWS = True
PRESENT_THR = 0

# Genus level:
# Bacteria;Phylum;Class;Order;Family;Genus
# Intermediate ranks may be missing; the genus itself must be informative.
TAXONOMY_LEVEL = 6

ALPHA = 0.05
PREV_HIGH = 0.20
PREV_LOW = 0.02
MIN_LOG2FC = 0.5
CONSISTENCY = 0.8


def clean_taxon_name(taxon):
    """Preserve rank positions, including missing ranks, when cleaning names."""
    ranks = ("k", "p", "c", "o", "f", "g", "s")
    ranked = {}
    plain = []
    for part in re.split(r"[;,|]", str(taxon).strip()):
        part = part.strip()
        match = re.match(r"^(?:D_([0-6])__|([dkpcofgs])(?:__|:))\s*(.*)$",
                         part, flags=re.IGNORECASE)
        name = match.group(3) if match else part
        name = re.sub(r"\s*\((?:0(?:\.\d+)?|1(?:\.0+)?)\)\s*$", "", name).strip()
        if name.lower() in {"nan", "none"}:
            name = ""
        if match:
            rank = ranks[int(match.group(1))] if match.group(1) is not None else match.group(2).lower()
            ranked["k" if rank == "d" else rank] = name
        else:
            plain.append(name)
    return ";".join(ranked.get(rank, "") for rank in ranks) if ranked else ";".join(plain)


def keep_taxonomy_level(taxon, level=TAXONOMY_LEVEL):
    """Match exact genus names; retain database suffixes such as _A and _B."""
    parts = [p.strip() for p in str(taxon).split(";")]
    if len(parts) < level or not parts[level - 1]:
        return None
    return parts[level - 1]


PROJECT_SUMMARY = []
EXCLUDED_SAMPLES = []


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

    input_samples = len(df)

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

    # A sample without classified genus abundance is not biological absence.
    no_genus = df.sum(axis=1) <= 0
    project = Path(path).stem
    EXCLUDED_SAMPLES.extend(
        {"project": project, "group": group, "sample_id": sample,
         "reason": "no_informative_genus_abundance"}
        for sample in df.index[no_genus]
    )
    df = df.loc[~no_genus].copy()
    PROJECT_SUMMARY.append({
        "project": project, "group": group, "input_samples": input_samples,
        "excluded_samples": int(no_genus.sum()), "retained_samples": len(df),
        "informative_genera": len(df.columns),
    })

    # Normalize over all informative genera, before restricting to matched genera.
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

print("Input projects:")
for group, paths in (("Pregnant", PREGNANT_FILES), ("Control", CONTROL_FILES)):
    for path in paths:
        print(f"  {group}: {path}")

preg_data = [load_project(f, "Pregnant") for f in PREGNANT_FILES]
ctrl_data = [load_project(f, "Control") for f in CONTROL_FILES]

preg_taxa_all_projects = set.intersection(*[x[2] for x in preg_data])
ctrl_taxa_all_projects = set.intersection(*[x[2] for x in ctrl_data])

# Require an observed exact genus match in every project. Missing columns
# cannot distinguish biological absence from incompatible taxonomy labels.
taxa_to_test = sorted(preg_taxa_all_projects & ctrl_taxa_all_projects)
all_projects = preg_data + ctrl_data
all_genera = sorted(set.union(*[item[2] for item in all_projects]))
coverage = pd.DataFrame({
    summary["project"]: [genus in item[2] for genus in all_genera]
    for summary, item in zip(PROJECT_SUMMARY, all_projects)
}, index=pd.Index(all_genera, name="genus"))
coverage["included_in_tests"] = coverage.index.isin(taxa_to_test)
coverage.to_csv(OUTDIR / "genus_matching_summary.csv")
for summary in PROJECT_SUMMARY:
    summary["matched_genera"] = len(taxa_to_test)
    summary["unmatched_genera"] = summary["informative_genera"] - len(taxa_to_test)
pd.DataFrame(PROJECT_SUMMARY).to_csv(OUTDIR / "project_summary.csv", index=False)
pd.DataFrame(EXCLUDED_SAMPLES, columns=["project", "group", "sample_id", "reason"]).to_csv(
    OUTDIR / "excluded_samples.csv", index=False)
(OUTDIR / "analysis_notes.txt").write_text(
    "Input files for this run:\n"
    + "".join(f"Pregnant: {path}\n" for path in PREGNANT_FILES)
    + "".join(f"Control: {path}\n" for path in CONTROL_FILES)
    + "\n"
    "Exploratory comparison using exact genus names shared and observed in every project.\n"
    "Higher taxonomy ranks are ignored; genus suffixes are preserved.\n"
    "Unmatched genera are excluded, not filled with biological zeros.\n"
    "Samples with no informative genus abundance are excluded.\n"
    "Abundances are normalized over all informative genera before matching.\n"
    "Taxonomy and project differences remain potential confounders.\n"
    "Appears/disappears categories describe pooled prevalence thresholds among matched genera, "
    "not absolute presence/absence or a causal pregnancy effect.\n"
)
print(pd.DataFrame(PROJECT_SUMMARY).to_string(index=False))

print("\nTaxa in all pregnancy projects:", len(preg_taxa_all_projects))
print("Taxa in all control projects:", len(ctrl_taxa_all_projects))
print("Total taxa to test:", len(taxa_to_test))

if len(taxa_to_test) == 0:
    print("\nNo taxa to test.")
    raise SystemExit


dfs = []
metas = []

for df, meta, _ in preg_data + ctrl_data:
    dfs.append(df.loc[:, taxa_to_test])
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
print("\nTaxa by category:")
for category in (
    "appears_in_pregnancy",
    "disappears_in_pregnancy",
    "higher_in_pregnancy",
    "lower_in_pregnancy",
    "not_candidate",
):
    taxa = sorted(res.loc[res["category"] == category, "taxon"])
    print(f"\n{category} ({len(taxa)}):")
    if taxa:
        for taxon in taxa:
            print(f"  - {taxon}")
    else:
        print("  (none)")

print("Saved to:", OUTDIR)
