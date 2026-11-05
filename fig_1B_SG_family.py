
from pathlib import Path
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
import numpy as np

# ===== SETTINGS =====
BASE_DIR = Path(__file__).resolve().parent
input_folder = BASE_DIR / "datasets_after_yamas/SG"
output_folder = BASE_DIR / "1B_plots/SG_family"

sample_id_col = "SampleID"
group_col = "Group"
top_n = 12   # keep top families per dataset, merge the rest into "Other" in the stacked plots


def extract_family(taxon_name):
    """Return the full lineage only for an explicitly family-level lineage.

    Accept comma, pipe or semicolon separators; never merge taxonomic synonyms.
    Kingdom totals and descendant ranks are deliberately excluded for SG tables.
    """
    import re
    parts = [p.strip() for p in re.split(r"[,;|]", str(taxon_name))]
    if len(parts) != 5:
        return None
    if not parts[0].startswith(("k__", "d__")) or not all(
        part.startswith(prefix)
        for part, prefix in zip(parts[1:], ("p__", "c__", "o__", "f__"))
    ):
        return None
    names = [part[3:].strip() for part in parts]
    if not all(names):
        return None
    return ";".join(names)


def prepare_relative_family_table(df, sample_id_col, group_col, file_name=""):
    """Reference-equivalent level-5 mean grouping, no transform, then row %.

    Implements the explicit reference operations locally (no MIPMLP dependency):
    taxonomy_level=5, taxnomy_group='mean', normalization='none', then closure.
    This is NOT a claim of equivalence to version-dependent MIPMLP defaults
    such as rare-feature filtering. No implicit filtering, log, z-score or PCA.
    SG hierarchical data must first be restricted to explicit family columns;
    averaging parents and descendants would reweight families by feature counts.
    Duplicate canonical family profiles must agree; conflicting profiles fail
    rather than being summed or silently averaged. Unknown metadata also fails.
    """
    import warnings
    if df.empty or len(df.columns) < 2:
        raise ValueError(f"Empty abundance table: {file_name}")
    df = df.copy()
    if not df.columns.is_unique:
        raise ValueError(f"Duplicate input column labels: {file_name}")
    if sample_id_col not in df:
        df = df.rename(columns={df.columns[0]: sample_id_col})
    ids = df[sample_id_col]
    if ids.isna().any() or ids.astype(str).str.strip().eq("").any():
        raise ValueError(f"Missing sample ID: {file_name}")
    df[sample_id_col] = ids.astype(str).str.strip()
    # A reference-style taxonomy row is annotation, never a biological sample.
    df = df.loc[df[sample_id_col].str.casefold() != "taxonomy"].copy()
    if df[sample_id_col].duplicated().any():
        raise ValueError(f"Repeated sample IDs (not silently pooled): {file_name}")
    if group_col not in df:
        labels = [g for g in ("Control", "Pregnant")
                  if g.casefold() in file_name.casefold().split("_")]
        if len(labels) != 1:
            raise ValueError(f"Cannot infer group from {file_name}")
        df[group_col] = labels[0]
    df = df.set_index(sample_id_col)
    groups = df[group_col].astype(str).str.strip().str.casefold().map(
        {"control": "Control", "pregnant": "Pregnant"})
    if groups.isna().any():
        raise ValueError(f"Missing/unsupported group labels: {file_name}")
    # All abundance fields, including descendants, must be numeric and valid.
    features = df.drop(columns=group_col)
    import re
    unknown = [c for c in features if not re.match(r"^[kd]__", str(c).strip())]
    if unknown:
        raise ValueError(f"Unexpected metadata/feature columns in {file_name}: {unknown[:5]}")
    abundance = features.apply(pd.to_numeric, errors="raise")
    if not np.isfinite(abundance.to_numpy(dtype=float)).all() or (abundance < 0).any().any():
        raise ValueError(f"Missing, infinite or negative abundance: {file_name}")
    mapping = {c: extract_family(c) for c in abundance}
    mapping = {c: p for c, p in mapping.items() if p is not None}
    if not mapping:
        raise ValueError(f"No explicit family-level columns in {file_name}; no rank fallback")
    selected = abundance[list(mapping)].copy()
    selected.columns = [mapping[c] for c in selected]
    for taxon in selected.columns.unique():
        block = selected.loc[:, selected.columns == taxon].to_numpy(dtype=float)
        if block.shape[1] > 1 and not np.allclose(block, block[:, :1], rtol=0, atol=1e-10):
            raise ValueError(f"Conflicting duplicate family {taxon} in {file_name}")
    # Same mean grouping as the reference, applied to non-overlapping ranks.
    processed = selected.T.groupby(level=0, sort=True).mean().T
    totals = processed.sum(axis=1)
    valid = totals > 0
    if (~valid).any():
        warnings.warn(f"{file_name}: excluded {(~valid).sum()} zero-family-sum samples", stacklevel=2)
    processed, groups = processed.loc[valid], groups.loc[valid]
    if processed.empty:
        raise ValueError(f"No positive-family-sum samples in {file_name}")
    rel = processed.div(processed.sum(axis=1), axis=0) * 100
    # Reference's second closure is numerically redundant; keep it explicitly.
    rel = rel.div(rel.sum(axis=1), axis=0) * 100
    if not np.allclose(rel.sum(axis=1), 100, atol=1e-8):
        raise AssertionError("Relative abundance closure failed")
    assert rel.index.equals(groups.index)
    print(f"{file_name}: {len(df)} samples -> {len(rel)} retained; "
          f"{len(mapping)} explicit family features -> {rel.shape[1]} families")
    return rel, groups


def main(input_folder=input_folder, output_folder=output_folder):
    os.makedirs(output_folder, exist_ok=True)
    csv_files = sorted(glob.glob(os.path.join(input_folder, "*_SG_formatted.csv")))
    if not csv_files:
        raise FileNotFoundError(f"No formatted abundance tables found in {input_folder}")
    print(f"Found {len(csv_files)} CSV files")

    # ===== BUILD GLOBAL COLOR MAP =====
    all_top_taxa = []

    for file_path in csv_files:
        file_name = os.path.splitext(os.path.basename(file_path))[0]
        print(f"Scanning for colors: {file_name}")

        df = pd.read_csv(file_path)
        rel, groups = prepare_relative_family_table(df, sample_id_col, group_col, file_name)

        if rel is None:
            print(f"Skipping {file_name} in color-map stage")
            continue

        family_means = rel.mean(axis=0).sort_values(ascending=False)
        top_families = family_means.head(top_n).index.tolist()
        all_top_taxa.extend(top_families)

    # unique taxa in stable order
    global_taxa_order = list(dict.fromkeys(all_top_taxa))

    colors = (
        list(plt.cm.tab20.colors) +
        list(plt.cm.Set3.colors) +
        list(plt.cm.Pastel1.colors)
    )

    global_color_map = {
        taxon: colors[i % len(colors)]
        for i, taxon in enumerate(global_taxa_order)
    }

    # Make "Other" always gray
    global_color_map["Other"] = "#BDBDBD"

    # ===== MAIN LOOP =====
    all_records = []

    for file_path in csv_files:
        file_name = os.path.splitext(os.path.basename(file_path))[0]
        title_name = file_name

        print(f"Processing: {file_name}")

        df = pd.read_csv(file_path)
        rel, groups = prepare_relative_family_table(df, sample_id_col, group_col, file_name)

        if rel is None:
            print(f"Skipping {file_name}: could not prepare relative family table")
            continue

        # ===== KEEP TOP FAMILIES FOR THIS DATASET =====
        family_means = rel.mean(axis=0).sort_values(ascending=False)
        top_families = family_means.head(top_n).index.tolist()

        # --------- collect records for combined violin ----------
        for taxon in global_taxa_order:
            if taxon not in rel.columns:
                continue  # Unmeasured taxa are missing, never fabricated zeros.
            tmp = pd.DataFrame({
                "dataset": file_name,
                "sample_id": rel.index,
                "group": groups.values,
                "taxon": taxon,
                "value": rel[taxon].values
            })
            all_records.append(tmp)

        # ===== STACKED PLOT FOR THIS DATASET =====
        rel_top = rel[top_families].copy()

        other_families = [c for c in rel.columns if c not in top_families]
        if other_families:
            rel_top["Other"] = rel[other_families].sum(axis=1)

        # Sort taxa by mean abundance
        rel_top = rel_top[rel_top.mean(axis=0).sort_values(ascending=False).index]

        # Sort samples by dominant family
        dominant_family = family_means.index[0]
        if dominant_family in rel_top.columns:
            rel_top = rel_top.sort_values(by=dominant_family, ascending=False)

        # Final order
        final_order = rel_top.mean(axis=0).sort_values(ascending=False).index.tolist()
        rel_top = rel_top[final_order]

        # ===== COLORS =====
        color_map = {taxon: global_color_map[taxon] for taxon in rel_top.columns}

        # ===== STACKED PLOT =====
        fig, axes = plt.subplots(
            1, 2,
            figsize=(16, 6),
            gridspec_kw={"width_ratios": [5, 1]}
        )

        ax_samples, ax_mean = axes

        # Left panel: samples
        x_positions = range(len(rel_top))
        bottom = [0] * len(rel_top)

        for taxon in rel_top.columns:
            values = rel_top[taxon].values
            ax_samples.bar(
                x_positions,
                values,
                bottom=bottom,
                width=0.9,
                color=color_map[taxon],
                edgecolor="none"
            )
            bottom = [b + v for b, v in zip(bottom, values)]

        ax_samples.set_ylabel("Relative abundance (%)", fontsize=12)
        ax_samples.set_ylim(0, 100)
        ax_samples.set_xticks([])
        ax_samples.set_title("Samples", fontsize=13)
        ax_samples.spines["top"].set_visible(False)
        ax_samples.spines["right"].set_visible(False)

        # Right panel: mean
        mean_values = rel_top.mean(axis=0)
        bottom = 0

        for taxon in rel_top.columns:
            value = mean_values[taxon]
            ax_mean.bar(
                [0],
                [value],
                bottom=bottom,
                width=0.9,
                color=color_map[taxon],
                edgecolor="none"
            )
            bottom += value

        ax_mean.set_ylim(0, 100)
        ax_mean.set_xticks([0])
        ax_mean.set_xticklabels(["Mean"], fontsize=11)
        ax_mean.set_title("Average", fontsize=13)
        ax_mean.spines["top"].set_visible(False)
        ax_mean.spines["right"].set_visible(False)
        ax_mean.spines["left"].set_visible(False)
        ax_mean.yaxis.set_visible(False)

        # Legend
        legend_handles = [
            Rectangle((0, 0), 1, 1, color=color_map[taxon], label=taxon.rsplit(";", 1)[-1])
            for taxon in rel_top.columns
        ]
        fig.legend(
            handles=legend_handles,
            loc="lower center",
            bbox_to_anchor=(0.5, -0.12),
            ncol=4,
            frameon=False,
            fontsize=9
        )

        fig.suptitle(title_name, fontsize=16, fontweight="bold")
        plt.tight_layout(rect=[0, 0.12, 1, 1])

        out_path = os.path.join(output_folder, f"{file_name}_family_composition.png")
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close()

        print(f"Saved stacked plot to: {out_path}")

    # ===== COMBINED SPLIT VIOLIN PLOT =====
    if len(all_records) == 0:
        print("No valid records collected for combined violin plot.")
    else:
        violin_df = pd.concat(all_records, ignore_index=True)
        if violin_df.duplicated(["dataset", "sample_id", "taxon"]).any():
            raise AssertionError("Repeated sample/taxon observations")
        violin_df.to_csv(Path(output_folder) / "violin_observations_SG.csv", index=False)

        violin_df["group"] = violin_df["group"].astype(str).str.strip()
        violin_df["group"] = violin_df["group"].replace({
            "pregnant": "Pregnant",
            "Pregnant": "Pregnant",
            "control": "Control",
            "Control": "Control"
        })

        violin_df = violin_df[violin_df["group"].isin(["Control", "Pregnant"])]

        missing_groups = {"Control", "Pregnant"} - set(violin_df["group"])
        if missing_groups:
            import warnings
            warnings.warn("Combined violin has no observations for: " +
                          ", ".join(sorted(missing_groups)) +
                          "; no two-group comparison is possible.", stacklevel=2)

        if violin_df.empty:
            print("No Control/Pregnant rows found for violin plot.")
        else:
            taxon_order = (
                violin_df.groupby("taxon")["value"]
                .mean()
                .sort_values(ascending=False)
                .index
                .tolist()
            )

            fig, ax = plt.subplots(figsize=(max(14, len(taxon_order) * 0.8), 7))

            control_color = "#4C72B0"
            pregnant_color = "#DD8452"

            for i, taxon in enumerate(taxon_order):
                control_vals = violin_df[
                    (violin_df["taxon"] == taxon) &
                    (violin_df["group"] == "Control")
                ]["value"].dropna().values

                pregnant_vals = violin_df[
                    (violin_df["taxon"] == taxon) &
                    (violin_df["group"] == "Pregnant")
                ]["value"].dropna().values

                center = i
            # KDE is undefined for singleton/constant distributions.
            for vals, side, color in ((control_vals, -1, control_color),
                                      (pregnant_vals, 1, pregnant_color)):
                if len(vals) and (len(vals) == 1 or np.ptp(vals) == 0):
                    ax.plot([center, center + side * 0.32], [vals[0], vals[0]],
                            color=color, linewidth=2)


                # Control = left half
                if len(control_vals) > 1 and np.ptp(control_vals) > 0:
                    parts = ax.violinplot(
                        [control_vals],
                        positions=[center],
                        widths=0.8,
                        showmeans=False,
                        showmedians=True,
                        showextrema=False
                    )

                    for body in parts["bodies"]:
                        body.set_facecolor(control_color)
                        body.set_edgecolor("black")
                        body.set_alpha(0.7)

                        verts = body.get_paths()[0].vertices
                        verts[:, 0] = np.minimum(verts[:, 0], center)

                    if "cmedians" in parts:
                        med = parts["cmedians"]
                        med.set_color("black")
                        med.set_linewidth(1.2)

                # Pregnant = right half
                if len(pregnant_vals) > 1 and np.ptp(pregnant_vals) > 0:
                    parts = ax.violinplot(
                        [pregnant_vals],
                        positions=[center],
                        widths=0.8,
                        showmeans=False,
                        showmedians=True,
                        showextrema=False
                    )

                    for body in parts["bodies"]:
                        body.set_facecolor(pregnant_color)
                        body.set_edgecolor("black")
                        body.set_alpha(0.7)

                        verts = body.get_paths()[0].vertices
                        verts[:, 0] = np.maximum(verts[:, 0], center)

                    if "cmedians" in parts:
                        med = parts["cmedians"]
                        med.set_color("black")
                        med.set_linewidth(1.2)

            ax.set_xticks(range(len(taxon_order)))
            ax.set_xticklabels([taxon.rsplit(";", 1)[-1] for taxon in taxon_order],
                               rotation=45, ha="right")
            ax.set_ylabel("Relative abundance (%)")
            ax.set_title("Top families across SG datasets: Control vs Pregnant")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            legend_handles = [
                Patch(facecolor=control_color, edgecolor="black", label="Control"),
                Patch(facecolor=pregnant_color, edgecolor="black", label="Pregnant")
            ]
            ax.legend(handles=legend_handles, frameon=False)

            plt.tight_layout()

            violin_out_path = os.path.join(output_folder, "violin_plot_SG.png")
            plt.savefig(violin_out_path, dpi=300, bbox_inches="tight")
            plt.close()

            print(f"Saved combined split violin plot to: {violin_out_path}")

    print("Done.")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-folder", type=Path, default=input_folder)
    parser.add_argument("--output-folder", type=Path, default=output_folder)
    args = parser.parse_args()
    main(args.input_folder, args.output_folder)
