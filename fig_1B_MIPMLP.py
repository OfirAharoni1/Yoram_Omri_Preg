import os
import glob
import warnings
import traceback
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# ===== PATCH FOR OLD MIPMLP CODE WITH NEW PANDAS =====

_original_series_getitem = pd.Series.__getitem__

def _patched_series_getitem(self, key):
    if isinstance(key, int) and key not in self.index:
        return self.iloc[key]
    return _original_series_getitem(self, key)

pd.Series.__getitem__ = _patched_series_getitem


_original_index_str = pd.Index.str

class _SafeStringAccessor:
    def __get__(self, obj, cls=None):
        if obj is None:
            return self
        return _original_index_str.__get__(obj.astype(str), cls)

pd.Index.str = _SafeStringAccessor()

import MIPMLP


warnings.filterwarnings("ignore")

# ===== SETTINGS =====
input_folder = "/home/aharonox/Yoram_Omri_Preg/datasets_after_yamas/16S"
output_folder = "/home/aharonox/Yoram_Omri_Preg/1B_plots/16S"
os.makedirs(output_folder, exist_ok=True)

sort_taxon = "Bacteria;Verrucomicrobia"


# ===== COLOR PALETTE =====
def create_color_palette():
    return {
        "Archaea;Euryarchaeota": "#535473",
        "Bacteria;Actinobacteria": "#BE95DC",
        "Bacteria;Bacteroidetes": "#EFA2D2",
        "Bacteria;Chloroflexi": "#FFCBE1",
        "Bacteria;Cyanobacteria": "#F9C6AB",
        "Bacteria;Deferribacteres": "#faa007",
        "Bacteria;Firmicutes": "#C8E6C9",
        "Bacteria;Fusobacteria": "#A8D8EA",
        "Bacteria;Proteobacteria": "#90EE90",
        "Bacteria;SR1": "#8b7da8",
        "Bacteria;Spirochaetes": "#db9ec1",
        "Bacteria;Synergistetes": "#fc5d83",
        "Bacteria;TM7": "#EFA2D2",
        "Bacteria;Tenericutes": "#BE95DC",
        "Bacteria;Verrucomicrobia": "#6B8AD7",
    }


# ===== PROCESS ONE FILE =====
def process_single_file(file_path: str) -> pd.DataFrame:
    df = pd.read_csv(file_path)

    # Ensure columns are strings for MIPMLP compatibility
    df.columns = df.columns.astype(str)
    

    print("Before MIPMLP preprocess:")
    print("Shape:", df.shape)
    print("Columns first 5:", df.columns[:5].tolist())
    print("Index first 5:", df.index[:5].tolist())


    df = MIPMLP.preprocess(
        df,
        taxnomy_group="mean",
        normalization="none",
        taxonomy_level=2
    )

    # Keep numeric columns only
    df = df.select_dtypes(include="number")

    # Remove zero-sum samples
    row_sums = df.sum(axis=1)
    df = df.loc[row_sums > 0]

    if df.empty:
        raise ValueError("All samples have zero total abundance after preprocessing.")

    # Convert to relative abundance (%)
    df = df.div(df.sum(axis=1), axis=0) * 100

    # Keep taxonomy level 2 columns only
    taxa_list = [taxa for taxa in df.columns if len(str(taxa).split(";")) == 2]
    df = df[taxa_list]

    if df.shape[1] == 0:
        raise ValueError("No taxonomy level 2 taxa found after preprocessing.")

    return df


# ===== SORT SAMPLES =====
def sort_samples(df: pd.DataFrame, target_taxon: str) -> pd.DataFrame:
    if target_taxon in df.columns:
        sorted_idx = df[target_taxon].sort_values(ascending=False).index
        return df.loc[sorted_idx]

    fallback_taxon = df.mean(axis=0).sort_values(ascending=False).index[0]
    print(f"Sort taxon not found: {target_taxon}")
    print(f"Sorting instead by: {fallback_taxon}")

    sorted_idx = df[fallback_taxon].sort_values(ascending=False).index
    return df.loc[sorted_idx]


# ===== PLOT ONE FILE =====
def plot_single_dataset(df: pd.DataFrame, file_name: str, output_folder: str) -> str:
    color_map = create_color_palette()

    extra_colors = (
        list(plt.cm.tab20.colors) +
        list(plt.cm.Set3.colors) +
        list(plt.cm.Pastel1.colors)
    )

    missing_taxa = [t for t in df.columns if t not in color_map]
    for i, taxa in enumerate(missing_taxa):
        color_map[taxa] = extra_colors[i % len(extra_colors)]

    fig, (ax1, ax2) = plt.subplots(
        1, 2,
        figsize=(18, 7),
        gridspec_kw={"width_ratios": [5, 1]}
    )

    taxa_list = df.columns.tolist()

    # Left panel: individual samples
    bottom = np.zeros(len(df))
    x_pos = np.arange(len(df))

    for taxa in taxa_list:
        values = df[taxa].values
        ax1.bar(
            x_pos,
            values,
            bottom=bottom,
            color=color_map[taxa],
            width=0.8,
            edgecolor="none"
        )
        bottom += values

    ax1.set_ylabel("Relative abundance (%)", fontsize=14)
    ax1.set_ylim(0, 100)
    ax1.set_xlim(-0.5, len(df) - 0.5)
    ax1.set_xticks([])
    ax1.set_title("Samples", fontsize=14)

    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.spines["bottom"].set_visible(False)

    # Right panel: average
    mean_values = df.mean(axis=0)
    bottom_mean = 0

    for taxa in taxa_list:
        value = mean_values[taxa]
        ax2.bar(
            [0],
            [value],
            bottom=bottom_mean,
            color=color_map[taxa],
            width=0.6,
            edgecolor="none"
        )
        bottom_mean += value

    ax2.set_ylim(0, 100)
    ax2.set_xlim(-0.5, 0.5)
    ax2.set_xticks([0])
    ax2.set_xticklabels(["Average"], fontsize=12)
    ax2.set_title("Average", fontsize=14)

    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    # Legend
    handles = [Rectangle((0, 0), 1, 1, color=color_map[taxa]) for taxa in taxa_list]

    fig.legend(
        handles,
        taxa_list,
        bbox_to_anchor=(0.5, -0.02),
        loc="lower center",
        ncol=3,
        fontsize=10,
        frameon=True,
        columnspacing=1.0,
        handletextpad=0.5
    )

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.24)

    out_path = os.path.join(output_folder, f"{file_name}_composition.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()

    return out_path


# ===== MAIN =====
def main():
    csv_files = sorted(glob.glob(os.path.join(input_folder, "*.csv")))
    print(f"Found {len(csv_files)} CSV files")

    for file_path in csv_files:
        file_name = os.path.splitext(os.path.basename(file_path))[0]
        print(f"\nProcessing: {file_name}")

        try:
            df = process_single_file(file_path)
            df = sort_samples(df, sort_taxon)
            out_path = plot_single_dataset(df, file_name, output_folder)
            print(f"Saved: {out_path}")

        except Exception as e:
            print(f"Skipping {file_name}: {repr(e)}")
            traceback.print_exc()

    print("\nDone.")


if __name__ == "__main__":
    main()