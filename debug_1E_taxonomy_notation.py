import re
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd


BASE_DIR = Path("/home/aharonox/Yoram_Omri_Preg")

FILES = {
    "Control_PRJNA669650": {
        "type": "yamas",
        "folder": BASE_DIR / "datasets_after_yamas/16S/PRJNA669650",
    },
    "Pregnant_PRJNA1254708": {
        "type": "yamas",
        "folder": BASE_DIR / "datasets_after_yamas/16S/PRJNA1254708",
    },
    "Pregnant_Omri_stool": {
        "type": "preprocessed",
        "path": BASE_DIR / "datasets_after_MIPMLP/omri/omri_stool_Pregnant_16S.csv",
    },
}


def normalize_taxonomy_name(tax, split_pipe=False):
    tax = str(tax).strip()

    if tax.lower() in ["nan", "none", "", "0", "unassigned"]:
        return None

    # current script splits on ; and ,
    # pipe-aware version also splits on |
    splitter = r"[;,|]\s*" if split_pipe else r"[;,]\s*"
    parts = re.split(splitter, tax)

    level_order = ["k", "p", "c", "o", "f", "g", "s"]
    level_map = {
        "d": "k",
        "k": "k",
        "p": "p",
        "c": "c",
        "o": "o",
        "f": "f",
        "g": "g",
        "s": "s",
    }

    tax_dict = {}

    for p in parts:
        p = p.strip()

        if not p or p.lower() in ["nan", "none", "unassigned"]:
            continue

        p = re.sub(r"^D_0__", "k__", p)
        p = re.sub(r"^D_1__", "p__", p)
        p = re.sub(r"^D_2__", "c__", p)
        p = re.sub(r"^D_3__", "o__", p)
        p = re.sub(r"^D_4__", "f__", p)
        p = re.sub(r"^D_5__", "g__", p)
        p = re.sub(r"^D_6__", "s__", p)

        m = re.match(r"^([dkpcfgos]):(.+)$", p)
        if m:
            code = level_map.get(m.group(1))
            value = m.group(2).strip()
            p = f"{code}__{value}"

        if p.startswith("d__"):
            p = "k__" + p.replace("d__", "", 1)

        m = re.match(r"^([kpcofgs])__(.*)$", p)
        if m:
            code = m.group(1)
            value = m.group(2).strip()

            if value == "":
                value = "unclassified"

            value = value.strip("[]")
            tax_dict[code] = value

    if not tax_dict:
        plain_parts = [
            p.strip().strip("[]")
            for p in parts
            if p.strip() and p.strip().lower() not in ["nan", "none", "unassigned"]
        ]
        for code, value in zip(level_order, plain_parts):
            tax_dict[code] = value

    if not tax_dict:
        return None

    last_value = "unclassified"
    keep = []

    for code in level_order:
        value = tax_dict.get(code)

        if value is None or value == "":
            value = f"{last_value}_unclassified"
        else:
            last_value = value

        keep.append(f"{code}__{value}")

    return ";".join(keep)


def cut_level(tax, level):
    return ";".join(str(tax).split(";")[:level])


def raw_taxa_from_source(info):
    if info["type"] == "yamas":
        folder = info["folder"]
        tax_path = next(folder.glob("taxonomy_*.csv"))
        tax = pd.read_csv(tax_path, index_col=0)
        return tax["Taxon"].dropna().astype(str).tolist()

    path = info["path"]
    df = pd.read_csv(path, nrows=3, index_col=0)
    cols = [str(c).strip() for c in df.columns]
    cols = [c for c in cols if c.lower() not in ["taxonomy", "taxon", "id", "sampleid", "sample_id", "sample"]]
    cols = [c for c in cols if c not in ["Bacteria", "Unassigned"]]
    return cols


def load_table(info, split_pipe=False):
    if info["type"] == "yamas":
        folder = info["folder"]
        otu_path = next(folder.glob("otu_*.csv"))
        tax_path = next(folder.glob("taxonomy_*.csv"))

        otu = pd.read_csv(otu_path, index_col=0)
        tax = pd.read_csv(tax_path, index_col=0)
        tax = tax.reindex(otu.index)

        new_cols = [normalize_taxonomy_name(x, split_pipe=split_pipe) for x in tax["Taxon"]]
        keep = [x is not None for x in new_cols]

        otu = otu.loc[keep]
        new_cols = [x for x in new_cols if x is not None]

        df = otu.T
        df.columns = new_cols
        df.index.name = "SampleID"

        df = df[~df.index.astype(str).str.endswith("_2")]
        df.index = df.index.astype(str).str.replace(r"_1$", "", regex=True)

    else:
        path = info["path"]
        df = pd.read_csv(path, index_col=0)

        df.columns = df.columns.astype(str).str.strip()
        df.columns = df.columns.str.replace(r"\s*;\s*", ";", regex=True)

        df = df.drop(columns=["Bacteria", "Unassigned"], errors="ignore")
        df = df[~df.index.astype(str).str.lower().str.contains("tax")]

        new_cols = [normalize_taxonomy_name(c, split_pipe=split_pipe) for c in df.columns]
        keep = [c is not None for c in new_cols]

        df = df.loc[:, keep]
        df.columns = [c for c in new_cols if c is not None]

    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
    df = df.T.groupby(level=0).sum().T

    return df


def collapse_to_level(df, level):
    out = df.copy()
    out.columns = [cut_level(c, level) for c in out.columns]
    out = out.T.groupby(level=0).sum().T
    return out


def summarize_raw_notation():
    print("\n" + "=" * 100)
    print("RAW TAXONOMY NOTATION SUMMARY")
    print("=" * 100)

    for name, info in FILES.items():
        taxa = raw_taxa_from_source(info)

        print("\n" + "-" * 100)
        print(name)
        print("num raw taxa:", len(taxa))
        print("first 10 raw taxa:")
        for t in taxa[:10]:
            print(" ", repr(t))

        print("\nnotation counts:")
        print("contains ';'      :", sum(";" in t for t in taxa))
        print("contains ','      :", sum("," in t for t in taxa))
        print("contains '|'      :", sum("|" in t for t in taxa))
        print("contains '__'     :", sum("__" in t for t in taxa))
        print("contains 'k__'    :", sum("k__" in t for t in taxa))
        print("contains 'd__'    :", sum("d__" in t for t in taxa))
        print("contains 'D_0__'  :", sum("D_0__" in t for t in taxa))
        print("contains 'p__'    :", sum("p__" in t for t in taxa))

        semicolon_levels = pd.Series([len(str(t).split(";")) for t in taxa]).value_counts().sort_index()
        pipe_levels = pd.Series([len(str(t).split("|")) for t in taxa]).value_counts().sort_index()

        print("\nsemicolon level counts:")
        print(semicolon_levels)

        print("\npipe level counts:")
        print(pipe_levels)


def report_intersections(split_pipe=False):
    title = "PIPE-AWARE NORMALIZER" if split_pipe else "CURRENT SCRIPT NORMALIZER"
    print("\n" + "=" * 100)
    print(title)
    print("=" * 100)

    tables = {name: load_table(info, split_pipe=split_pipe) for name, info in FILES.items()}

    print("\nloaded table shapes after normalization/collapse:")
    for name, df in tables.items():
        print(f"{name}: {df.shape}")

    for level, level_name in [
        (7, "species/full"),
        (6, "genus"),
        (5, "family"),
        (4, "order"),
        (3, "class"),
        (2, "phylum"),
    ]:
        level_tables = {name: collapse_to_level(df, level) for name, df in tables.items()}
        sets = {name: set(df.columns) for name, df in level_tables.items()}

        print("\n" + "-" * 100)
        print(f"LEVEL {level} ({level_name})")

        for a, b in combinations(sets.keys(), 2):
            print(f"{a} ∩ {b}: {len(sets[a] & sets[b])}")

        common = sorted(set.intersection(*sets.values()))
        print("ALL 3 shared:", len(common))

        if len(common) == 0:
            continue

        print("first 10 shared:")
        for t in common[:10]:
            print(" ", t)

        print("\nzero rows after keeping ALL-3 shared taxa:")
        for name, df in level_tables.items():
            sub = df.loc[:, common]
            row_sum = sub.sum(axis=1)
            zero_rows = int((row_sum == 0).sum())
            total = len(row_sum)
            print(
                f"{name}: zero_rows={zero_rows}/{total} "
                f"({zero_rows / total:.1%}), "
                f"mean_row_sum={row_sum.mean():.6g}, "
                f"median_row_sum={row_sum.median():.6g}"
            )


def main():
    print("Checking files:")
    for name, info in FILES.items():
        if info["type"] == "yamas":
            folder = info["folder"]
            print(name, "folder exists:", folder.exists(), folder)
            print("  otu files:", list(folder.glob("otu_*.csv")))
            print("  taxonomy files:", list(folder.glob("taxonomy_*.csv")))
        else:
            path = info["path"]
            print(name, "file exists:", path.exists(), path)

    summarize_raw_notation()
    report_intersections(split_pipe=False)
    report_intersections(split_pipe=True)


if __name__ == "__main__":
    main()
