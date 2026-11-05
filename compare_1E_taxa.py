"""Compare existing pooled/pairwise 1E leaf results; no statistical reruns.

Run: python compare_1E_taxa.py (mipmlp_env).
Writes CSV matrices, detailed source records and a searchable HTML table per type.
"""
from html import escape
from pathlib import Path
import xml.etree.ElementTree as ET

import pandas as pd

BASE = Path(__file__).resolve().parent
CONFIG = {
    "16S": (0.05, "Control_VS_Pregnant_five_projects"),
    "SG": (0.01, "Control_VS_Pregnant_four_projects_p0.01"),
}
HIGH = "Higher in pregnancy"
LOW = "Lower in pregnancy"
NEUTRAL = "Significant (no direction)"
NOT_SELECTED = "Not selected"
NOT_INPUT = "Not an input leaf"
SIGNIFICANT = {HIGH, LOW, NEUTRAL}
RANKS = dict(zip("kpcofgs", ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]))


def read_comparison(folder, alpha):
    # SAMBA removes spaces from names. Retain the full path and database suffixes.
    columns = pd.read_csv(folder / "processed_for_mimic.csv", index_col=0, nrows=0).columns
    original = {column.replace(" ", ""): column for column in columns}
    if len(original) != len(columns):
        raise ValueError(f"Ambiguous SAMBA taxonomy keys: {folder}")
    results = pd.read_csv(folder / "df_corrs.csv", index_col=0)
    if not results.index.is_unique:
        raise ValueError(f"Duplicate result keys: {folder}")
    leaves = results.loc[results.index.isin(original) & (results.p < alpha)]
    # Verify that the significant input leaves correspond to the actual figure.
    texts = ET.parse(folder / "plots/correlations_tree.svg").getroot().findall(
        ".//{http://www.w3.org/2000/svg}text")
    displayed = sorted(text.text.rsplit(" (", 1)[0].replace(" ", "") for text in texts)
    expected = sorted(key.split(";")[-1][3:] for key in leaves.index)
    if displayed != expected:
        raise ValueError(f"Figure labels differ from saved leaf results: {folder}")
    return original, leaves


def write_html(path, kind, alpha, table, pvalues, columns):
    colors = {HIGH: "higher", LOW: "lower", NEUTRAL: "neutral",
              NOT_SELECTED: "not-selected", NOT_INPUT: "not-input"}
    rows = []
    for key, row in table.iterrows():
        cells = [f'<td title="{escape(row.taxonomy_path, quote=True)}">{escape(row.taxon)}</td>',
                 f"<td>{row['rank']}</td>", f"<td>{row.significant_pairs}</td>"]
        for column in columns:
            status = row[column]
            p = pvalues.loc[key, column]
            detail = f"<small>p = {p:.3g}</small>" if pd.notna(p) else ""
            cells.append(f'<td class="{colors[status]}">{status}{detail}</td>')
        rows.append("<tr>" + "".join(cells) + "</tr>")
    header = "".join(f"<th>{escape(c)}</th>" for c in ["Taxon", "Rank", "Significant pairs", *columns])
    path.write_text(f'''<!doctype html>
<html lang="en"><meta charset="utf-8"><title>1E {kind}: pooled vs pairs</title>
<style>
body{{font:14px Arial,sans-serif;margin:24px;color:#222}} input{{padding:9px;width:360px}}
.scroll{{overflow:auto;max-height:75vh;margin-top:16px}} table{{border-collapse:separate;border-spacing:0}}
th,td{{padding:9px;border-right:1px solid #ddd;border-bottom:1px solid #ddd;min-width:140px}}
th{{position:sticky;top:0;background:#253449;color:white;z-index:2}}
td:first-child{{position:sticky;left:0;background:white;min-width:240px;z-index:1}}
small{{display:block;margin-top:5px}} .higher{{background:#dbeafe}} .lower{{background:#fee2e2}}
.neutral{{background:#fef3c7}} .not-selected{{background:#f1f5f9}} .not-input{{color:#888}}
.guide{{background:#f8fafc;border:1px solid #dbe3ec;border-radius:8px;padding:16px 20px;
margin:16px 0;line-height:1.7;max-width:1000px}}
.guide h2{{margin:0 0 8px;font-size:18px}} .guide ul{{margin:8px 0;padding-inline-start:24px}}
.guide .badge{{padding:2px 6px;border-radius:4px;display:inline-block}}
</style><h1>1E {kind}: pooled vs individual project pairs</h1>
<section class="guide" lang="en" dir="ltr" aria-labelledby="guide-title">
<h2 id="guide-title">How to read this table</h2>
<p>This table compares taxa shown in the pooled figure and the individual project-pair figures.
Each row represents a taxon displayed in at least one figure.
There are {len(table)} taxa, using a significance threshold of {alpha:g}.</p>
<ul>
<li><b>Taxon</b> — taxon name. Hover over it to see the complete taxonomic path.</li>
<li><b>Rank</b> — taxonomic level, such as species, genus or family.</li>
<li><b>Pooled</b> — the result from the combined analysis of all projects.</li>
<li><b>Pair columns</b> — the result for each project pair. Column headings list the control
project before “vs” and the pregnancy project after it.</li>
<li><b>Significant pairs</b> — the number of pairs in which the taxon was significant,
excluding the pooled comparison. This does not require the same direction across pairs.</li>
</ul>
<p><b>Cell meanings:</b></p>
<ul>
<li><span class="badge higher">Higher in pregnancy</span> — significant, with higher abundance in pregnancy.</li>
<li><span class="badge lower">Lower in pregnancy</span> — significant, with lower abundance in pregnancy.</li>
<li><span class="badge neutral">Significant (no direction)</span> — significant, but the saved direction score is zero.</li>
<li><span class="badge not-selected">Not selected</span> — included as an input leaf, but not selected as significant for the figure.</li>
<li><span class="badge not-input">Not an input leaf</span> — not included as a leaf in that comparison's input.
It may still be an internal tree node; this does not imply biological absence.</li>
</ul>
<p>The displayed p-value is the saved miMic result, including its corrections within each comparison.
Missing values are not zero or one. No global correction across comparisons or cross-project consistency filter is applied.</p>
<details><summary>Which file should I use?</summary>
<p>Use this page for searching and comparing taxa, or open
<a href="taxa_comparison.csv">taxa_comparison.csv</a> in Excel.
<a href="all_input_taxa.csv">all_input_taxa.csv</a> also includes taxa never selected in any figure;
<a href="p_values.csv">p_values.csv</a> contains the saved significance values;
<a href="comparisons.csv">comparisons.csv</a> lists comparison details and figure paths;
<a href="results_long.csv">results_long.csv</a> provides one record per input taxon and comparison.</p>
<p>In the CSV files, <b>input_leaf_pairs</b> counts pairs that included the taxon as an input leaf;
<b>taxonomy_path</b> is the complete taxonomic path, and
<b>taxonomy_key</b> is the matching identifier with spaces removed.</p>
</details>
</section>
<input id="search" placeholder="Search taxon or full taxonomy" aria-label="Search taxa">
<span id="count">{len(table)} rows</span><div class="scroll"><table>
<thead><tr>{header}</tr></thead><tbody>{''.join(rows)}</tbody></table></div>
<script>
const rows = [...document.querySelectorAll('tbody tr')];
document.querySelector('#search').addEventListener('input', event => {{
 const query = event.target.value.toLowerCase(); let count = 0;
 for (const row of rows) {{
  const match = (row.cells[0].title + ' ' + row.cells[0].textContent).toLowerCase().includes(query);
  row.hidden = !match; if (match) count++;
 }}
 document.querySelector('#count').textContent = count + ' rows';
}});
</script></html>''', encoding="utf-8")


def build_tables(kind, alpha, pooled):
    root = BASE / f"1E_{kind}"
    pairs = pd.read_csv(root / f"pairs_p{alpha:g}" / "comparison_summary.csv")
    if not pairs.status.eq("complete").all() or not pairs.alpha.eq(alpha).all():
        raise ValueError(f"Incomplete comparisons or mixed thresholds in {kind}")
    comparisons = [("Pooled", root / pooled)] + [
        (f"{row.control_project} vs {row.pregnancy_project}", Path(row.output_dir))
        for row in pairs.itertuples()]
    records, metadata, originals = [], [], {}
    for label, folder in comparisons:
        original, significant = read_comparison(folder, alpha)
        originals.update(original)
        for key, path in original.items():
            status, p, score = NOT_SELECTED, float("nan"), float("nan")
            if key in significant.index:
                p, score = significant.loc[key, ["p", "scc"]]
                status = HIGH if score > 0 else LOW if score < 0 else NEUTRAL
            records.append({"taxonomy_key": key, "comparison": label, "status": status,
                            "p": p, "scc": score})
        metadata.append({"comparison": label, "alpha": alpha, "input_leaves": len(original),
                         "plotted_taxa": len(significant), "results_file": str(folder / "df_corrs.csv"),
                         "figure_pdf": str(folder / "plots/correlations_tree.pdf")})
    details = pd.DataFrame(records)
    columns = [label for label, _ in comparisons]
    status = details.pivot(index="taxonomy_key", columns="comparison", values="status").reindex(columns=columns)
    status = status.fillna(NOT_INPUT)
    pvalues = details.pivot(index="taxonomy_key", columns="comparison", values="p").reindex(columns=columns)
    # The main table is the union of the plotted taxa, not just the intersection.
    selected = status.isin(SIGNIFICANT).any(axis=1)
    taxa = pd.DataFrame(index=status.index)
    taxa["taxon"] = [originals[key].split(";")[-1][3:] for key in status.index]
    taxa["rank"] = [RANKS[key.split(";")[-1][0]] for key in status.index]
    taxa["taxonomy_path"] = [originals[key] for key in status.index]
    taxa["significant_pairs"] = status[columns[1:]].isin(SIGNIFICANT).sum(axis=1)
    taxa["input_leaf_pairs"] = status[columns[1:]].ne(NOT_INPUT).sum(axis=1)
    all_taxa = taxa.join(status).sort_values(["significant_pairs", "taxon"], ascending=[False, True])
    table = all_taxa.loc[selected.reindex(all_taxa.index)]
    out = root / f"comparison_tables_p{alpha:g}"
    out.mkdir(parents=True, exist_ok=True)
    table.to_csv(out / "taxa_comparison.csv", encoding="utf-8-sig")
    all_taxa.to_csv(out / "all_input_taxa.csv", encoding="utf-8-sig")
    taxa.loc[table.index, ["taxon", "rank", "taxonomy_path"]].join(pvalues).to_csv(
        out / "p_values.csv", encoding="utf-8-sig")
    details.to_csv(out / "results_long.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(metadata).to_csv(out / "comparisons.csv", index=False, encoding="utf-8-sig")
    write_html(out / "taxa_comparison.html", kind, alpha, table, pvalues, columns)
    print(f"{kind}: {len(table)} plotted taxa across {len(columns)} comparisons; {out}")


if __name__ == "__main__":
    for kind, (alpha, pooled) in CONFIG.items():
        build_tables(kind, alpha, pooled)
