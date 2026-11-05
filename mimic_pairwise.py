"""Run independent miMic pregnancy/control comparisons and index their outputs."""
from contextlib import redirect_stdout, redirect_stderr
from itertools import product
import traceback
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
import pandas as pd


def run_pairwise(run_comparison, projects, output_root, alpha):
    output_root.mkdir(parents=True, exist_ok=True)
    summary = []
    for control, pregnant in product(projects["Control"], projects["Pregnant"]):
        folder = output_root / f"{control}_Control_VS_{pregnant}_Pregnant"
        folder.mkdir(parents=True, exist_ok=True)
        row = {"control_project": control, "pregnancy_project": pregnant,
               "alpha": alpha, "output_dir": str(folder)}
        print(f"Running {control} (Control) vs {pregnant} (Pregnant)", flush=True)
        with (folder / "run.log").open("w") as log, redirect_stdout(log), redirect_stderr(log):
            try:
                row.update(run_comparison({"Control": [control], "Pregnant": [pregnant]}, folder))
                if row["starting_taxonomy"] in (None, "nosignificant"):
                    row.update(status="no_significant_taxa", plotted_taxa=0, figure_pdf="")
                else:
                    tree = folder / "plots/correlations_tree.svg"
                    labels = ET.parse(tree).getroot().findall(".//{http://www.w3.org/2000/svg}text")
                    row.update(status="complete", plotted_taxa=len(labels),
                               figure_pdf=str(tree.with_suffix(".pdf")))
            except Exception as error:
                traceback.print_exc()
                row.update(status="error", error=str(error))
            finally:
                plt.close("all")
        summary.append(row)
        pd.DataFrame(summary).to_csv(output_root / "comparison_summary.csv", index=False)
        print(f"  {row['status']}: {row.get('plotted_taxa', '?')} plotted taxa; {folder}", flush=True)
    if any(row["status"] == "error" for row in summary):
        raise RuntimeError(f"Some comparisons failed; see {output_root / 'comparison_summary.csv'}")
