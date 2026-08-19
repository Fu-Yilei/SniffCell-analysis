#!/usr/bin/env python3
"""Replot Figure 4g-h from the sanitized donor-normalized source rows."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "data" / "matched_blood_donor_locus_rows.tsv.gz"
COLORS = {"Neuron": "#E67E22", "Oligodendrocyte": "#1B9E77"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-prefix", type=Path, default=ROOT / "results" / "figure_4gh_reproduced")
    args = parser.parse_args()
    data = pd.read_csv(args.input, sep="\t")
    required = {"locus", "donor", "change_cell", "neuron_gain_bp", "oligo_gain_bp", "is_aag_like", "motif"}
    if missing := required - set(data.columns):
        raise ValueError(f"Missing columns: {sorted(missing)}")

    fig, axes = plt.subplots(1, 2, figsize=(18, 14))
    for ax, cell in zip(axes, ["Neuron", "Oligodendrocyte"]):
        part = data[data["change_cell"].eq(cell)].copy()
        order = (part.groupby("locus")
                    .apply(lambda x: (x["neuron_gain_bp"] - x["oligo_gain_bp"]).mean(), include_groups=False)
                    .sort_values(ascending=cell == "Oligodendrocyte").index.tolist())
        ymap = {locus: i for i, locus in enumerate(order)}
        for row in part.itertuples(index=False):
            y = ymap[row.locus]
            ax.plot([row.neuron_gain_bp, row.oligo_gain_bp], [y, y], color="#A0A0A0", linewidth=1)
            ax.scatter(row.neuron_gain_bp, y, color=COLORS["Neuron"], s=24, zorder=3)
            ax.scatter(row.oligo_gain_bp, y, color=COLORS["Oligodendrocyte"], s=24, zorder=3)
        labels = []
        for locus in order:
            row = part[part["locus"].eq(locus)].iloc[0]
            gene = "" if pd.isna(row.get("gene_name")) else f"{row['gene_name']} "
            motif = str(row["motif"])
            motif = motif if len(motif) <= 10 else f"{motif[:6]}…{motif[-3:]}"
            labels.append(f"{gene}{motif}\n{locus}")
        ax.set_yticks(range(len(order)), labels, fontsize=7)
        ax.axvline(0, color="#666666", linewidth=1.4)
        ax.set_xscale("symlog", linthresh=100)
        ax.set_xlabel("Repeat length - matched blood length (bp)")
        ax.set_title(f"{cell}-expanded TRs", color=COLORS[cell], fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="x", color="#EAEAEA", linewidth=.6)
    fig.subplots_adjust(left=.21, right=.98, top=.94, bottom=.08, wspace=.65)
    args.output_prefix.parent.mkdir(exist_ok=True)
    fig.savefig(args.output_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(args.output_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame([{
        "n_donor_locus_rows": len(data), "n_loci": data["locus"].nunique(),
        "n_donors": data["donor"].nunique(), "n_neuron_rows": data["change_cell"].eq("Neuron").sum(),
        "n_oligodendrocyte_rows": data["change_cell"].eq("Oligodendrocyte").sum(),
    }]).to_csv(ROOT / "results" / "matched_blood_source_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
