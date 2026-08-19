#!/usr/bin/env python3
"""Plot leaf-resolved read-level disagreement with matched FANS references."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


PRETTY = {
    "Neuron": "Neuron",
    "Oligodendrocyte": "Oligodendrocyte",
    "Monocyte": "Monocyte",
    "T-cell": "T-cell",
    "NK-cell": "NK-cell",
    "B-cell": "B-cell",
    "Lymph-multi": "Lymphocyte (multi/unresolved)",
}
TICK = {
    "Neuron": "Neuron",
    "Oligodendrocyte": "Oligo",
    "Monocyte": "Mono",
    "T-cell": "T-cell",
    "NK-cell": "NK-cell",
    "B-cell": "B-cell",
    "Lymph-multi": "Lymph-multi",
}
COLOR = {
    "Neuron": "#4C72B0",
    "Oligodendrocyte": "#55A868",
    "Monocyte": "#C44E52",
    "T-cell": "#8172B3",
    "NK-cell": "#CCB974",
    "B-cell": "#DA8BC3",
    "Lymph-multi": "#999999",
}
BRAIN_NAME = {"b1": "bcontrol1", "b2": "bmsa", "b3": "bparkinsons", "b4": "bcontrol2"}
BRAIN_ORDER = ["b1", "b4", "b2", "b3"]


def parse_args() -> argparse.Namespace:
    analysis_dir = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=analysis_dir / "data/read_level_fans_disagreement.tsv")
    parser.add_argument("--output", type=Path, default=analysis_dir / "results/supplementary_figure_5")
    return parser.parse_args()


def pooled(rows: list[tuple[str, str, float, int, int]]) -> float:
    observations = sum(row[3] for row in rows)
    disagreements = sum(row[2] / 100 * row[3] for row in rows)
    return 100 * disagreements / observations


def main() -> int:
    args = parse_args()
    rows = []
    with args.input.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append(
                (
                    row["sample"],
                    row["group"],
                    float(row["pct_gt0.4"].rstrip("%")),
                    int(row["n_obs"]),
                    int(row["total_reads"]),
                )
            )

    brain = [row for sample in BRAIN_ORDER for row in rows if row[0] == sample]
    pbmc = [row for row in rows if row[0].startswith("p")]
    brain_pool, pbmc_pool, overall = pooled(brain), pooled(pbmc), pooled(rows)

    fig, (axp, axb) = plt.subplots(
        1, 2, figsize=(13, 4.2), gridspec_kw={"width_ratios": [len(pbmc), len(brain)]}
    )
    for axis, data, title in (
        (axp, pbmc, f"PBMC (pooled {pbmc_pool:.1f}%)"),
        (axb, brain, f"Brain (pooled {brain_pool:.1f}%)"),
    ):
        labels = [f"{BRAIN_NAME.get(sample, sample)}\n{TICK[group]}" for sample, group, *_ in data]
        values = [value for _, _, value, _, _ in data]
        total_reads = [total for _, _, _, _, total in data]
        colors = [COLOR[group] for _, group, *_ in data]
        x = range(len(data))
        axis.bar(x, values, color=colors, edgecolor="black", linewidth=0.5)
        for xi, value, total in zip(x, values, total_reads):
            axis.text(xi, value + 0.15, f"{value:.1f}%\nn={total / 1e6:.2f}M", ha="center", va="bottom", fontsize=6.5)
        axis.set_xticks(list(x))
        axis.set_xticklabels(labels, fontsize=7)
        axis.set_title(title, fontsize=11)
        axis.set_ylim(0, 13)
        axis.grid(axis="y", alpha=0.3)

    axp.set_ylabel("Classified reads disagreeing with FANS\n(|Δβ| > 0.4 at ctDMR, %)")
    handles = [Patch(facecolor=COLOR[group], edgecolor="black", label=PRETTY[group]) for group in PRETTY]
    fig.legend(handles=handles, loc="upper center", ncol=4, fontsize=8, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout(rect=[0, 0, 1, 0.91])
    fig.subplots_adjust(wspace=0.12)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        fig.savefig(args.output.with_suffix(f".{extension}"), dpi=300)
    plt.close(fig)
    print(f"wrote {args.output}.png/.pdf (brain {brain_pool:.1f}%, PBMC {pbmc_pool:.1f}%)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
