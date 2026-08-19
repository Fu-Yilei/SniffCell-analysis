#!/usr/bin/env python3
"""Plot Figure 2d and Supplementary Figure 8 genomic-region coverage."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from statistics import fmean

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
REGION_INPUT = ANALYSIS_DIR / "data" / "region_base_coverage_summary.tsv"
GTEX_INPUT = ANALYSIS_DIR / "data" / "gtex_tissue_specific_gene_coverage_summary.tsv"
COUNT_INPUT = ANALYSIS_DIR / "data" / "gtex_tissue_specific_gene_counts.tsv"
RESULTS = ANALYSIS_DIR / "results"

TISSUES = ["Brain", "PBMC"]
GENERAL_CATEGORIES = ["GENCODE\ngenes", "MRG", "cCREs", "All other\ngenome"]
MAIN_CATEGORIES = ["GTEx\ntissue-spec", *GENERAL_CATEGORIES]
COLORS = {
    "GTEx\ntissue-spec": "#2F6F98",
    "GENCODE\ngenes": "#4C78A8",
    "MRG": "#59A14F",
    "cCREs": "#C58C3A",
    "All other\ngenome": "#8C8C8C",
}


def load_values() -> tuple[dict[tuple[str, str], float], dict[str, int]]:
    dataset_to_tissue = {"brain": "Brain", "pbmc": "PBMC"}
    region_to_category = {
        "GENCODE genes": "GENCODE\ngenes",
        "Medically relevant genes": "MRG",
        "cCREs": "cCREs",
        "All other genome": "All other\ngenome",
    }
    grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
    with REGION_INPUT.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            grouped[(dataset_to_tissue[row["dataset"]], region_to_category[row["region_set"]])].append(
                float(row["pct_region_bases_ge_5x"])
            )
    values = {key: fmean(items) for key, items in grouped.items()}

    gtex_counts: dict[str, int] = {}
    with GTEX_INPUT.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["gene_class"] != "GTEx tissue-specific":
                continue
            tissue = {"brain_cereb": "Brain", "pbmc": "PBMC"}[row["catalog"]]
            values[(tissue, "GTEx\ntissue-spec")] = float(row["mean_pct_bases_ge_5x"])
            gtex_counts[tissue] = int(row["n_genes"])

    expected_counts = {}
    with COUNT_INPUT.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["catalog"] in {"brain_cereb", "pbmc"}:
                tissue = {"brain_cereb": "Brain", "pbmc": "PBMC"}[row["catalog"]]
                expected_counts[tissue] = int(row["n_matched_tissue_specific_genes"])
    if gtex_counts != expected_counts:
        raise ValueError(f"GTEx gene counts disagree across inputs: {gtex_counts} vs {expected_counts}")
    return values, gtex_counts


def render(
    values: dict[tuple[str, str], float],
    gtex_counts: dict[str, int],
    categories: list[str],
    outbase: Path,
    show_gtex_counts: bool,
) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 15,
            "axes.titlesize": 20,
            "axes.labelsize": 18,
            "xtick.labelsize": 12,
            "ytick.labelsize": 14,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 8), sharey=True)
    x = list(range(len(categories)))
    for ax, tissue in zip(axes, TISSUES):
        y = [values[(tissue, category)] for category in categories]
        labels = categories.copy()
        if show_gtex_counts:
            labels[0] = f"GTEx\ntissue-spec\n(n={gtex_counts[tissue]:,})"
        bars = ax.bar(x, y, color=[COLORS[c] for c in categories], width=0.72)
        ax.set_title(tissue, weight="bold")
        ax.set_xticks(x, labels)
        ax.set_ylim(0, 100)
        ax.set_yticks(range(0, 101, 20))
        ax.grid(axis="y", color="#D2D2D2", linewidth=1)
        ax.grid(axis="x", visible=False)
        ax.spines[["top", "right"]].set_visible(False)
        for bar, value in zip(bars, y):
            ax.text(bar.get_x() + bar.get_width() / 2, value + 1.4, f"{value:.0f}", ha="center", va="bottom")
    axes[0].set_ylabel("Bases covered at >=5x (%)")
    fig.tight_layout()
    outbase.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outbase.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(outbase.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)


def main() -> None:
    values, gtex_counts = load_values()
    render(
        values,
        gtex_counts,
        MAIN_CATEGORIES,
        RESULTS / "figure_2d_genomic_region_coverage",
        show_gtex_counts=True,
    )
    render(
        values,
        gtex_counts,
        GENERAL_CATEGORIES,
        RESULTS / "supplementary_figure_8_genomic_region_coverage",
        show_gtex_counts=False,
    )
    print(f"GTEx tissue-specific genes: Brain n={gtex_counts['Brain']:,}; Blood/PBMC n={gtex_counts['PBMC']:,}")
    print(f"Wrote figures under {RESULTS}")


if __name__ == "__main__":
    main()
