#!/usr/bin/env python3
"""Plot GTEx expression and gene-body ctDMR enrichment summaries."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
GENES = ANALYSIS_DIR / "data" / "gene_level_expression_ctdmr_density.tsv.gz"
SUMMARY = ANALYSIS_DIR / "data" / "expression_class_ctdmr_coverage_summary.tsv"
OUTBASE = ANALYSIS_DIR / "results" / "supplementary_figure_9"

TISSUE_ORDER = ["Brain (cerebellum + frontal cortex)", "PBMC/whole blood"]
TISSUE_LABEL = {
    "Brain (cerebellum + frontal cortex)": "Brain",
    "PBMC/whole blood": "Blood/PBMC",
}
CLASS_ORDER = ["Expressed", "Non-expressed"]
COLORS = {"Expressed": "#4C78A8", "Non-expressed": "#B8B8B8"}


def spearman(x: pd.Series, y: pd.Series) -> float:
    keep = x.notna() & y.notna()
    return float(x[keep].rank().corr(y[keep].rank()))


def main() -> None:
    genes = pd.read_csv(GENES, sep="\t")
    summary = pd.read_csv(SUMMARY, sep="\t")
    required = {"gene_id", "Tissue", "Expression value", "Expression class", "ctDMR density"}
    if not required.issubset(genes.columns):
        raise ValueError(f"Gene table is missing: {sorted(required - set(genes.columns))}")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    # A: fraction of expressed and non-expressed genes with at least one ctDMR.
    ax = axes[0, 0]
    x = np.arange(len(TISSUE_ORDER))
    width = 0.34
    for offset, expression_class in zip((-width / 2, width / 2), CLASS_ORDER):
        rows = summary.set_index(["Tissue", "Expression class"])
        values = [float(rows.loc[(t, expression_class), "pct_with_ctdmr"]) for t in TISSUE_ORDER]
        counts = [int(rows.loc[(t, expression_class), "n_genes"]) for t in TISSUE_ORDER]
        bars = ax.bar(x + offset, values, width, label=expression_class, color=COLORS[expression_class])
        for bar, value, count in zip(bars, values, counts):
            ax.text(bar.get_x() + bar.get_width() / 2, value + 1.2, f"{value:.1f}%\n(n={count:,})", ha="center", fontsize=7.3)
    ax.set_xticks(x, [TISSUE_LABEL[t] for t in TISSUE_ORDER])
    ax.set_ylabel("Genes with >=1 gene-body ctDMR (%)")
    ax.set_ylim(0, 82)
    ax.legend(frameon=False)
    ax.set_title("A  Gene-body ctDMR coverage", loc="left", weight="bold")

    # B: density distributions among genes that contain a ctDMR.
    ax = axes[0, 1]
    density_groups = []
    density_labels = []
    density_colors = []
    for tissue in TISSUE_ORDER:
        for expression_class in CLASS_ORDER:
            subset = genes.loc[
                (genes["Tissue"] == tissue)
                & (genes["Expression class"] == expression_class)
                & (genes["ctDMR density"] > 0),
                "ctDMR density",
            ]
            density_groups.append(np.log10(subset.to_numpy() + 1))
            density_labels.append(f"{TISSUE_LABEL[tissue]}\n{expression_class}")
            density_colors.append(COLORS[expression_class])
    boxes = ax.boxplot(density_groups, patch_artist=True, showfliers=False, medianprops={"color": "black"})
    for patch, color in zip(boxes["boxes"], density_colors):
        patch.set_facecolor(color)
    ax.set_xticks(range(1, len(density_labels) + 1), density_labels, fontsize=7.5)
    ax.set_ylabel("log10(ctDMR count per Mb + 1)")
    ax.set_title("B  Density among genes with ctDMRs", loc="left", weight="bold")

    # C-D: expression versus density for all genes.
    for ax, tissue, letter in zip(axes[1], TISSUE_ORDER, ["C", "D"]):
        subset = genes.loc[genes["Tissue"] == tissue].copy()
        xval = subset["Expression value"]
        yval = np.log10(subset["ctDMR density"] + 1)
        hb = ax.hexbin(xval, yval, gridsize=55, bins="log", mincnt=1, cmap="Blues")
        rho = spearman(xval, subset["ctDMR density"])
        ax.text(0.98, 0.96, f"Spearman rho = {rho:.2f}\nn = {len(subset):,} genes", transform=ax.transAxes, ha="right", va="top", fontsize=8)
        ax.set_xlabel("GTEx mean log2(TPM + 1)")
        ax.set_ylabel("log10(ctDMR count per Mb + 1)")
        ax.set_title(f"{letter}  {TISSUE_LABEL[tissue]}", loc="left", weight="bold")
        fig.colorbar(hb, ax=ax, label="log10(gene count)", fraction=0.046, pad=0.04)

    for ax in axes.flat:
        ax.grid(axis="y", color="#E0E0E0", linewidth=0.6)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    OUTBASE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTBASE.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(OUTBASE.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)
    print(f"Wrote {OUTBASE.with_suffix('.png')}")
    print(f"Wrote {OUTBASE.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
