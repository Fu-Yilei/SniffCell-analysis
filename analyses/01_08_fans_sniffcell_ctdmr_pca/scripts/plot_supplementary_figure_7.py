#!/usr/bin/env python3
"""Plot joint PCA scores for FACS/FANS and SniffCell-assigned methylation."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
SCORES = ANALYSIS_DIR / "data" / "pca_scores.tsv"
VARIANCE = ANALYSIS_DIR / "data" / "pca_explained_variance.tsv"
DISTANCES = ANALYSIS_DIR / "data" / "pca_pair_distances.tsv"
OUTBASE = ANALYSIS_DIR / "results" / "supplementary_figure_7"

COLORS = {
    "Neuron": "#2B83BA",
    "Oligodendrocyte": "#FDAE61",
    "Monocyte": "#1A9850",
    "T-cell_NK-cell_B-cell": "#D7191C",
}
DISPLAY = {
    "Neuron": "Neuron",
    "Oligodendrocyte": "Oligodendrocyte",
    "Monocyte": "Monocyte",
    "T-cell_NK-cell_B-cell": "T/NK/B",
}
MARKERS = {"FANS": "o", "SniffCell": "X"}


def main() -> None:
    scores = pd.read_csv(SCORES, sep="\t")
    variance = pd.read_csv(VARIANCE, sep="\t")
    distances = pd.read_csv(DISTANCES, sep="\t")
    required = {"PC1", "PC2", "source", "group", "pair_key"}
    if not required.issubset(scores.columns):
        raise ValueError(f"PCA table is missing: {sorted(required - set(scores.columns))}")

    pc1 = 100 * float(variance.loc[variance["PC"] == "PC1", "explained_variance_ratio"].iloc[0])
    pc2 = 100 * float(variance.loc[variance["PC"] == "PC2", "explained_variance_ratio"].iloc[0])

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
    fig, ax = plt.subplots(figsize=(8.5, 6.2))

    for _, pair in scores.groupby("pair_key", sort=False):
        if set(pair["source"]) == {"FANS", "SniffCell"}:
            pair = pair.set_index("source")
            ax.plot(
                [pair.loc["FANS", "PC1"], pair.loc["SniffCell", "PC1"]],
                [pair.loc["FANS", "PC2"], pair.loc["SniffCell", "PC2"]],
                color="#9E9E9E",
                linewidth=0.8,
                alpha=0.55,
                zorder=1,
            )

    for _, row in scores.iterrows():
        ax.scatter(
            row["PC1"],
            row["PC2"],
            s=75,
            marker=MARKERS[row["source"]],
            color=COLORS[row["group"]],
            edgecolor="black",
            linewidth=0.7,
            zorder=3,
        )

    ax.axhline(0, color="#BDBDBD", linewidth=0.8)
    ax.axvline(0, color="#BDBDBD", linewidth=0.8)
    ax.grid(color="#E1E1E1", linewidth=0.6)
    ax.set_xlabel(f"PC1 ({pc1:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({pc2:.1f}% variance)")
    ax.text(
        0.02,
        0.98,
        "Matched source pairs\n"
        f"median PC1-PC2 distance = {distances['distance_PC1_PC2'].median():.1f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#CCCCCC"},
    )

    source_handles = [
        plt.Line2D(
            [0], [0], marker=marker, linestyle="", markerfacecolor="#777777",
            markeredgecolor="black", label=("FACS/FANS" if source == "FANS" else source), markersize=9,
        )
        for source, marker in MARKERS.items()
    ]
    group_handles = [
        plt.Line2D(
            [0], [0], marker="o", linestyle="", markerfacecolor=color,
            markeredgecolor="black", label=DISPLAY[group], markersize=8,
        )
        for group, color in COLORS.items()
    ]
    first = ax.legend(handles=source_handles, title="Read source", loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.add_artist(first)
    ax.legend(handles=group_handles, title="Cell group", loc="upper left", bbox_to_anchor=(1.01, 0.68))

    OUTBASE.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(OUTBASE.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(OUTBASE.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUTBASE.with_suffix('.png')}")
    print(f"Wrote {OUTBASE.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
