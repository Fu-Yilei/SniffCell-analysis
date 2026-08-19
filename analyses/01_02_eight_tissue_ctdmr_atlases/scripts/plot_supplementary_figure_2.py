#!/usr/bin/env python3
"""Reproduce SniffCell Supplementary Figure 2 (panels A-D)."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


TISSUES = {
    "Brain": ("brain_cereb", "brain_cereb.combined_dt04.ctdmr.tsv"),
    "PBMC": ("pbmc", "pbmc.combined_dt04.ctdmr.tsv"),
    "Breast": ("breast", "breast.combined_dt04.ctdmr.tsv"),
    "Colon": ("colon", "colon.combined_dt04.ctdmr.tsv"),
    "Kidney": ("kidney", "kidney.combined_dt04.ctdmr.tsv"),
    "Liver": ("liver", "liver.combined_dt04.ctdmr.tsv"),
    "Lung": ("lung", "lung.combined_dt04.ctdmr.tsv"),
    "Pancreas": ("pancreas", "pancreas.combined_dt04.ctdmr.tsv"),
}
TISSUE_ORDER = list(TISSUES)
TISSUE_COLORS = {
    "Brain": "#1f77b4",
    "PBMC": "#ff7f0e",
    "Breast": "#e377c2",
    "Colon": "#8c564b",
    "Kidney": "#17becf",
    "Liver": "#bcbd22",
    "Lung": "#9467bd",
    "Pancreas": "#2ca02c",
}
CCRE_ORDER = [
    "No overlap",
    "PLS",
    "pELS",
    "dELS",
    "CA-CTCF",
    "CA-H3K4me3",
    "CA-TF",
    "CA",
    "TF",
    "Multi-class overlap",
]
CCRE_COLORS = {
    "No overlap": "#d9d9d9",
    "PLS": "#c0392b",
    "pELS": "#e67e22",
    "dELS": "#f1c40f",
    "CA-CTCF": "#2980b9",
    "CA-H3K4me3": "#8e44ad",
    "CA-TF": "#16a085",
    "CA": "#27ae60",
    "TF": "#7f8c8d",
    "Multi-class overlap": "#34495e",
}


def parse_args() -> argparse.Namespace:
    analysis_dir = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--atlas-dir",
        type=Path,
        default=analysis_dir / "data/sniffcell_find_atlases",
    )
    parser.add_argument(
        "--ccre-summary",
        type=Path,
        default=analysis_dir / "data/figure_ready/tissue_overview_ccre_breakdown.tsv",
    )
    parser.add_argument(
        "--lofyer-summary",
        type=Path,
        default=analysis_dir / "data/figure_ready/lofyer_overlap_summary.tsv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=analysis_dir / "results",
    )
    return parser.parse_args()


def setup_style() -> None:
    plt.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 10,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 8,
        }
    )
    sns.set_theme(
        style="ticks",
        context="paper",
        rc={
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#444444",
            "grid.color": "#e0e0e0",
            "grid.linewidth": 0.45,
        },
    )


def fmt_k(value: float) -> str:
    if value >= 1e6:
        return f"{value / 1e6:.1f}M"
    if value >= 1e3:
        return f"{value / 1e3:.0f}k"
    return str(int(value))


def panel_letter(ax, letter: str, x: float = -0.12, y: float = 1.05) -> None:
    ax.text(
        x,
        y,
        letter,
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        va="top",
        ha="left",
    )


def atlas_path(atlas_dir: Path, tissue_key: str, filename: str) -> Path:
    nested = atlas_dir / tissue_key / filename
    flat = atlas_dir / filename
    if nested.is_file():
        return nested
    if flat.is_file():
        return flat
    raise FileNotFoundError(f"missing atlas: {nested}")


def load_catalog_metrics(
    atlas_dir: Path,
) -> tuple[dict[str, int], dict[str, float], dict[str, float]]:
    counts: dict[str, int] = {}
    coverage_mb: dict[str, float] = {}
    median_len: dict[str, float] = {}
    for tissue, (key, filename) in TISSUES.items():
        frame = pd.read_csv(atlas_path(atlas_dir, key, filename), sep="\t", usecols=["bp_len"])
        lengths = pd.to_numeric(frame["bp_len"], errors="raise")
        counts[tissue] = len(frame)
        coverage_mb[tissue] = float(lengths.sum()) / 1e6
        median_len[tissue] = float(lengths.median())
    return counts, coverage_mb, median_len


def plot_panel_a(ax, counts: dict[str, int]) -> None:
    x = np.arange(len(TISSUE_ORDER))
    ax.bar(
        x,
        [counts[t] for t in TISSUE_ORDER],
        width=0.66,
        color=[TISSUE_COLORS[t] for t in TISSUE_ORDER],
        edgecolor="white",
        linewidth=0.4,
        zorder=3,
    )
    ax.set_xticks(x)
    ax.set_xticklabels(TISSUE_ORDER, rotation=32, ha="right")
    ax.set_ylabel("Number of ctDMRs")
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda value, _: fmt_k(value)))
    ax.yaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    ax.set_ylim(0, max(counts.values()) * 1.08)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    panel_letter(ax, "A")


def plot_panel_b(
    ax,
    coverage_mb: dict[str, float],
    median_len: dict[str, float],
) -> None:
    x = np.arange(len(TISSUE_ORDER))
    ax.bar(
        x,
        [coverage_mb[t] for t in TISSUE_ORDER],
        width=0.66,
        color=[TISSUE_COLORS[t] for t in TISSUE_ORDER],
        edgecolor="white",
        linewidth=0.4,
        zorder=3,
    )
    ax.set_xticks(x)
    ax.set_xticklabels(TISSUE_ORDER, rotation=32, ha="right")
    ax.set_ylabel("Total coverage (Mb)")
    ax.yaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    ax.set_ylim(0, max(coverage_mb.values()) * 1.15)
    ax.spines["top"].set_visible(False)

    ax2 = ax.twinx()
    medians = [median_len[t] for t in TISSUE_ORDER]
    ax2.plot(
        x,
        medians,
        color="#555555",
        marker="D",
        markersize=4.2,
        linewidth=1.2,
        linestyle="--",
        zorder=4,
        label="Median size",
    )
    ax2.set_ylabel("Median size (bp)", color="#555555")
    ax2.tick_params(axis="y", colors="#555555", labelsize=8)
    ax2.set_ylim(max(0, min(medians) * 0.82), max(medians) * 1.22)
    ax2.spines["top"].set_visible(False)
    ax2.legend(loc="upper left", frameon=False, handlelength=1.5)
    panel_letter(ax, "B")


def plot_panel_c(ax, ccre_path: Path, counts: dict[str, int]) -> None:
    ccre = pd.read_csv(ccre_path, sep="\t")
    expected = pd.Series(counts, name="n_total")
    observed = ccre.groupby("tissue", sort=False)["n_total"].first()
    if not observed.reindex(TISSUE_ORDER).equals(expected.reindex(TISSUE_ORDER)):
        raise ValueError("cCRE summary totals do not match the atlas row counts")
    pivot = (
        ccre.pivot_table(
            index="tissue",
            columns="exclusive_category",
            values="fraction",
            aggfunc="sum",
        )
        .fillna(0)
        .reindex(TISSUE_ORDER)
    )
    present = [category for category in CCRE_ORDER if category in pivot.columns]
    bottoms = np.zeros(len(TISSUE_ORDER))
    x = np.arange(len(TISSUE_ORDER))
    text_color = {
        category: (
            "#333333" if category in {"No overlap", "dELS", "pELS"} else "white"
        )
        for category in CCRE_ORDER
    }
    for category in present:
        values = pivot[category].to_numpy()
        ax.bar(
            x,
            values,
            width=0.66,
            bottom=bottoms,
            color=CCRE_COLORS[category],
            label=category,
            edgecolor="none",
            zorder=3,
        )
        for index, (value, bottom) in enumerate(zip(values, bottoms)):
            if value >= 0.10:
                ax.text(
                    index,
                    bottom + value / 2,
                    f"{value:.0%}",
                    ha="center",
                    va="center",
                    fontsize=7,
                    color=text_color[category],
                    fontweight="bold",
                )
        bottoms += values
    ax.set_xticks(x)
    ax.set_xticklabels([f"{t}\n(n={fmt_k(counts[t])})" for t in TISSUE_ORDER])
    ax.set_ylabel("Fraction of ctDMRs")
    ax.set_ylim(0, 1.0)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda value, _: f"{value:.0%}"))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(
        title=None,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.11),
        ncol=10,
        frameon=False,
        columnspacing=0.55,
        handlelength=0.8,
        handletextpad=0.3,
        fontsize=7.4,
    )
    panel_letter(ax, "C", y=1.11)


def load_panel_d(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    if frame["tissue"].tolist() != sorted(
        frame["tissue"].tolist(), key=lambda value: TISSUE_ORDER.index(value)
    ):
        raise ValueError("Loyfer summary tissue order is not the manuscript order")
    return frame


def plot_panel_d(ax, panel_d: pd.DataFrame) -> None:
    y = np.arange(len(panel_d))
    colors = [TISSUE_COLORS[tissue] for tissue in panel_d["tissue"]]
    ax.hlines(
        y,
        60,
        panel_d["coverage_percent"],
        color="#d7d7d7",
        linewidth=0.8,
        zorder=1,
    )
    ax.scatter(
        panel_d["coverage_percent"],
        y,
        s=30,
        c=colors,
        edgecolor="#333333",
        linewidth=0.35,
        zorder=3,
    )
    ax.axvline(100, color="#777777", linewidth=0.6, linestyle=":")
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{row.tissue}: {row.sniffcell_group}" for row in panel_d.itertuples()],
        fontsize=7.2,
    )
    ax.invert_yaxis()
    ax.set_xlim(60, 102)
    ax.set_xticks([60, 70, 80, 90, 100])
    ax.set_xlabel("Loyfer ctDMRs overlapped by SniffCell (%)")
    ax.xaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for yi, percent in zip(y, panel_d["coverage_percent"]):
        ax.text(
            min(percent + 0.8, 101.0),
            yi,
            f"{percent:.1f}",
            va="center",
            ha="left",
            fontsize=6.5,
        )
    n_concordant = int(panel_d["n_concordant_intervals"].sum())
    n_overlap = int(panel_d["n_loyfer_regions_with_any_sniffcell_overlap"].sum())
    if (n_concordant, n_overlap) != (6050, 6050):
        raise ValueError("unexpected directional-concordance totals")
    ax.text(
        0.0,
        1.03,
        f"Direction concordance: {n_concordant}/{n_overlap} overlapped intervals "
        "(U atlas -> SniffCell hypo)",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8,
    )
    panel_letter(ax, "D", x=-0.06, y=1.10)


def main() -> int:
    args = parse_args()
    setup_style()
    counts, coverage_mb, median_len = load_catalog_metrics(args.atlas_dir)
    panel_d = load_panel_d(args.lofyer_summary)

    from matplotlib.gridspec import GridSpec

    figure = plt.figure(figsize=(12.8, 9.8))
    grid = GridSpec(
        3,
        2,
        figure=figure,
        height_ratios=[0.9, 1.08, 1.55],
        left=0.155,
        right=0.985,
        top=0.975,
        bottom=0.075,
        hspace=0.38,
        wspace=0.26,
    )
    plot_panel_a(figure.add_subplot(grid[0, 0]), counts)
    plot_panel_b(figure.add_subplot(grid[0, 1]), coverage_mb, median_len)
    plot_panel_c(figure.add_subplot(grid[1, :]), args.ccre_summary, counts)
    plot_panel_d(figure.add_subplot(grid[2, :]), panel_d)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / "supplementary_figure_2.png"
    figure.savefig(output, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.04)
    figure.savefig(output.with_suffix(".pdf"), facecolor="white", bbox_inches="tight", pad_inches=0.04)
    plt.close(figure)
    print(f"saved {output}")
    print(f"saved {output.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
