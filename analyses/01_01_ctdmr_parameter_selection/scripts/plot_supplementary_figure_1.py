#!/usr/bin/env python3
"""Reproduce SniffCell Supplementary Figure 1 from aggregate public data."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter, FuncFormatter


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ANALYSIS_DIR / "data" / "supplementary_figure_1_parameter_summary.tsv"
DEFAULT_OUTPUT_DIR = ANALYSIS_DIR / "results"

REQUIRED_COLUMNS = {
    "parameter",
    "value",
    "brain_similarity",
    "pbmc_similarity",
    "brain_ctdmr_count",
    "pbmc_ctdmr_count",
    "total_ctdmr_count",
}
EXPECTED_VALUES = {
    "min_cpgs": [2.0, 3.0, 4.0, 5.0],
    "max_gap_bp": [250.0, 500.0, 1000.0, 2000.0],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot the minimum-CpG and maximum-gap ctDMR benchmark panels."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="Path to the aggregate parameter-summary TSV.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for supplementary_figure_1.png and .pdf.",
    )
    return parser.parse_args()


def load_data(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path, sep="\t")
    missing = REQUIRED_COLUMNS.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")

    if data[list(REQUIRED_COLUMNS - {"parameter"})].isna().any().any():
        raise ValueError("The input contains missing numeric values.")

    for parameter, expected in EXPECTED_VALUES.items():
        observed = sorted(
            data.loc[data["parameter"].eq(parameter), "value"].astype(float).tolist()
        )
        if observed != expected:
            raise ValueError(
                f"Unexpected {parameter} values: observed {observed}, expected {expected}"
            )

    calculated_total = data["brain_ctdmr_count"] + data["pbmc_ctdmr_count"]
    if not np.allclose(calculated_total, data["total_ctdmr_count"], rtol=0, atol=1e-8):
        raise ValueError("total_ctdmr_count does not equal brain + PBMC counts.")

    return data


def format_thousands(value: float) -> str:
    return f"{value / 1000.0:.1f}k"


def configure_style() -> None:
    plt.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 14,
            "axes.titlesize": 15,
            "axes.labelsize": 14,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "legend.fontsize": 12,
        }
    )
    sns.set_theme(
        style="ticks",
        context="paper",
        rc={
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#444444",
            "grid.color": "#d9d9d9",
            "grid.linewidth": 0.6,
        },
    )


def plot_panel(
    ax: plt.Axes,
    data: pd.DataFrame,
    parameter: str,
    x_label: str,
    panel_letter: str,
    show_region_axis: bool,
) -> tuple[list[object], list[str]]:
    panel = data.loc[data["parameter"].eq(parameter)].sort_values("value")
    x = panel["value"].to_numpy(dtype=float)
    brain = panel["brain_similarity"].to_numpy(dtype=float)
    pbmc = panel["pbmc_similarity"].to_numpy(dtype=float)
    total = panel["total_ctdmr_count"].to_numpy(dtype=float)

    brain_color = "#1f77b4"
    pbmc_color = "#ff7f0e"
    region_color = "#3b3b3b"

    region_ax = ax.twinx()
    region_ax.plot(
        x,
        total,
        color=region_color,
        marker="s",
        markersize=4,
        linewidth=1.7,
        alpha=0.95,
        label="Total regions",
        zorder=2,
    )

    dense_x = np.linspace(float(x.min()), float(x.max()), 300)
    ax.plot(
        dense_x,
        np.interp(dense_x, x, brain),
        color=brain_color,
        linestyle="--",
        linewidth=2.1,
        label="Brain similarity",
        zorder=3,
    )
    ax.plot(
        dense_x,
        np.interp(dense_x, x, pbmc),
        color=pbmc_color,
        linestyle="--",
        linewidth=2.1,
        label="PBMC similarity",
        zorder=3,
    )
    ax.scatter(x, brain, color=brain_color, edgecolor="white", linewidth=0.5, s=54, zorder=4)
    ax.scatter(x, pbmc, color=pbmc_color, edgecolor="white", linewidth=0.5, s=54, zorder=4)

    label_indexes = {0, len(x) - 1}
    if parameter == "max_gap_bp":
        label_indexes.add(2)
    x_span = float(max(x.max() - x.min(), 1.0))
    for index, (x_value, region_count) in enumerate(zip(x, total)):
        if index not in label_indexes:
            continue
        x_offset = -0.025 * x_span if index == len(x) - 1 else 0.0
        horizontal_alignment = "right" if index == len(x) - 1 else "left"
        region_ax.text(
            x_value + x_offset,
            region_count,
            format_thousands(region_count),
            fontsize=9.4,
            color=region_color,
            ha=horizontal_alignment,
            va="bottom",
            bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none", "pad": 0.2},
        )

    ax.set_xticks(EXPECTED_VALUES[parameter])
    ax.text(
        0.0,
        1.03,
        f"{panel_letter})",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_xlabel(x_label, fontsize=14.5, labelpad=5)
    ax.set_ylabel("Similarity", fontsize=14.5, labelpad=1)
    ax.set_ylim(0.84, 1.005)
    ax.set_yticks(np.arange(0.84, 1.001, 0.02))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.margins(x=0.05)
    ax.grid(alpha=0.12, axis="y", linewidth=0.6)
    ax.spines["top"].set_visible(False)
    region_ax.spines["top"].set_visible(False)
    ax.tick_params(axis="both", labelsize=12)

    region_ax.set_ylim(0, float(total.max() * 1.10))
    region_ax.yaxis.set_major_formatter(
        FuncFormatter(lambda value, _position: "" if abs(value) < 1e-6 else f"{int(round(value / 1000.0))}k")
    )
    region_ax.grid(False)
    if show_region_axis:
        region_ax.set_ylabel(
            "Total regions (brain+PBMC)", color=region_color, fontsize=14
        )
        region_ax.tick_params(axis="y", colors=region_color, labelsize=11.2)
        region_ax.spines["right"].set_color(region_color)
    else:
        region_ax.tick_params(axis="y", right=False, labelright=False)
        region_ax.spines["right"].set_visible(False)

    handles_left, labels_left = ax.get_legend_handles_labels()
    handles_right, labels_right = region_ax.get_legend_handles_labels()
    return handles_left + handles_right, labels_left + labels_right


def make_figure(data: pd.DataFrame, output_dir: Path) -> None:
    configure_style()
    figure = plt.figure(figsize=(13.4, 4.23), constrained_layout=False)
    grid = figure.add_gridspec(
        1,
        2,
        left=0.065,
        right=0.954,
        top=0.91,
        bottom=0.24,
        wspace=0.16,
    )
    axes = [figure.add_subplot(grid[0, 0]), figure.add_subplot(grid[0, 1])]

    handles, labels = plot_panel(
        axes[0], data, "min_cpgs", "min CpGs", "A", show_region_axis=False
    )
    plot_panel(
        axes[1], data, "max_gap_bp", "max gap (bp)", "B", show_region_axis=True
    )

    unique_handles: dict[str, object] = {}
    for handle, label in zip(handles, labels):
        if label and label not in unique_handles:
            unique_handles[label] = handle
    figure.legend(
        list(unique_handles.values()),
        list(unique_handles.keys()),
        loc="lower center",
        bbox_to_anchor=(0.5, 0.015),
        ncol=3,
        fontsize=10.5,
        frameon=False,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = output_dir / "supplementary_figure_1"
    figure.savefig(
        output_prefix.with_suffix(".png"), dpi=400, bbox_inches="tight", facecolor="white"
    )
    figure.savefig(
        output_prefix.with_suffix(".pdf"),
        bbox_inches="tight",
        facecolor="white",
        metadata={
            "Creator": "SniffCell-analysis",
            "Producer": "Matplotlib",
            "CreationDate": None,
            "ModDate": None,
        },
    )
    plt.close(figure)


def main() -> int:
    args = parse_args()
    data = load_data(args.input)
    make_figure(data, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
