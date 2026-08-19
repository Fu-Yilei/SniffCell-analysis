#!/usr/bin/env python3
"""Render Supplementary Figure 3 from path-free directional-effect arrays."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd


ORDER = ["pbmc1", "pbmc2", "b1", "b4", "b2", "b3"]
PRETTY = {
    "pbmc1": "PBMC 1",
    "pbmc2": "PBMC 2",
    "b1": "bcontrol1",
    "b4": "bcontrol2",
    "b2": "bmsa",
    "b3": "bparkinsons",
}


def parse_args() -> argparse.Namespace:
    analysis_dir = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=analysis_dir / "data/supplementary_figure_3_panel_data.npz")
    parser.add_argument("--output", type=Path, default=analysis_dir / "results/supplementary_figure_3")
    return parser.parse_args()


def plot_panel(axis: plt.Axes, data: pd.DataFrame, title: str, x_label: str, y_label: str):
    axis.add_patch(Rectangle((0.0, 0.0), 1.0, 1.0, facecolor="#6fbf73", alpha=0.07, edgecolor="none", zorder=0))
    axis.add_patch(Rectangle((-1.0, -1.0), 1.0, 1.0, facecolor="#6fbf73", alpha=0.07, edgecolor="none", zorder=0))
    axis.axvspan(-0.4, 0.4, color="gray", alpha=0.08, lw=0, zorder=0)
    axis.plot([-1, 1], [-1, 1], linestyle="--", color="#b54b4b", linewidth=0.9, alpha=0.75)
    axis.axhline(0.0, color="#2e8b57", linewidth=1.15, alpha=0.78)
    axis.axvline(0.0, color="#2e8b57", linewidth=1.15, alpha=0.78)
    axis.set_xlim(-1, 1)
    axis.set_ylim(-1, 1)
    axis.set_aspect("equal", adjustable="box")
    axis.set_title(title, loc="left", fontweight="bold", pad=8)
    axis.set_xlabel(x_label)
    axis.set_ylabel(y_label)
    axis.grid(alpha=0.15, linewidth=0.6)
    hexbin = axis.hexbin(data["x"], data["y"], gridsize=70, mincnt=3, cmap="viridis", linewidths=0.0, edgecolors="none")
    correlation = data["x"].corr(data["y"])
    concordance = (np.sign(data["x"]) == np.sign(data["y"])).mean()
    axis.text(
        0.02,
        0.94,
        f"r={correlation:.3f}, directional concordance={concordance:.3f}",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=10.5,
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.75),
    )
    return hexbin


def main() -> int:
    args = parse_args()
    cache = np.load(args.data)
    panels = []
    for index, sample in enumerate(ORDER, start=1):
        pbmc = sample.startswith("pbmc")
        panels.append(
            (
                f"{chr(64 + index)}) {PRETTY[sample]}",
                pd.DataFrame({"x": cache[f"{sample}_x"], "y": cache[f"{sample}_y"]}),
                "Catalog directional effect (best group - other group)" if pbmc else "Catalog directional effect (Neuron - Oligodendrocyte)",
                "Sorted reference effect (best group - other group)" if pbmc else "Sorted reference effect (Neuron - Oligodendrocyte)",
            )
        )

    figure, axes = plt.subplots(2, 3, figsize=(18.6, 11.4), constrained_layout=True, squeeze=False)
    hexbins = []
    for index, (title, data, x_label, y_label) in enumerate(panels):
        hexbins.append(plot_panel(axes[index // 3, index % 3], data, title, x_label, y_label))
    maximum = max(float(item.get_array().max()) for item in hexbins)
    maximum = max(2.0, maximum)
    for item in hexbins:
        item.set_norm(LogNorm(vmin=1.0, vmax=maximum))
    colorbar = figure.colorbar(hexbins[0], ax=axes.ravel().tolist(), fraction=0.026, pad=0.015)
    colorbar.set_label("Region count (log scale)", fontsize=13)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    figure.savefig(args.output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
