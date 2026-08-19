#!/usr/bin/env python3
"""Plot the CARD GAA-family permutation-test supplement."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colormaps
from matplotlib.colors import Normalize, to_rgb
from matplotlib.patches import Patch, Rectangle


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    ROOT / "data/card_gaa_motif_family_recurrence_permutation_tests_155.tsv"
)
DEFAULT_OUTPUT_PREFIX = (
    ROOT / "results/supplementary_figure_18a_reproduced"
)
FAMILY_ORDER = [
    "pure_GAA",
    "GAAA_family",
    "GGAA_family",
    "other_4bp_GAA",
    "other_5_10bp_GAA",
]
THRESHOLDS = [1, 2, 3, 5]
NEURON_COLOR = "#EB7F18"
OLIGO_COLOR = "#1B9E77"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-prefix", type=Path, default=DEFAULT_OUTPUT_PREFIX)
    return parser.parse_args()


def blend_with_white(color: str, intensity: float) -> tuple[float, float, float]:
    rgb = np.asarray(to_rgb(color))
    return tuple(1 - intensity * (1 - rgb))


def format_p(value: float) -> str:
    if value < 0.001:
        return f"{value:.1e}"
    if value < 0.01:
        return f"{value:.4f}"
    return f"{value:.3f}"


def plot_count_matrix(axis: plt.Axes, table: pd.DataFrame) -> None:
    maximum = max(
        table["n_neuron_recurrent_loci"].max(),
        table["n_oligodendrocyte_recurrent_loci"].max(),
    )
    log_maximum = np.log1p(maximum)
    for row_index, family in enumerate(FAMILY_ORDER):
        for column_index, threshold in enumerate(THRESHOLDS):
            record = table.loc[
                table["motif_family"].eq(family)
                & table["min_recurrent_donors"].eq(threshold)
            ].iloc[0]
            for side, count, color in (
                ("left", int(record.n_neuron_recurrent_loci), NEURON_COLOR),
                (
                    "right",
                    int(record.n_oligodendrocyte_recurrent_loci),
                    OLIGO_COLOR,
                ),
            ):
                x = column_index if side == "left" else column_index + 0.5
                intensity = 0.08 + 0.92 * np.log1p(count) / log_maximum
                axis.add_patch(
                    Rectangle(
                        (x, row_index),
                        0.5,
                        1,
                        facecolor=blend_with_white(color, intensity),
                        edgecolor="white",
                        linewidth=1.5,
                    )
                )
                text_color = "white" if intensity >= 0.62 else "#202020"
                axis.text(
                    x + 0.25,
                    row_index + 0.5,
                    str(count),
                    ha="center",
                    va="center",
                    fontsize=12,
                    fontweight="bold",
                    color=text_color,
                )
    axis.set_xlim(0, len(THRESHOLDS))
    axis.set_ylim(len(FAMILY_ORDER), 0)
    axis.set_xticks(
        np.arange(len(THRESHOLDS)) + 0.5,
        [f"≥{value}" for value in THRESHOLDS],
    )
    family_labels = (
        table.drop_duplicates("motif_family")
        .set_index("motif_family")["family_label"]
        .reindex(FAMILY_ORDER)
    )
    axis.set_yticks(
        np.arange(len(FAMILY_ORDER)) + 0.5,
        family_labels,
    )
    axis.set_xlabel("Minimum donors with a differential expansion")
    axis.set_title(
        "A  Observed recurrent loci",
        loc="left",
        fontsize=14,
        fontweight="bold",
    )
    axis.spines[["top", "right"]].set_visible(False)


def plot_p_matrix(axis: plt.Axes, table: pd.DataFrame) -> None:
    p_values = table["two_sided_permutation_p"].clip(lower=1e-5)
    scores = -np.log10(p_values)
    normalization = Normalize(vmin=0, vmax=max(2.5, float(scores.max())))
    color_map = colormaps["magma"]

    for row_index, family in enumerate(FAMILY_ORDER):
        for column_index, threshold in enumerate(THRESHOLDS):
            record = table.loc[
                table["motif_family"].eq(family)
                & table["min_recurrent_donors"].eq(threshold)
            ].iloc[0]
            p_value = float(record.two_sided_permutation_p)
            score = -np.log10(max(p_value, 1e-5))
            facecolor = color_map(normalization(score))
            axis.add_patch(
                Rectangle(
                    (column_index, row_index),
                    1,
                    1,
                    facecolor=facecolor,
                    edgecolor="white",
                    linewidth=1.5,
                )
            )
            marker = "†" if p_value < 0.05 else ""
            text_color = "#202020" if score > 1.65 else "white"
            axis.text(
                column_index + 0.5,
                row_index + 0.5,
                f"P={format_p(p_value)}{marker}",
                ha="center",
                va="center",
                fontsize=9.2,
                fontweight="bold" if marker else "normal",
                color=text_color,
            )

    axis.set_xlim(0, len(THRESHOLDS))
    axis.set_ylim(len(FAMILY_ORDER), 0)
    axis.set_xticks(
        np.arange(len(THRESHOLDS)) + 0.5,
        [f"≥{value}" for value in THRESHOLDS],
    )
    axis.set_yticks([])
    axis.set_xlabel("Minimum donors with a differential expansion")
    axis.set_title(
        "B  Within-donor permutation tests",
        loc="left",
        fontsize=14,
        fontweight="bold",
    )
    axis.spines[["top", "right", "left"]].set_visible(False)

    scalar_mappable = plt.cm.ScalarMappable(
        norm=normalization,
        cmap=color_map,
    )
    colorbar = plt.colorbar(
        scalar_mappable,
        ax=axis,
        fraction=0.045,
        pad=0.03,
    )
    colorbar.set_label("−log₁₀(two-sided permutation P)")


def main() -> None:
    args = parse_args()
    table = pd.read_csv(args.input, sep="\t")
    expected_rows = len(FAMILY_ORDER) * len(THRESHOLDS)
    if len(table) != expected_rows:
        raise ValueError(f"Expected {expected_rows} result rows; found {len(table)}")

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 13,
            "axes.labelsize": 10.5,
            "pdf.fonttype": 42,
        }
    )
    figure, axes = plt.subplots(
        1,
        2,
        figsize=(13.5, 6.2),
        gridspec_kw={"width_ratios": [1.1, 1]},
    )
    plot_count_matrix(axes[0], table)
    plot_p_matrix(axes[1], table)
    figure.legend(
        handles=[
            Patch(facecolor=NEURON_COLOR, label="Left: Neuron"),
            Patch(facecolor=OLIGO_COLOR, label="Right: Oligodendrocyte"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.31, 0.885),
        ncol=2,
        frameon=False,
    )
    figure.subplots_adjust(
        left=0.19,
        right=0.96,
        bottom=0.10,
        top=0.78,
        wspace=0.13,
    )

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        output = args.output_prefix.with_suffix(f".{suffix}")
        figure.savefig(
            output,
            dpi=400 if suffix == "png" else None,
            bbox_inches="tight",
        )
        print(f"Wrote {output}")
    plt.close(figure)


if __name__ == "__main__":
    main()
