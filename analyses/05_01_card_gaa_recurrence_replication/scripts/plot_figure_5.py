#!/usr/bin/env python3
"""Compose manuscript Figure 5 from public aggregate result tables."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm, to_rgb
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "data"
FIGURES = ROOT / "results"

RECURRENCE = TABLES / "card_gaa_motif_family_recurrence_count_matrix.tsv"
PAIRWISE = TABLES / "cross_cohort_neuron_pairwise_replication.tsv"
THREEWAY = TABLES / "cross_cohort_threeway_gaa_loci.tsv"

OUT_STEM = "figure_5_reproduced"

NEURON_COLOR = "#E67E22"
OLIGO_COLOR = "#1B9E77"
FAMILY_ORDER = [
    "pure_GAA",
    "GAAA_family",
    "GGAA_family",
    "other_4bp_GAA",
    "other_5_10bp_GAA",
]
THRESHOLDS = [1, 2, 3, 5]
COHORTS = ["CARD", "SMaHT", "4-brain benchmark"]
COHORT_LABELS = {
    "CARD": "CARD\n(155 donors)",
    "SMaHT": "SMaHT\n(13 donors)",
    "4-brain benchmark": "4-brain\nbenchmark",
}


def lighten(root_color: str, count: int, maximum: int) -> tuple[float, float, float]:
    strength = (
        0.04
        if count == 0
        else 0.20 + 0.80 * np.log1p(count) / np.log1p(maximum)
    )
    root = np.array(to_rgb(root_color))
    return tuple(1 - (1 - root) * strength)


def parse_pair(value: str) -> tuple[int, int]:
    neuron, oligo = value.split("/")
    return int(neuron.strip()), int(oligo.strip())


def panel_a(ax: plt.Axes, recurrence: pd.DataFrame) -> None:
    recurrence = recurrence.set_index("motif_family").loc[FAMILY_ORDER]
    all_counts = []
    parsed: dict[tuple[str, int], tuple[int, int]] = {}
    for family in FAMILY_ORDER:
        for threshold in THRESHOLDS:
            counts = parse_pair(recurrence.loc[family, f">={threshold}"])
            parsed[(family, threshold)] = counts
            all_counts.extend(counts)
    maximum = max(all_counts)

    group_centers = np.arange(len(THRESHOLDS), dtype=float) * 1.24
    half_width = 0.41
    for row_index, family in enumerate(FAMILY_ORDER):
        for column_index, threshold in enumerate(THRESHOLDS):
            center = group_centers[column_index]
            neuron, oligo = parsed[(family, threshold)]
            for cell_type, count, x_start, color in [
                ("Neuron", neuron, center - half_width, NEURON_COLOR),
                ("Oligodendrocyte", oligo, center, OLIGO_COLOR),
            ]:
                ax.add_patch(
                    Rectangle(
                        (x_start, row_index - 0.5),
                        half_width,
                        1,
                        facecolor=lighten(color, count, maximum),
                        edgecolor="white",
                        linewidth=1.5,
                    )
                )
                strength = (
                    0
                    if count == 0
                    else 0.20 + 0.80 * np.log1p(count) / np.log1p(maximum)
                )
                ax.text(
                    x_start + half_width / 2,
                    row_index,
                    str(count),
                    ha="center",
                    va="center",
                    fontsize=11,
                    fontweight="bold",
                    color="white" if strength > 0.64 else "#222222",
                )

    ax.set_xlim(group_centers[0] - half_width, group_centers[-1] + half_width)
    ax.set_ylim(len(FAMILY_ORDER) - 0.5, -0.5)
    ax.set_xticks(group_centers, [f"≥{value}" for value in THRESHOLDS])
    ax.set_yticks(np.arange(len(FAMILY_ORDER)), recurrence["family_label"])
    ax.set_xlabel("Minimum donors with an expansion")
    ax.set_title("A  CARD GAA-family recurrence", loc="left", fontweight="bold")
    ax.legend(
        handles=[
            Patch(facecolor=NEURON_COLOR, label="Left: Neuron"),
            Patch(facecolor=OLIGO_COLOR, label="Right: Oligodendrocyte"),
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, 1.01),
        frameon=False,
        ncol=2,
        fontsize=9,
    )
    ax.tick_params(axis="both", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)


def panel_b(ax: plt.Axes, pairwise: pd.DataFrame) -> None:
    pairwise = pairwise.loc[pairwise["scope"].eq("any_neuron")]
    values = np.full((3, 3), np.nan)
    annotations: dict[tuple[int, int], str] = {}
    for i, source in enumerate(COHORTS):
        for j, target in enumerate(COHORTS):
            if source == target:
                continue
            row = pairwise.loc[
                pairwise["source_cohort"].eq(source)
                & pairwise["target_cohort"].eq(target)
            ].iloc[0]
            values[i, j] = row["gaa_replication_fraction"]
            annotations[(i, j)] = (
                f"{int(row['replicated_gaa_element_loci'])}/"
                f"{int(row['source_gaa_element_loci'])}\n"
                f"{100 * row['gaa_replication_fraction']:.1f}%"
            )

    cmap = LinearSegmentedColormap.from_list(
        "replication_orange", ["#fbf7f3", NEURON_COLOR]
    )
    ax.imshow(values, vmin=0, vmax=np.nanmax(values), cmap=cmap)
    for i in range(3):
        for j in range(3):
            if i == j:
                ax.add_patch(
                    Rectangle(
                        (j - 0.5, i - 0.5),
                        1,
                        1,
                        facecolor="#e4e4e4",
                        edgecolor="white",
                        linewidth=1.5,
                    )
                )
                ax.text(j, i, "—", ha="center", va="center", color="#777777")
            else:
                ax.text(
                    j,
                    i,
                    annotations[(i, j)],
                    ha="center",
                    va="center",
                    fontsize=9,
                    fontweight="bold",
                    color="#222222",
                )

    ax.set_xticks(range(3), [COHORT_LABELS[c] for c in COHORTS], fontsize=8.5)
    ax.set_yticks(range(3), [COHORT_LABELS[c] for c in COHORTS], fontsize=8.5)
    ax.set_xlabel("Replication cohort")
    ax.set_ylabel("Source cohort")
    ax.set_title(
        "B  Cross-cohort replication of GAA-element Neuron loci",
        loc="left",
        fontweight="bold",
    )
    ax.set_xticks(np.arange(-0.5, 3, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, 3, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)


def locus_label(row: pd.Series) -> str:
    gene = row["gene"] if isinstance(row["gene"], str) else ""
    prefix = f"{gene}\n" if gene else ""
    return f"{prefix}{row['site_key']}\n({row['motif']})"


def panel_c(ax: plt.Axes, loci: pd.DataFrame) -> None:
    # Match the requested wide layout: six loci as columns, cohorts as rows.
    row_order = ["4-brain benchmark", "SMaHT", "CARD"]
    prefixes = {
        "CARD": "card",
        "SMaHT": "smaht",
        "4-brain benchmark": "benchmark",
    }
    matrix = np.array(
        [
            [row[f"{prefixes[cohort]}_neuron_donors"] for _, row in loci.iterrows()]
            for cohort in row_order
        ],
        dtype=float,
    )
    image = ax.imshow(
        matrix,
        aspect="auto",
        cmap="Oranges",
        norm=LogNorm(vmin=1, vmax=max(2, matrix.max())),
    )
    for i, cohort in enumerate(row_order):
        prefix = prefixes[cohort]
        for j, (_, row) in enumerate(loci.iterrows()):
            n_neuron = int(row[f"{prefix}_neuron_donors"])
            n_oligo = int(row[f"{prefix}_oligodendrocyte_donors"])
            ax.text(
                j,
                i,
                f"N {n_neuron} / O {n_oligo}",
                ha="center",
                va="center",
                fontsize=14,
                fontweight="bold",
                color="white" if n_neuron >= 4 else "#222222",
            )

    ax.set_xticks(
        range(len(loci)),
        [locus_label(row) for _, row in loci.iterrows()],
        rotation=0,
        ha="center",
        va="top",
        rotation_mode="default",
        fontsize=10.5,
    )
    ax.set_yticks(
        range(3), [COHORT_LABELS[c] for c in row_order], fontsize=13
    )
    ax.set_title(
        "C  GAA-element Neuron loci found in all three cohorts",
        loc="left",
        fontweight="bold",
        fontsize=17,
    )
    ax.set_xlabel("")
    ax.text(
        1.0,
        1.015,
        "N/O = donors with Neuron/Oligodendrocyte expansion",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=12,
        color="#555555",
    )
    ax.set_xticks(np.arange(-0.5, len(loci), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, 3, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.tick_params(axis="x", pad=5)
    for spine in ax.spines.values():
        spine.set_visible(False)
    image.set_clim(1, max(2, matrix.max()))


def reconcile(
    recurrence: pd.DataFrame, pairwise: pd.DataFrame, loci: pd.DataFrame
) -> None:
    checks = [
        ("Panel A motif families", len(recurrence), 5),
        ("Panel B directional comparisons", len(pairwise.loc[pairwise["scope"].eq("any_neuron")]), 6),
        ("Panel C three-way GAA loci", len(loci), 6),
        (
            "Panel C SH3RF3 CARD Neuron donors",
            int(
                loci.loc[loci["site_key"].eq("chr2:109199301-109199876"), "card_neuron_donors"].iloc[0]
            ),
            17,
        ),
        (
            "Panel B CARD to SMaHT replicated GAA loci",
            int(
                pairwise.loc[
                    pairwise["scope"].eq("any_neuron")
                    & pairwise["source_cohort"].eq("CARD")
                    & pairwise["target_cohort"].eq("SMaHT"),
                    "replicated_gaa_element_loci",
                ].iloc[0]
            ),
            11,
        ),
    ]
    result = pd.DataFrame(checks, columns=["check", "observed", "expected"])
    result["status"] = np.where(
        result["observed"].eq(result["expected"]), "PASS", "FAIL"
    )
    result.to_csv(
        FIGURES / "figure_5_reconciliation.tsv", sep="\t", index=False
    )
    if result["status"].ne("PASS").any():
        raise RuntimeError("Figure 5 reconciliation failed")


def main() -> None:
    recurrence = pd.read_csv(RECURRENCE, sep="\t")
    pairwise = pd.read_csv(PAIRWISE, sep="\t")
    loci = pd.read_csv(THREEWAY, sep="\t", keep_default_na=False)
    reconcile(recurrence, pairwise, loci)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
        }
    )
    fig = plt.figure(figsize=(15.5, 11.5), constrained_layout=False)
    grid = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.35, 0.95],
        height_ratios=[1.0, 0.92],
        left=0.08,
        right=0.985,
        top=0.955,
        bottom=0.17,
        wspace=0.30,
        hspace=0.36,
    )
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, :])

    panel_a(ax_a, recurrence)
    panel_b(ax_b, pairwise)
    panel_c(ax_c, loci)

    for suffix in ("png", "pdf"):
        fig.savefig(
            FIGURES / f"{OUT_STEM}.{suffix}",
            dpi=400 if suffix == "png" else None,
            bbox_inches="tight",
        )
    plt.close(fig)
    print(FIGURES / f"{OUT_STEM}.png")
    print(FIGURES / f"{OUT_STEM}.pdf")


if __name__ == "__main__":
    main()
