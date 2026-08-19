#!/usr/bin/env python3
"""Recompute cohort N50 summaries and plot Supplementary Figure 17."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/read_n50_all_samples.tsv"
RESULTS = ROOT / "results"


def main() -> None:
    samples = pd.read_csv(INPUT, sep="\t")
    required = {"cohort", "sample_id", "n50_bp"}
    missing = required.difference(samples.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    summary = (
        samples.groupby("cohort", sort=False)["n50_bp"]
        .agg(
            sample_count="size",
            median_n50_bp="median",
            mean_n50_bp="mean",
            min_n50_bp="min",
            max_n50_bp="max",
        )
        .reset_index()
    )
    summary.to_csv(RESULTS / "read_n50_cohort_summary.tsv", sep="\t", index=False)

    hbcc = samples.loc[samples["cohort"].eq("HBCC"), "n50_bp"]
    comparisons = []
    for cohort in ("brain_benchmark_bulk", "smaht_ont_brain"):
        other = samples.loc[samples["cohort"].eq(cohort), "n50_bp"]
        test = mannwhitneyu(hbcc, other, alternative="greater", method="auto")
        comparisons.append(
            {
                "comparison": f"HBCC > {cohort}",
                "hbcc_n": len(hbcc),
                "other_n": len(other),
                "hbcc_median_n50_bp": float(hbcc.median()),
                "other_median_n50_bp": float(other.median()),
                "median_difference_bp": float(hbcc.median() - other.median()),
                "hbcc_fraction_above_other_median": float(
                    (hbcc > other.median()).mean()
                ),
                "mann_whitney_u": float(test.statistic),
                "one_sided_pvalue": float(test.pvalue),
            }
        )
    pd.DataFrame(comparisons).to_csv(
        RESULTS / "read_n50_comparisons.tsv", sep="\t", index=False
    )

    order = ["brain_benchmark_bulk", "smaht_ont_brain", "HBCC"]
    labels = [
        "FANS validation\nbrain cohort",
        "multi-region\nbrain cohort",
        "large-scale\nbrain cohort",
    ]
    colors = ["#0072B2", "#009E73", "#D55E00"]
    values = [
        samples.loc[samples["cohort"].eq(cohort), "n50_bp"] / 1000
        for cohort in order
    ]

    rng = np.random.default_rng(17)
    figure, axis = plt.subplots(figsize=(8.6, 5.2))
    boxes = axis.boxplot(
        values,
        positions=range(1, 4),
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.5},
    )
    for patch, color in zip(boxes["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.35)
    for position, cohort_values, color in zip(range(1, 4), values, colors):
        jitter = rng.uniform(-0.16, 0.16, size=len(cohort_values))
        axis.scatter(
            position + jitter,
            cohort_values,
            s=22,
            color=color,
            alpha=0.75,
            edgecolors="white",
            linewidths=0.35,
        )
    axis.set_xticks(range(1, 4), labels)
    axis.set_ylabel("Read length N50 (kb)")
    axis.grid(axis="y", alpha=0.25)
    figure.tight_layout()
    for suffix in ("png", "pdf"):
        figure.savefig(
            RESULTS / f"supplementary_figure_17_reproduced.{suffix}",
            dpi=300 if suffix == "png" else None,
        )
    plt.close(figure)


if __name__ == "__main__":
    main()
