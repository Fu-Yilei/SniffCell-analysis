#!/usr/bin/env python3
"""Plot Figure 4b-c aggregate regional prevalence and permutation results."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
REGIONS = ["Frontal Lobe", "Temporal Lobe", "Hippocampus", "Cerebellum"]
COLORS = {"Neuron": "#E67E22", "Oligodendrocyte": "#1B9E77"}


def main() -> None:
    counts = pd.read_csv(DATA / "regional_counts.tsv", sep="\t")
    observed = pd.read_csv(DATA / "observed_regional_prevalence.tsv", sep="\t")
    tests = pd.read_csv(DATA / "regional_permutation_tests.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.7))

    for subset, marker in [("All exclusive loci", "o"), ("Recurrent exclusive loci", "s")]:
        part = counts[counts["subset"].eq(subset)]
        for cell in ["Neuron", "Oligodendrocyte"]:
            values = part[part["change_cell"].eq(cell)].set_index("brain_region").reindex(REGIONS)
            axes[0].plot(REGIONS, 100 * values["fraction"], marker=marker, color=COLORS[cell],
                         linestyle="-" if subset.startswith("All") else "--",
                         label=f"{cell}; {subset.replace(' exclusive loci', '')}")
    axes[0].set_ylabel("Loci observed in region (%)")
    axes[0].tick_params(axis="x", rotation=25)
    axes[0].set_title("Regional breadth of exclusive loci")
    axes[0].legend(frameon=False, fontsize=8)

    width = 0.36
    x = np.arange(len(REGIONS))
    for offset, cell in [(-width / 2, "Neuron"), (width / 2, "Oligodendrocyte")]:
        part = observed[observed["cell_type"].eq(cell)].set_index("brain_region").reindex(REGIONS)
        axes[1].bar(x + offset, part["donor_locus_prevalence_percent"], width, color=COLORS[cell], label=cell)
    axes[1].set_xticks(x, REGIONS, rotation=25)
    axes[1].set_ylabel("Donor-locus prevalence (%)")
    axes[1].set_title("Donor-normalized regional prevalence")
    axes[1].legend(frameon=False)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    RESULTS.mkdir(exist_ok=True)
    fig.savefig(RESULTS / "figure_4bc_reproduced.png", dpi=300, bbox_inches="tight")
    fig.savefig(RESULTS / "figure_4bc_reproduced.pdf", bbox_inches="tight")
    plt.close(fig)
    tests.to_csv(RESULTS / "regional_permutation_tests.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
