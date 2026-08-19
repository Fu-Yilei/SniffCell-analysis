#!/usr/bin/env python3
"""Plot Figure 4e donor-by-region counts of globally specific TR loci."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
REGIONS = ["Cerebellum", "Frontal Lobe", "Hippocampus", "Temporal Lobe"]


def main() -> None:
    data = pd.read_csv(DATA / "globally_specific_locus_counts_by_donor_region.tsv", sep="\t")
    donors = sorted(data["donor"].unique())
    fig, axes = plt.subplots(2, 1, figsize=(13, 5.8), sharex=True, constrained_layout=True)
    for ax, column, title, cmap in [
        (axes[0], "neuron_specific_loci", "Neuron-specific loci", "Oranges"),
        (axes[1], "oligodendrocyte_specific_loci", "Oligodendrocyte-specific loci", "BuGn"),
    ]:
        matrix = data.pivot(index="brain_region", columns="donor", values=column).reindex(index=REGIONS, columns=donors)
        assayed = data.pivot(index="brain_region", columns="donor", values="has_sample").reindex(index=REGIONS, columns=donors)
        shown = matrix.to_numpy(float); shown[~assayed.to_numpy(bool)] = np.nan
        cm = plt.get_cmap(cmap).copy(); cm.set_bad("#D9D9D9")
        ax.imshow(shown, aspect="auto", cmap=cm, vmin=0)
        ax.set_yticks(range(len(REGIONS)), REGIONS)
        ax.set_title(title, loc="left", fontweight="bold")
        for i in range(len(REGIONS)):
            for j in range(len(donors)):
                if assayed.iloc[i, j]: ax.text(j, i, str(int(matrix.iloc[i, j])), ha="center", va="center", fontsize=7)
    axes[1].set_xticks(range(len(donors)), donors, rotation=90)
    fig.savefig(RESULTS / "figure_4e_reproduced.png", dpi=300, bbox_inches="tight")
    fig.savefig(RESULTS / "figure_4e_reproduced.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
