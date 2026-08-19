#!/usr/bin/env python3
"""Plot Figure 4d and Supplementary Figure 16 from aggregate overlap tables."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
REGIONS = ["Frontal Lobe", "Temporal Lobe", "Hippocampus", "Cerebellum"]
GROUPS = [("All", "All specific loci", "#5F5F5F"), ("Neuron", "Neuron-specific loci", "#E67E22"), ("Oligodendrocyte", "Oligodendrocyte-specific loci", "#1B9E77")]


def main() -> None:
    pair = pd.read_csv(DATA / "figure6_region_agreement_specific_loci_source.tsv", sep="\t")
    heat = pd.read_csv(DATA / "figure6_region_agreement_heatmap_source.tsv", sep="\t")
    all_rows = pair[pair["group"].eq("All")].copy()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    colors = np.where(all_rows["region_pair"].str.contains("Cerebellum"), "#A0A0A0", "#4C78A8")
    axes[0].barh(all_rows["region_pair"], all_rows["jaccard_percent"], color=colors)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Shared / union (%)")
    axes[0].set_title("Pairwise regional overlap")
    for y, row in enumerate(all_rows.itertuples(index=False)):
        axes[0].text(row.jaccard_percent + .5, y, f"{row.n_shared}/{row.n_union}", va="center", fontsize=8)
    lobe = all_rows.iloc[0]
    axes[1].bar(["Any cerebral lobe", "Also in cerebellum"], [lobe.n_seen_in_any_lobe, lobe.n_lobe_seen_also_cerebellum], color=["#4C78A8", "#A0A0A0"])
    axes[1].set_ylabel("Specific loci")
    axes[1].set_title(f"Cerebellum capture: {lobe.label} ({lobe.percent_lobe_seen_also_cerebellum:.1f}%)")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    RESULTS.mkdir(exist_ok=True)
    fig.savefig(RESULTS / "figure_4d_reproduced.png", dpi=300, bbox_inches="tight")
    fig.savefig(RESULTS / "figure_4d_reproduced.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), constrained_layout=True)
    image = None
    for ax, (group, title, _) in zip(axes, GROUPS):
        key = "Oligo" if group == "Oligodendrocyte" else group
        part = heat[heat["group"].eq(key)]
        matrix = part.pivot(index="region_a", columns="region_b", values="jaccard_percent").reindex(index=REGIONS, columns=REGIONS)
        shared = part.pivot(index="region_a", columns="region_b", values="n_shared").reindex(index=REGIONS, columns=REGIONS)
        union = part.pivot(index="region_a", columns="region_b", values="n_union").reindex(index=REGIONS, columns=REGIONS)
        display = matrix.to_numpy(float)
        display[np.triu_indices_from(display)] = np.nan
        cmap = plt.get_cmap("YlGnBu").copy(); cmap.set_bad("white")
        image = ax.imshow(display, cmap=cmap, vmin=0, vmax=35)
        ax.set_xticks(range(4), [r.replace(" Lobe", "") for r in REGIONS], rotation=35, ha="right")
        ax.set_yticks(range(4), [r.replace(" Lobe", "") for r in REGIONS])
        ax.set_title(title)
        for i in range(4):
            for j in range(i):
                ax.text(j, i, f"{matrix.iloc[i,j]:.0f}%\n{int(shared.iloc[i,j])}/{int(union.iloc[i,j])}", ha="center", va="center", fontsize=8)
    fig.colorbar(image, ax=axes, label="Shared / union (%)", shrink=.8)
    fig.savefig(RESULTS / "supplementary_figure_16_reproduced.png", dpi=300, bbox_inches="tight")
    fig.savefig(RESULTS / "supplementary_figure_16_reproduced.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
