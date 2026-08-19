#!/usr/bin/env python3
"""Reproduce the aggregate Figure 4a burden and recurrence panels."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
COLORS = {"Neuron": "#E67E22", "Oligodendrocyte": "#1B9E77"}


def main() -> None:
    directional = pd.read_csv(DATA / "directional_recurrence_summary.tsv", sep="\t")
    exclusive = pd.read_csv(DATA / "exclusive_recurrence_summary.tsv", sep="\t")
    motifs = pd.read_csv(DATA / "recurrent_exclusive_aag_like_summary.tsv", sep="\t")
    motifs = motifs.set_index("change_cell").loc[["Neuron", "Oligodendrocyte"]].reset_index()
    table = motifs[["aag_like_loci", "other_loci"]].to_numpy(dtype=int)
    odds_ratio, p_value = fisher_exact(table)

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.2))
    for ax, frame, title in zip(
        axes[:2],
        [directional, exclusive],
        ["All directional events", "Cell-type-exclusive loci"],
    ):
        x = np.arange(len(frame))
        ax.bar(x, frame["nonrecurrent_events"], color="#D9D9D9", label="One donor")
        ax.bar(
            x,
            frame["recurrent_events"],
            bottom=frame["nonrecurrent_events"],
            color=[COLORS[v] for v in frame["change_cell"]],
            label="Recurrent",
        )
        for xx, total, recurrent in zip(x, frame["total_events"], frame["recurrent_events"]):
            ax.text(xx, total + 1.2, f"{int(total)} total\n{int(recurrent)} recurrent", ha="center", fontsize=9)
        ax.set_xticks(x, frame["change_cell"])
        ax.set_ylabel("TR directional events")
        ax.set_title(title)
        ax.spines[["top", "right"]].set_visible(False)

    x = np.arange(len(motifs))
    axes[2].bar(x, 100 * motifs["aag_like_fraction"], color=[COLORS[v] for v in motifs["change_cell"]])
    for xx, row in zip(x, motifs.itertuples(index=False)):
        axes[2].text(xx, 100 * row.aag_like_fraction + 2, f"{row.aag_like_loci}/{row.total_loci}", ha="center")
    axes[2].set_xticks(x, motifs["change_cell"])
    axes[2].set_ylabel("AAG-like recurrent exclusive loci (%)")
    axes[2].set_ylim(0, 52)
    axes[2].set_title(f"Motif enrichment\nFisher P={p_value:.3g}")
    axes[2].spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    RESULTS.mkdir(exist_ok=True)
    fig.savefig(RESULTS / "figure_4a_reproduced.png", dpi=300, bbox_inches="tight")
    fig.savefig(RESULTS / "figure_4a_reproduced.pdf", bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(
        [{"comparison": "Neuron vs Oligodendrocyte AAG-like recurrent exclusive loci", "odds_ratio": odds_ratio, "fisher_exact_two_sided_p": p_value}]
    ).to_csv(RESULTS / "motif_fisher_test.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
