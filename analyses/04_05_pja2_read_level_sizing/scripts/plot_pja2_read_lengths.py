#!/usr/bin/env python3
"""Replot the Figure 4f PJA2 locus from sanitized per-read lengths."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "data" / "pja2_read_lengths.tsv.gz"
COLORS = {"Blood": "#111111", "Neuron": "#9C6B62", "Oligodendrocyte": "#2CBBC3"}
MARKERS = {"Blood-ONT": "^", "Blood-PB": "v", "SniffCell-ONT": "o", "SniffCell-PB": "s"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-prefix", type=Path, default=ROOT / "results" / "figure_4f_reproduced")
    args = parser.parse_args()
    data = pd.read_csv(args.input, sep="\t")
    required = {"donor", "tissue_name", "source", "cell_type", "read_length"}
    if missing := required - set(data.columns):
        raise ValueError(f"Missing columns: {sorted(missing)}")

    donors = sorted(data["donor"].unique())
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharey=True)
    rng = np.random.default_rng(20260818)
    for ax, donor in zip(axes.ravel(), donors):
        part = data[data["donor"].eq(donor)].copy()
        part["track"] = part["tissue_name"] + " | " + part["cell_type"] + " | " + part["source"]
        tracks = part[["track", "tissue_name", "cell_type", "source"]].drop_duplicates().sort_values(["tissue_name", "cell_type", "source"])
        order = tracks["track"].tolist()
        for xx, track in enumerate(order):
            rows = part[part["track"].eq(track)]
            source = rows["source"].iloc[0]
            cell = rows["cell_type"].iloc[0]
            jitter = rng.uniform(-0.18, 0.18, len(rows))
            ax.scatter(np.full(len(rows), xx) + jitter, rows["read_length"], s=13,
                       color=COLORS[cell], marker=MARKERS[source], alpha=.75, linewidths=0)
        ax.set_xticks(range(len(order)), order, rotation=90, fontsize=6)
        ax.set_title(donor, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="#EAEAEA", linewidth=.6)
    for ax in axes[:, 0]:
        ax.set_ylabel("Repeat length (bp)")
    fig.tight_layout()
    args.output_prefix.parent.mkdir(exist_ok=True)
    fig.savefig(args.output_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(args.output_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    (data.groupby(["donor", "tissue_name", "source", "cell_type"], as_index=False)
         .agg(n_reads=("read_length", "size"), min_read_length=("read_length", "min"), median_read_length=("read_length", "median"), max_read_length=("read_length", "max"))
         .to_csv(ROOT / "results" / "pja2_read_length_summary.tsv", sep="\t", index=False))


if __name__ == "__main__":
    main()
