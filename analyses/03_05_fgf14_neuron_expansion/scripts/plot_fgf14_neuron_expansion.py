#!/usr/bin/env python3
"""Plot FGF14 neuron-vs-oligodendrocyte quantile distributions."""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

HERE = Path(__file__).resolve().parents[1]
DATA = HERE / "data"
OUT = HERE / "results" / "figure_3d_fgf14_neuron_expansion"
SAMPLES = ["b1", "b4", "b2", "b3"]
DISPLAY = {"b1": "bcontrol1", "b4": "bcontrol2", "b2": "bmsa", "b3": "bparkinsons"}
MARKERS = {"b1": "o", "b4": "s", "b2": "^", "b3": "D"}
NEURON, OLIGO = "#E67E22", "#1B9E77"
THRESHOLD = 320.0


def read(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


reads = read(DATA / "fgf14_per_read.tsv")
pairs = read(DATA / "fgf14_quantile_pairs.tsv")
if "read_name" in reads[0] or "run_id" in reads[0]:
    raise SystemExit("Public per-read input contains an excluded identifier column")
for row in reads:
    expected = float(row["consensus_motif_count"]) + (float(row["read_allele_length_bp"]) - float(row["consensus_allele_length_bp"])) / 3.0
    if abs(expected - float(row["estimated_GAA_repeats"])) > 1e-8:
        raise SystemExit("Stored FGF14 estimate does not match its documented formula")


def smooth(values, upper, bins=100):
    edges = np.linspace(0, upper, bins + 1)
    density, _ = np.histogram(values, bins=edges, density=True)
    x = np.arange(-7, 8)
    kernel = np.exp(-0.5 * (x / 2.2) ** 2)
    kernel /= kernel.sum()
    return (edges[:-1] + edges[1:]) / 2, np.convolve(density, kernel, mode="same")


fig = plt.figure(figsize=(9.1, 9.5))
grid = fig.add_gridspec(2, 2, width_ratios=(5.2, 1.15), height_ratios=(1.15, 5.2), left=0.11, right=0.94, bottom=0.09, top=0.96, wspace=0.04, hspace=0.04)
top = fig.add_subplot(grid[0, 0])
ax = fig.add_subplot(grid[1, 0])
right = fig.add_subplot(grid[1, 1])
for sample in SAMPLES:
    subset = [row for row in pairs if row["sample"] == sample]
    ax.scatter([float(row["oligodendrocyte_estimated_GAA"]) for row in subset],
               [float(row["neuron_estimated_GAA"]) for row in subset],
               marker=MARKERS[sample], s=42, facecolor=NEURON, edgecolor=OLIGO,
               linewidth=1.15, alpha=0.72, label=DISPLAY[sample])
ax.plot([0, 400], [0, 400], color="#9CA3AF", lw=1, ls="--")
ax.axvline(THRESHOLD, color=OLIGO, lw=1.7, ls=(0, (4, 3)))
ax.axhline(THRESHOLD, color=NEURON, lw=1.7, ls=(0, (4, 3)))
ax.text(312, 333, "320 GAA pathogenic threshold", ha="right", va="bottom", color="#374151")
ax.set(xlim=(0, 400), ylim=(0, 500), xlabel="Oligodendrocyte per-read TR length (GAA repeats)", ylabel="Neuron per-read TR length (GAA repeats)")
ax.xaxis.label.set_color(OLIGO)
ax.yaxis.label.set_color(NEURON)
ax.tick_params(axis="x", colors=OLIGO)
ax.tick_params(axis="y", colors=NEURON)
ax.grid(color="#E5E7EB", lw=0.7)
ax.spines[["top", "right"]].set_visible(False)

oligo = np.array([float(row["estimated_GAA_repeats"]) for row in reads if row["cell_type"] == "Oligodendrocyte"])
neuron = np.array([float(row["estimated_GAA_repeats"]) for row in reads if row["cell_type"] == "Neuron"])
x, d = smooth(oligo, 400)
top.fill_between(x, d, color=OLIGO, alpha=0.65)
top.plot(x, d, color=OLIGO)
top.axvline(THRESHOLD, color=OLIGO, ls=(0, (4, 3)))
top.set_xlim(0, 400)
top.axis("off")
y, d = smooth(neuron, 500)
right.fill_betweenx(y, 0, d, color=NEURON, alpha=0.65)
right.plot(d, y, color=NEURON)
right.axhline(THRESHOLD, color=NEURON, ls=(0, (4, 3)))
right.set_ylim(0, 500)
right.axis("off")

sample_handles = [Line2D([0], [0], marker=MARKERS[s], linestyle="", markerfacecolor=NEURON, markeredgecolor=OLIGO, label=DISPLAY[s]) for s in SAMPLES]
cell_handles = [Line2D([0], [0], color=NEURON, lw=4, label="Neuron (y / point fill)"), Line2D([0], [0], color=OLIGO, lw=4, label="Oligodendrocyte (x / point edge)")]
legend = ax.legend(handles=sample_handles, title="Brain sample", frameon=False, loc="upper left", ncol=2)
ax.add_artist(legend)
ax.legend(handles=cell_handles, frameon=False, loc="lower right")
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT.with_suffix(".png"), dpi=240, bbox_inches="tight")
fig.savefig(OUT.with_suffix(".pdf"), bbox_inches="tight")
plt.close(fig)
