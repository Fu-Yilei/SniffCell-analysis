#!/usr/bin/env python3
"""Plot SGIP1 per-read repeat lengths across four brain donors."""

import csv
import random
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parents[1]
INPUT = HERE / "data" / "sgip1_per_read_lengths.tsv"
OUT = HERE / "results" / "supplementary_figure_15_sgip1"
DONORS = {"b1": "bcontrol1", "C2": "bcontrol2", "b2": "bmsa", "P2": "bparkinsons"}
DONOR_ORDER = ["b1", "C2", "b2", "P2"]
CELLS = ["Neuron", "Oligodendrocyte"]
SOURCES = ["SniffCell", "FANS"]
COLORS = {"Neuron": "#E67E22", "Oligodendrocyte": "#1B9E77"}
MARKERS = {"SniffCell": "o", "FANS": "^"}

with INPUT.open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
groups = defaultdict(list)
for row in rows:
    if row["donor"] in DONORS:
        groups[(row["donor"], row["source"], row["cell_type"])].append(float(row["flank_corrected_length_bp"]))
series = [(donor, source, cell) for donor in DONOR_ORDER for source in SOURCES for cell in CELLS]
if any(not groups[key] for key in series):
    raise SystemExit("Missing a donor/source/cell group")

fig, ax = plt.subplots(figsize=(11.5, 6.8))
rng = random.Random("SGIP1")
for x, key in enumerate(series):
    for value in groups[key]:
        ax.scatter(x + rng.uniform(-0.22, 0.22), value, s=34, alpha=0.75,
                   color=COLORS[key[2]], marker=MARKERS[key[1]], edgecolor="black", linewidth=0.2)
for x in range(len(series) - 1):
    if series[x][0] != series[x + 1][0]:
        ax.axvline(x + 0.5, color="#666666", lw=1.5)
ax.set_xticks(range(len(series)), ["" for _ in series])
ax.set_ylabel("TR sizes per read")
ax.grid(axis="y", color="#D9D9D9", lw=0.7)
ax.spines[["top", "right"]].set_visible(False)
for x, (donor, source, cell) in enumerate(series):
    ax.text(x, -0.08, "Oligo" if cell == "Oligodendrocyte" else cell,
            transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=8)
for donor in DONOR_ORDER:
    idx = [i for i, key in enumerate(series) if key[0] == donor]
    ax.text(sum(idx)/len(idx), -0.23, DONORS[donor], transform=ax.get_xaxis_transform(), ha="center", va="top", fontweight="bold")
    for source in SOURCES:
        sub = [i for i, key in enumerate(series) if key[:2] == (donor, source)]
        ax.text(sum(sub)/len(sub), -0.15, source, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=9)
for label, y in (("cell", -0.08), ("source", -0.15), ("donor", -0.23)):
    ax.text(-0.95, y, label, transform=ax.get_xaxis_transform(), ha="right", va="top", fontweight="bold")

cell_handles = [Line2D([0], [0], marker="o", linestyle="", color=color, label=cell) for cell, color in COLORS.items()]
source_handles = [Line2D([0], [0], marker=marker, linestyle="", color="#888888", label=source) for source, marker in MARKERS.items()]
legend = ax.legend(handles=cell_handles, title="cell type", loc="upper left", bbox_to_anchor=(1.02, 0.88))
ax.add_artist(legend)
ax.legend(handles=source_handles, title="source", loc="upper left", bbox_to_anchor=(1.02, 0.55))
fig.subplots_adjust(left=0.09, right=0.78, bottom=0.28, top=0.97)
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT.with_suffix(".png"), dpi=300)
fig.savefig(OUT.with_suffix(".pdf"))
plt.close(fig)

with OUT.with_name(OUT.name + "_summary.tsv").open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["donor", "source", "cell_type", "read_count", "minimum_bp", "median_bp", "maximum_bp"])
    for key in series:
        values = groups[key]
        writer.writerow([DONORS[key[0]], key[1], key[2], len(values), min(values), statistics.median(values), max(values)])
