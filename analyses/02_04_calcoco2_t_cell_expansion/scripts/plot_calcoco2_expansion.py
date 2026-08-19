#!/usr/bin/env python3
"""Plot the CALCOCO2 flank-corrected repeat-length distributions."""

from __future__ import annotations

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
INPUT = HERE / "data" / "calcoco2_per_read_lengths.tsv"
OUT = HERE / "results" / "figure_3c_calcoco2_expansion"
COLORS = {"Lymphocytes": "#D55E00", "T-cell": "#CC6677", "NK-cell": "#DDCC77", "B-cell": "#44AA99", "Monocyte": "#0072B2"}
MARKERS = {"SniffCell": "o", "FANS": "^"}
CELL_ORDER = {"Lymphocytes": 1, "T-cell": 2, "NK-cell": 3, "B-cell": 4, "Monocyte": 5}


def normalize_cell(source: str, cell: str) -> str:
    if source == "SniffCell":
        return "Monocyte" if cell == "Monocyte" else "Lymphocytes"
    return cell


with INPUT.open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
if not rows:
    raise SystemExit("No CALCOCO2 rows")
for row in rows:
    row["cell_type"] = normalize_cell(row["source"], row["cell_type"])
    row["length"] = float(row["flank_corrected_length_bp"])

groups: dict[tuple[str, str, str], list[float]] = defaultdict(list)
for row in rows:
    groups[(row["donor"], row["source"], row["cell_type"])].append(row["length"])
series = sorted(groups, key=lambda key: (key[0], 0 if key[1] == "SniffCell" else 1, CELL_ORDER[key[2]]))

fig, ax = plt.subplots(figsize=(11.6, 5.0))
rng = random.Random(48857948)
for x, key in enumerate(series):
    donor, source, cell = key
    for value in groups[key]:
        ax.scatter(x + rng.uniform(-0.16, 0.16), value / 5.0, s=28, alpha=0.78,
                   color=COLORS[cell], marker=MARKERS[source], edgecolor="#333333", linewidth=0.25)
for x in range(len(series) - 1):
    left, right = series[x], series[x + 1]
    ax.axvline(x + 0.5, color="#666666" if left[0] != right[0] else "#BDBDBD", lw=1.4 if left[0] != right[0] else 0.7)
ax.set_xlim(-1.25, len(series) - 0.2)
ax.set_xticks(range(len(series)), ["" for _ in series])
ax.set_ylabel("Read length minus 500 bp flanks / 5")
ax.grid(axis="y", color="#D9D9D9", lw=0.7)
ax.spines[["top", "right"]].set_visible(False)

labels = {"T-cell": "T", "NK-cell": "NK", "B-cell": "B", "Monocyte": "Mono", "Lymphocytes": "Lympho\ncytes"}
for x, (donor, source, cell) in enumerate(series):
    ax.text(x, -0.18, labels[cell], transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=8)
for donor in sorted({key[0] for key in series}):
    idx = [i for i, key in enumerate(series) if key[0] == donor]
    ax.text(sum(idx) / len(idx), -0.08, donor, transform=ax.get_xaxis_transform(), ha="center", va="top", fontweight="bold")
for donor, source in sorted({(key[0], key[1]) for key in series}):
    idx = [i for i, key in enumerate(series) if key[:2] == (donor, source)]
    ax.text(sum(idx) / len(idx), -0.13, "FACS" if source == "FANS" else source,
            transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=9)
ax.text(-0.95, -0.08, "donor", transform=ax.get_xaxis_transform(), ha="right", va="top", fontweight="bold", color="#555555")
ax.text(-0.95, -0.13, "source", transform=ax.get_xaxis_transform(), ha="right", va="top", fontweight="bold", color="#555555")
ax.text(-0.95, -0.18, "cell", transform=ax.get_xaxis_transform(), ha="right", va="top", fontweight="bold", color="#555555")

cell_handles = [Line2D([0], [0], marker="o", linestyle="", color=color, label=cell) for cell, color in COLORS.items()]
source_handles = [Line2D([0], [0], marker=marker, linestyle="", color="#777777", label="FACS" if source == "FANS" else source) for source, marker in MARKERS.items()]
legend = fig.legend(cell_handles, [h.get_label() for h in cell_handles], title="cell type", loc="upper left", bbox_to_anchor=(0.82, 0.82))
fig.add_artist(legend)
fig.legend(source_handles, [h.get_label() for h in source_handles], title="source", loc="upper left", bbox_to_anchor=(0.82, 0.49))
fig.subplots_adjust(bottom=0.29, right=0.78, left=0.10, top=0.96)
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT.with_suffix(".png"), dpi=300)
fig.savefig(OUT.with_suffix(".pdf"))
plt.close(fig)

with OUT.with_name(OUT.name + "_summary.tsv").open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["donor", "source", "cell_type", "read_count", "minimum_bp", "median_bp", "maximum_bp"])
    for key in series:
        values = groups[key]
        writer.writerow([*key, len(values), min(values), statistics.median(values), max(values)])
