#!/usr/bin/env python3
"""Plot aggregate TR motif length and sequence composition."""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
INPUT = HERE / "data" / "tr_motif_composition.tsv"
OUT = HERE / "results" / "tr_motif_composition"
COLORS = {"PBMC": "#D95F02", "Brain": "#0072B2"}
ORDER = {
    "Motif length (bp)": ["2", "3", "4", "5", "6–10", ">10"],
    "Sequence composition": ["AT-rich", "Mixed", "GC-rich"],
}

with INPUT.open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
for tissue in ("PBMC", "Brain"):
    expected = 10 if tissue == "PBMC" else 105
    for dimension in ORDER:
        total = sum(int(row["count"]) for row in rows if row["tissue"] == tissue and row["dimension"] == dimension)
        if total != expected:
            raise SystemExit(f"{tissue} {dimension} sums to {total}, expected {expected}")

fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw={"wspace": 0.34})
for ax, (dimension, categories) in zip(axes, ORDER.items()):
    y = list(range(len(categories)))
    for tissue, offset in (("PBMC", -0.18), ("Brain", 0.18)):
        lookup = {(row["tissue"], row["category"]): row for row in rows if row["dimension"] == dimension}
        counts = [int(lookup[(tissue, category)]["count"]) for category in categories]
        percents = [float(lookup[(tissue, category)]["percent"]) for category in categories]
        bars = ax.barh([value + offset for value in y], counts, height=0.34, color=COLORS[tissue], label=f"{tissue} (n={10 if tissue == 'PBMC' else 105})")
        for bar, count, percent in zip(bars, counts, percents):
            ax.text(max(count + 0.7, 0.7), bar.get_y() + bar.get_height()/2, f"{count} ({percent:.1f}%)", va="center", fontsize=10)
    ax.set_yticks(y, categories)
    ax.invert_yaxis()
    ax.set_xlabel("TR loci (count)")
    ax.set_title("Motif length" if dimension.startswith("Motif") else "Sequence composition", fontsize=15)
    ax.grid(axis="x", color="#E6E6E6")
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
axes[0].set_xlim(0, 70)
axes[1].set_xlim(0, 100)
axes[0].legend(frameon=False, loc="upper right")
fig.subplots_adjust(left=0.09, right=0.985, bottom=0.16, top=0.90, wspace=0.34)
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT.with_suffix(".png"), dpi=300)
fig.savefig(OUT.with_suffix(".pdf"))
plt.close(fig)
