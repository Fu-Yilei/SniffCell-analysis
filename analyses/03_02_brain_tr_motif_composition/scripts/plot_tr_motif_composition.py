#!/usr/bin/env python3
"""Plot aggregate TR motif length and sequence composition."""

import csv
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"pdf.fonttype": 42, "svg.fonttype": "none"})

HERE = Path(__file__).resolve().parents[1]
INPUT = HERE / "data" / "tr_motif_composition.tsv"
OUT = HERE / "results" / "tr_motif_composition"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, default=INPUT)
parser.add_argument("--outdir", type=Path, default=OUT.parent)
args = parser.parse_args()
INPUT, OUT = args.input, args.outdir / OUT.name
COLORS = {"PBMC": "#D95F02", "Brain": "#0072B2"}
ORDER = {
    "Motif length (bp)": ["2", "3", "4", "5", "6–10", ">10"],
    "Sequence composition": ["AT-rich", "Mixed", "GC-rich"],
}

with INPUT.open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
totals = {}
for tissue in ("PBMC", "Brain"):
    declared = {int(row["tissue_total"]) for row in rows if row["tissue"] == tissue}
    if len(declared) != 1:
        raise SystemExit(f"Inconsistent denominators for {tissue}")
    expected = totals[tissue] = declared.pop()
    for dimension in ORDER:
        total = sum(int(row["count"]) for row in rows if row["tissue"] == tissue and row["dimension"] == dimension)
        if total != expected:
            raise SystemExit(f"{tissue} {dimension} sums to {total}, expected {expected}")
        subset = [row for row in rows if row["tissue"] == tissue and row["dimension"] == dimension]
        if len(subset) != len(ORDER[dimension]) or {r["category"] for r in subset} != set(ORDER[dimension]):
            raise SystemExit("Missing or duplicate motif categories")
        for row in subset:
            percent = 100 * int(row["count"])/expected if expected else 0
            if int(row["count"]) < 0 or abs(float(row["percent"])-percent) > 0.051:
                raise SystemExit("Invalid count or percent")

fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw={"wspace": 0.34})
for ax, (dimension, categories) in zip(axes, ORDER.items()):
    y = list(range(len(categories)))
    for tissue, offset in (("PBMC", -0.18), ("Brain", 0.18)):
        lookup = {(row["tissue"], row["category"]): row for row in rows if row["dimension"] == dimension}
        counts = [int(lookup[(tissue, category)]["count"]) for category in categories]
        percents = [float(lookup[(tissue, category)]["percent"]) for category in categories]
        bars = ax.barh([value + offset for value in y], counts, height=0.34, color=COLORS[tissue], label=f"{tissue} (n={totals[tissue]})")
        for bar, count, percent in zip(bars, counts, percents):
            ax.text(max(count + 0.7, 0.7), bar.get_y() + bar.get_height()/2, f"{count} ({percent:.1f}%)", va="center", fontsize=10)
    ax.set_yticks(y, categories)
    ax.invert_yaxis()
    ax.set_xlabel("TR candidates (donor–locus rows)")
    ax.set_ylabel(dimension)
    ax.grid(axis="x", color="#E6E6E6")
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
for ax, dimension in zip(axes, ORDER):
    maximum = max(int(r["count"]) for r in rows if r["dimension"] == dimension)
    ax.set_xlim(0, max(10, maximum * 1.35))
axes[0].legend(frameon=False, loc="upper right")
fig.subplots_adjust(left=0.09, right=0.985, bottom=0.16, top=0.90, wspace=0.34)
OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT.with_suffix(".png"), dpi=300)
fig.savefig(OUT.with_suffix(".pdf"), metadata={"CreationDate": None, "ModDate": None})
plt.close(fig)
