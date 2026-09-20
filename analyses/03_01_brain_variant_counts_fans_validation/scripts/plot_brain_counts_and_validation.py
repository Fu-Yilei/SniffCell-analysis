#!/usr/bin/env python3
"""Plot aggregate brain candidate counts and paired-FANS validation."""

import csv
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"pdf.fonttype": 42, "svg.fonttype": "none"})

HERE = Path(__file__).resolve().parents[1]
DATA = HERE / "data"
RESULTS = HERE / "results"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--data-dir", type=Path, default=DATA)
parser.add_argument("--outdir", type=Path, default=RESULTS)
args = parser.parse_args()
DATA, RESULTS = args.data_dir, args.outdir


def read(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


candidate = read(DATA / "cell_type_candidate_counts.tsv")
brain = [row for row in candidate if row["tissue"] == "brain"]
validation = read(DATA / "fans_validation_denominators.tsv")
if {row["category"] for row in brain} != {"Neuron", "Oligodendrocyte"}:
    raise SystemExit("Unexpected brain categories")

fig, axes = plt.subplots(1, 2, figsize=(8.4, 5.3), gridspec_kw={"width_ratios": [1, 1.15]})
colors = {"Neuron": "#1B9E77", "Oligodendrocyte": "#0072B2"}
bars = axes[0].bar([row["category"] for row in brain], [int(row["count"]) for row in brain],
                   color=[colors[row["category"]] for row in brain])
axes[0].bar_label(bars, fontweight="bold", fontsize=12, padding=4)
axes[0].set_ylabel("Variant candidates")
axes[0].grid(axis="y", alpha=0.25)

tr = next(row for row in validation if row["tissue"] == "brain" and
          row["variant_class"] == "TR" and row["counting_unit"] == "unique_loci")
yes, no = int(tr["observed"]), int(tr["not_observed"])
bars = axes[1].bar(["Observed in\npaired FANS", "Not observed"], [yes, no], color=["#1B9E77", "#B8C2CC"])
axes[1].bar_label(bars, fontweight="bold", fontsize=12, padding=4)
axes[1].set_ylabel("Brain TR loci (unique)")
axes[1].grid(axis="y", alpha=0.25)
for ax in axes:
    ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
RESULTS.mkdir(parents=True, exist_ok=True)
fig.savefig(RESULTS / "brain_counts_and_fans_validation.png", dpi=300)
fig.savefig(RESULTS / "brain_counts_and_fans_validation.pdf", metadata={"CreationDate": None, "ModDate": None})
plt.close(fig)

with (RESULTS / "manuscript_count_reconciliation.tsv").open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["tissue", "counting_unit", "observed", "eligible", "note"])
    for row in validation:
        if row["variant_class"] == "TR":
            writer.writerow([row["tissue"], row["counting_unit"], row["observed"], row["eligible"],
                             "Manual FANS annotations; unique loci use any observed donor"])
