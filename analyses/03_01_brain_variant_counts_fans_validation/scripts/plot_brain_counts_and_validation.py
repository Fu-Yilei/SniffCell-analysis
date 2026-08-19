#!/usr/bin/env python3
"""Plot aggregate brain candidate counts and paired-FANS validation."""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
DATA = HERE / "data"
RESULTS = HERE / "results"


def read(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


candidate = read(DATA / "cell_type_candidate_counts.tsv")
brain = [row for row in candidate if row["tissue"] == "brain"]
validation = read(DATA / "fans_validation_counts.tsv")
if {row["category"] for row in brain} != {"Neuron", "Oligodendrocyte"}:
    raise SystemExit("Unexpected brain categories")

fig, axes = plt.subplots(1, 2, figsize=(8.4, 5.3), gridspec_kw={"width_ratios": [1, 1.15]})
colors = {"Neuron": "#1B9E77", "Oligodendrocyte": "#0072B2"}
bars = axes[0].bar([row["category"] for row in brain], [int(row["count"]) for row in brain],
                   color=[colors[row["category"]] for row in brain])
axes[0].bar_label(bars, fontweight="bold", fontsize=12, padding=4)
axes[0].set_ylabel("Variant candidates")
axes[0].set_title("Figure 13 full-depth candidates\nBrain (n=109)")
axes[0].grid(axis="y", alpha=0.25)

tr = {row["validation_status"]: int(row["count"]) for row in validation if row["variant_class"] == "TR"}
yes, no = tr.get("Yes", 0), tr.get("No", 0)
bars = axes[1].bar(["Observed in\npaired FANS", "Not observed"], [yes, no], color=["#1B9E77", "#B8C2CC"])
axes[1].bar_label(bars, fontweight="bold", fontsize=12, padding=4)
axes[1].set_ylabel("TR candidates")
axes[1].set_title(f"Paired FANS validation\n{yes}/{yes + no} observed")
axes[1].grid(axis="y", alpha=0.25)
for ax in axes:
    ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
RESULTS.mkdir(parents=True, exist_ok=True)
fig.savefig(RESULTS / "brain_counts_and_fans_validation.png", dpi=300)
fig.savefig(RESULTS / "brain_counts_and_fans_validation.pdf")
plt.close(fig)

reported = {row["metric"]: int(row["reported_count"]) for row in read(DATA / "manuscript_reported_counts.tsv")}
current = {
    "figure13_matched_neuron_all_classes": int(next(row["count"] for row in brain if row["category"] == "Neuron")),
    "figure13_matched_oligodendrocyte_all_classes": int(next(row["count"] for row in brain if row["category"] == "Oligodendrocyte")),
    "fans_tr_validated": yes,
    "fans_tr_total": yes + no,
}
with (RESULTS / "manuscript_count_reconciliation.tsv").open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["metric", "manuscript_reported_count", "current_source_count", "status", "note"])
    writer.writerow(["fans_tr_validated", reported["fans_tr_validated"], current["fans_tr_validated"], "match" if reported["fans_tr_validated"] == current["fans_tr_validated"] else "discrepancy", "current validation table"])
    writer.writerow(["fans_tr_total", reported["fans_tr_total"], current["fans_tr_total"], "match" if reported["fans_tr_total"] == current["fans_tr_total"] else "discrepancy", "current validation table"])
    writer.writerow(["brain_neuron_all_classes", "NA", current["figure13_matched_neuron_all_classes"], "different definition", "Figure 13-matched source combines TR and SV"])
    writer.writerow(["brain_oligodendrocyte_all_classes", "NA", current["figure13_matched_oligodendrocyte_all_classes"], "different definition", "Figure 13-matched source combines TR and SV"])
