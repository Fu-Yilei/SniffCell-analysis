#!/usr/bin/env python3
"""Summarize FGF14 5'-CFV allele lengths and consensus interruptions."""

import csv
from hashlib import sha256
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
DATA = HERE / "data"
RESULTS = HERE / "results"
ASSET = RESULTS / "extended_data_figure_1_source.png"
EXPECTED = "04f2786958a9d2c5762b72db171bf62632bca19de5a3f462abc88a3e57541e8b"
DONORS = {"b1": "bcontrol1", "C2": "bcontrol2", "b2": "bmsa", "P2": "bparkinsons"}
DONOR_ORDER = {"b1": 0, "C2": 1, "b2": 2, "P2": 3}
CELL_ORDER = {"Neuron": 0, "Oligodendrocyte": 1}
COLORS = {"5-CFV": "#2B6CB0", "Non-CFV": "#D95F4F"}


def read(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def interval(value):
    if not value or value in {"n/a", "(combined)"}:
        return None
    left, right = value.split("-")
    return float(left), float(right)


digest = sha256(ASSET.read_bytes()).hexdigest()
if digest != EXPECTED:
    raise SystemExit(f"Extended Data Figure 1 checksum mismatch: {digest}")

rows = read(DATA / "fgf14_cfv_benchmark_brains.tsv")
expected_donors = set(DONORS)
if {row["Donor"] for row in rows} != expected_donors:
    raise SystemExit("CFV table does not contain exactly the four benchmark donors")

fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.7), sharey=True)
for ax, source in zip(axes, ("SniffCell", "FANS")):
    subset = sorted(
        (row for row in rows if row["Source"] == source),
        key=lambda row: (DONOR_ORDER[row["Donor"]], CELL_ORDER[row["Cell type"]]),
    )
    labels = []
    for i, row in enumerate(subset):
        labels.append(f"{DONORS[row['Donor']]}\n{row['Cell type'].replace('Oligodendrocyte', 'Oligo')}")
        for offset, (label, field) in enumerate((("5-CFV", "5-CFV allele read len (bp)"), ("Non-CFV", "Non-CFV allele read len (bp)"))):
            parsed = interval(row[field])
            if parsed is None:
                continue
            low, high = parsed
            x = i + (-0.13 if label == "5-CFV" else 0.13)
            ax.errorbar(x, (low + high) / 2, yerr=(high - low) / 2, fmt="o", color=COLORS[label], capsize=3, label=label if i == 0 else None)
    ax.set_xticks(range(len(subset)), labels, rotation=45, ha="right", fontsize=8)
    ax.set_title(source if source == "SniffCell" else "FANS")
    ax.grid(axis="y", alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)
axes[0].set_ylabel("Flank-corrected FGF14 allele length (bp)")
axes[0].legend(frameon=False)
fig.tight_layout()
RESULTS.mkdir(parents=True, exist_ok=True)
fig.savefig(RESULTS / "fgf14_cfv_range_summary.png", dpi=300)
fig.savefig(RESULTS / "fgf14_cfv_range_summary.pdf")
plt.close(fig)

consensus = read(DATA / "fgf14_consensus_interruptions.tsv")
fans = read(DATA / "fgf14_fans_motif_summary.tsv")
with (RESULTS / "fgf14_interruption_summary.tsv").open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["evidence_source", "donor", "cell_type", "expanded_allele", "interruption_motifs", "supporting_reads", "total_reads"])
    for row in consensus:
        if int(row["allele_index"]) != 2:
            continue
        writer.writerow(["SniffCell TRGT consensus", row["display_sample"], row["cell_type"], row["allele_index"], row["interruption_motifs"], "NA", row["support_depth"]])
    for row in fans:
        interruption = row["expanded_consensus_interruptions"]
        if "AGG" in interruption:
            support = row["reads_with_double_AGG_signature"]
        else:
            support = row["reads_with_at_least_one_AAAGAAGAAG"]
        writer.writerow(["FANS consensus and reads", DONORS.get(row["donor"], row["donor"]), row["cell_type"], row["expanded_hap"], interruption, support, row["expanded_reads_mapped_to_poa"]])
print(f"OK  {digest}  {ASSET.name}")
