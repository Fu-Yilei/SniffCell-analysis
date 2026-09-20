#!/usr/bin/env python3
"""Run sections 3.1 and 3.2 from private inputs to aggregate tables and plots."""
import argparse
import csv
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
MOTIF = HERE.parents[1] / "03_02_brain_tr_motif_composition" / "scripts"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--fans-annotations", type=Path, required=True)
    parser.add_argument("--motif-candidates", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    counts, motif = args.outdir / "counts_and_validation", args.outdir / "motif"
    commands = [
        [HERE / "rebuild_candidate_counts.py", "--manifest", args.manifest, "--outdir", counts],
        [HERE / "summarize_fans_validation.py", "--annotations", args.fans_annotations, "--outdir", counts],
        [MOTIF / "derive_tr_motif_composition.py", "--candidates", args.motif_candidates, "--outdir", motif],
    ]
    for command in commands:
        subprocess.run([sys.executable, *map(str, command)], check=True)
    with (counts / "candidate_counts_by_sample.tsv").open() as handle:
        candidate_rows = list(csv.DictReader(handle, delimiter="\t"))
    with (motif / "tr_motif_composition.tsv").open() as handle:
        motif_rows = list(csv.DictReader(handle, delimiter="\t"))
    for tissue, label in [("brain", "Brain"), ("pbmc", "PBMC")]:
        selected = sum(int(row["count"]) for row in candidate_rows if row["tissue"] == tissue
                       and row["coverage"] == "full" and row["variant_class"] == "TR")
        motif_total = next(int(row["tissue_total"]) for row in motif_rows if row["tissue"] == label)
        if selected != motif_total:
            raise ValueError(f"{label}: selected TR rows ({selected}) differ from motif rows ({motif_total}); resolve source versions")
    for command in [
        [HERE / "plot_brain_counts_and_validation.py", "--data-dir", counts, "--outdir", counts],
        [MOTIF / "plot_tr_motif_composition.py", "--input", motif / "tr_motif_composition.tsv", "--outdir", motif],
    ]:
        subprocess.run([sys.executable, *map(str, command)], check=True)
    print("Completed candidate filtering, manual-validation denominators, motif analysis and plots.")


if __name__ == "__main__":
    main()
