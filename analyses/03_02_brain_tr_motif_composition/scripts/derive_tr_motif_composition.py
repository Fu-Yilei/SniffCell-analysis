#!/usr/bin/env python3
"""Derive motif statistics from private TR rows, not precomputed counts.

Adapted from build_figure13_tr_motif_composition.py. A row is a donor-locus
candidate. Motif length is the annotated string length, not minimal period.
"""
import argparse
import csv
import hashlib
from collections import Counter
from pathlib import Path


LENGTH_BINS = [("2", 2, 2), ("3", 3, 3), ("4", 4, 4), ("5", 5, 5),
               ("6–10", 6, 10), (">10", 11, float("inf"))]


def sequence_class(motif):
    motif = motif.upper()
    if not motif or set(motif) - set("ACGT"):
        raise ValueError("Missing or non-ACGT motif; annotate explicitly before analysis")
    gc = sum(base in "GC" for base in motif) / len(motif)
    return "AT-rich" if gc <= 0.20 else ("GC-rich" if gc >= 0.80 else "Mixed")


def derive(source, outdir):
    with Path(source).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"variant_type", "tissue", "donor", "chrom", "start", "end", "motif"}
        if required - set(reader.fieldnames or []):
            raise ValueError(f"Missing input columns: {sorted(required-set(reader.fieldnames or []))}")
        variants = list(reader)
    counts, totals, seen = Counter(), Counter(), set()
    for row in variants:
        if row["variant_type"] not in {"TR", "SV"}:
            raise ValueError("Unexpected variant type")
        if row["variant_type"] != "TR":
            continue
        if row["tissue"] not in {"pbmc", "brain"}:
            raise ValueError("Unexpected tissue")
        tissue = {"pbmc": "PBMC", "brain": "Brain"}[row["tissue"]]
        start, end = int(row["start"]), int(row["end"])
        if start < 0 or end < start:
            raise ValueError("Invalid candidate coordinates")
        key = (tissue, row["donor"], row["chrom"], start, end)
        if key in seen:
            raise ValueError("Duplicate donor/locus TR row; resolve before counting")
        seen.add(key)
        motif = row["motif"].strip().upper()
        composition = sequence_class(motif)
        length_bin = next((name for name, low, high in LENGTH_BINS if low <= len(motif) <= high), None)
        if length_bin is None:
            raise ValueError("One-base motifs are outside the paper's bins; do not silently omit them")
        totals[tissue] += 1
        counts[(tissue, "Motif length (bp)", length_bin)] += 1
        counts[(tissue, "Sequence composition", composition)] += 1
    output = []
    for tissue in ["PBMC", "Brain"]:
        for dimension, categories in [("Motif length (bp)", [x[0] for x in LENGTH_BINS]),
                                       ("Sequence composition", ["AT-rich", "Mixed", "GC-rich"])]:
            for category in categories:
                n = counts[(tissue, dimension, category)]
                output.append(dict(tissue=tissue, dimension=dimension, category=category,
                                   count=n, tissue_total=totals[tissue],
                                   percent=round(100*n/totals[tissue], 1) if totals[tissue] else 0))
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / "tr_motif_composition.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    digest = hashlib.sha256(Path(source).read_bytes()).hexdigest()
    with (outdir / "motif_input_provenance.tsv").open("w") as handle:
        handle.write("role\tsha256\tinput_rows\ttr_rows\n")
        handle.write(f"full_depth_annotated_candidates\t{digest}\t{len(variants)}\t{sum(totals.values())}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    derive(args.candidates, args.outdir)
