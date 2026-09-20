#!/usr/bin/env python3
"""Recompute historical brain TR selection from per-locus sequence evidence.

This is the historical Figure 2e protocol, not a reinterpretation of the current
Methods. Detailed output is private; public output contains only counts/hashes.
"""
import argparse
from collections import Counter
import json
import math
from pathlib import Path

from rebuild_candidate_counts import read_table, write_table, fingerprint, boolean


def numeric(value, *, optional=False):
    if optional and value in {"", "nan", "NA"}:
        return None
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("Non-finite evidence value")
    return number


def mechanistic_filter(row):
    support, delta = numeric(row["support"]), numeric(row["delta"])
    if support < 0 or support != int(support) or delta < 0:
        raise ValueError("Invalid support count or expansion length")
    repro = boolean(row["reproducible"])
    gc, entropy, top4, dinuc = [numeric(row[k], optional=True) for k in (
        "support_gc_median", "support_entropy1_median",
        "support_top4mer_frac_median", "support_dinuc_purity_median")]
    score = sum((gc is not None and gc < .25, entropy is not None and entropy < .86,
                 top4 is not None and top4 > .12, dinuc is not None and dinuc > .28))
    if support >= 5:
        keep, reason = True, "support_ge5"
    elif support in (3, 4):
        limit = 2000 if support == 3 else 4000
        keep = (delta <= limit or repro) and not (score >= 3 and not repro)
        reason = f"support{int(support)}_" + (
            "scaled_delta_or_repro_no_unrep_low_complex" if keep else "artifact_signature")
    elif support == 2:
        keep = delta <= 1000 and repro and score < 3 and (gc is None or gc >= .25)
        reason = "support2_modest_repro_complex" if keep else "support2_artifact_signature"
    else:
        keep, reason = False, "below_support"
    return keep, reason, score


def reproducibility_class(row, runs):
    sample, cell, locus = row["sample"], row["cell_type"], row["trid"]
    exact = [r for r in runs if r["trid"] == locus and r["cell_type"] == cell]
    depth30 = any(r["sample"] == sample and r["coverage"] == "30x" for r in exact)
    depth60 = any(r["sample"] == sample and r["coverage"] == "60x" for r in exact)
    cross = [r for r in runs if r["sample"] != sample and r["coverage"] == "full"]
    near = [r for r in cross if r["chrom"] == row["chrom"]
            and int(r["start"]) <= int(row["end"]) + 1000
            and int(r["end"]) >= int(row["start"]) - 1000]
    if any(r["sample"] != sample and r["coverage"] == "full" for r in exact):
        return "cross_donor_same_celltype"
    if any(r["cell_type"] == cell for r in near):
        return "near_cross_donor_same_celltype"
    if depth30 and depth60:
        return "same_sample_60x_30x_only"
    if depth60:
        return "same_sample_60x_only"
    if near or any(r["trid"] == locus and r["cell_type"] != cell for r in cross):
        return "cross_donor_diff_or_unstable"
    return "private_full_only"


def derive(manifest, outdir):
    manifest = Path(manifest)
    config = json.loads(manifest.read_text())
    provenance = []

    def source(role, value, sample="", coverage=""):
        path = Path(value)
        if not path.is_absolute():
            path = manifest.parent / path
        provenance.append(dict(role=role, sample=sample, coverage=coverage,
                               bytes=path.stat().st_size, sha256=fingerprint(path)))
        return path

    metrics = read_table(source("brain_sequence_metrics", config["brain_sequence_metrics"]),
        ["sample", "cell_type", "trid", "chrom", "start", "end", "support", "delta",
         "reproducible", "support_gc_median", "support_entropy1_median",
         "support_top4mer_frac_median", "support_dinuc_purity_median"])
    runs, run_keys = [], set()
    for run in config["brain_reproducibility_runs"]:
        key = run["sample"], run["coverage"]
        if key in run_keys or key[1] not in {"full", "30x", "60x"}:
            raise ValueError("Duplicate or invalid brain run")
        run_keys.add(key)
        for row in read_table(source("brain_reproducibility", run["tr_changes"], *key),
                              ["chrom", "start", "end", "trid", "change_group"]):
            cell = row["change_group"].rsplit(".", 1)[-1]
            if cell not in {"Neuron", "Oligodendrocyte"}:
                raise ValueError("Invalid brain cell type")
            runs.append(dict(row, sample=key[0], coverage=key[1], cell_type=cell))
    samples = {r["sample"] for r in metrics}
    if run_keys != {(s, c) for s in samples for c in ("full", "30x", "60x")}:
        raise ValueError("Reproducibility requires all three depths for every brain donor")
    rescues = config.get("curated_rescues", [])
    rescue_keys = {(r["sample"], r["cell_type"], r["trid"]) for r in rescues}
    if len(rescue_keys) != len(rescues) or any(not r.get("reason") for r in rescues):
        raise ValueError("Rescues need unique keys and explicit reasons")
    evidence, decisions, seen = [], [], set()
    summary = Counter()
    for row in metrics:
        key = row["sample"], row["cell_type"], row["trid"]
        if key in seen or row["cell_type"] not in {"Neuron", "Oligodendrocyte"}:
            raise ValueError("Duplicate or invalid sample/locus/cell-type evidence")
        seen.add(key)
        mechanical, mechanical_reason, score = mechanistic_filter(row)
        repro = reproducibility_class(row, runs)
        main = mechanical and score < 3 and not (
            repro == "same_sample_60x_only" and numeric(row["support"]) < 5)
        rescue = key in rescue_keys
        if rescue and not mechanical:
            raise ValueError("Historical rescue is defined only after mechanistic selection")
        keep = main or rescue
        reason = ("main_keep" if main else "curated_rescue" if rescue else
                  "drop_mechanistic" if not mechanical else "drop_residual_low_complex"
                  if score >= 3 else "drop_same_sample_60x_only_support_lt5")
        record = dict(row, low_complex_score=score, mechanistic_keep=mechanical,
                      mechanistic_reason=mechanical_reason, repro_class=repro,
                      recommended_v2_main_keep=main, curated_rescue=rescue,
                      recommended_v2_keep_with_rescue=keep, recommended_v2_reason=reason)
        evidence.append(record)
        if mechanical:
            decisions.append(record)
        summary[(row["sample"], reason)] += 1
    if rescue_keys - seen:
        raise ValueError("Rescue key absent from evidence")
    outdir = Path(outdir)
    private = outdir / "private_inputs"
    private.mkdir(parents=True, exist_ok=True)
    # Never publish these detailed tables: they contain candidate coordinates.
    fields = list(evidence[0])
    write_table(private / "brain_filter_all_loci.tsv", fields, evidence)
    write_table(private / "brain_filter_decisions.tsv", fields, decisions)
    write_table(outdir / "brain_filter_stage_counts.tsv", ["sample", "decision", "count"],
                [dict(sample=s, decision=d, count=n) for (s, d), n in sorted(summary.items())])
    write_table(outdir / "brain_filter_input_provenance.tsv",
                ["role", "sample", "coverage", "bytes", "sha256"], provenance)
    return private / "brain_filter_decisions.tsv"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    derive(args.manifest, args.outdir)
