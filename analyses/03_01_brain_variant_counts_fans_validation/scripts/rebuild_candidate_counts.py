#!/usr/bin/env python3
"""Rebuild candidate counts from private SniffCell outputs; export no calls.

Adapted from build_selfcontained_sv_germline_filtered_figures.py and
build_figure2E_top_panel_fig13_matched.py. See README for input provenance.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import Counter
from pathlib import Path


GROUPS = {
    "brain": ({"Neuron"}, {"Oligodendrocyte"}),
    "pbmc": ({"T-cell", "NK-cell", "B-cell"}, {"Monocyte"}),
}


def read_table(path, required):
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing input columns: {sorted(missing)}")
        return list(reader)


def write_table(path, fields, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fingerprint(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def number(value):
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("Non-finite numeric input")
    return result


def boolean(value):
    if value not in {"True", "False"}:
        raise ValueError("Boolean columns must contain True or False")
    return value == "True"


def cell_set(value):
    return {cell for group in value.split("|") for part in group.split("+")
            for cell in part.strip().split(".")[-1].split("_") if cell}


def category(tissue, value):
    cells = cell_set(value)
    a, b = GROUPS[tissue]
    if not cells or cells - (a | b) or bool(cells & a) == bool(cells & b):
        raise ValueError("Cell-type assignment is missing, unknown, or spans both groups")
    return ("Neuron" if cells & a else "Oligodendrocyte") if tissue == "brain" else (
        "Lymphocyte group" if cells & a else "Monocyte")


def read_shared(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError("Malformed shared-SV VCF row")
            info = dict(item.split("=", 1) for item in fields[7].split(";") if "=" in item)
            rows.append((fields[2], fields[0], int(fields[1]), info["SVTYPE"],
                         abs(int(number(info.get("SVLEN", "0").split(",")[0])))))
    return rows


def shared_match(row, shared):
    pos, size = int(number(row["sv_pos"])), abs(int(number(row["sv_len"])))
    for sid, chrom, other_pos, kind, other_size in shared:
        if row["id"] not in {"", "."} and row["id"] == sid:
            return True
        largest = max(size, other_size)
        if (row["sv_chr"] == chrom and row["sv_type"] == kind
                and abs(pos - other_pos) <= max(1000, int(0.5 * largest))
                and (not largest or abs(size - other_size) / largest <= 0.5)):
            return True
    return False


def rebuild(manifest_path, outdir, sitecov=15):
    manifest_path = Path(manifest_path)
    config = json.loads(manifest_path.read_text())
    provenance = []

    def source(role, value, sample="", coverage=""):
        path = Path(value)
        if not path.is_absolute():
            path = manifest_path.parent / path
        if not path.is_file():
            raise ValueError(f"Missing required input for {role}, {sample}, {coverage}")
        provenance.append(dict(role=role, sample=sample, coverage=coverage,
                               bytes=path.stat().st_size, sha256=fingerprint(path)))
        return path

    # These are recorded adjudication decisions, not a replacement for the
    # unavailable original code that calculated the upstream artifact metrics.
    decision_rows = read_table(source("brain_v2_decisions", config["brain_v2_decisions"]),
                              ["sample", "cell_type", "trid", "recommended_v2_keep_with_rescue",
                               "recommended_v2_reason"])
    selected = {}
    decisions = Counter()
    seen = set()
    for row in decision_rows:
        key = (row["sample"], row["trid"], category("brain", row["cell_type"]))
        if key in seen:
            raise ValueError("Duplicate brain sample/locus/cell-type adjudication")
        seen.add(key)
        keep = boolean(row["recommended_v2_keep_with_rescue"])
        decisions[(row["sample"], "kept" if keep else "removed", row["recommended_v2_reason"])] += 1
        if keep:
            selected.setdefault(row["sample"], set()).add(key[1:])

    counts, qc = [], []
    run_keys = set()
    for run in config["runs"]:
        sample, tissue, coverage = run["sample"], run["tissue"], run["coverage"]
        if tissue not in GROUPS or coverage not in {"full", "30x", "60x"}:
            raise ValueError("Unknown tissue or coverage")
        if (sample, coverage) in run_keys:
            raise ValueError("Duplicate sample/coverage in manifest")
        run_keys.add((sample, coverage))
        by_class = Counter()
        audit = Counter()
        if tissue == "brain" and coverage == "full":
            if sample not in {r["sample"] for r in decision_rows}:
                raise ValueError("Brain sample missing from adjudication table")
            for _, cell in selected.get(sample, set()):
                by_class[("TR", cell)] += 1
            audit["TR_adjudicated_kept"] = sum(by_class.values())
        else:
            rows = read_table(source("tr_changes", run["tr_changes"], sample, coverage),
                              ["trid", "change_group", "n_change_support_reads",
                               "n_change_reads", "n_baseline_reads"])
            recovered = set()
            for row in rows:
                audit["TR_input_rows"] += 1
                cell = category(tissue, row["change_group"])
                if tissue == "brain":
                    key = (row["trid"], cell)
                    if key not in selected.get(sample, set()):
                        audit["TR_not_in_full_depth_selected_set"] += 1
                        continue
                    if key in recovered:
                        audit["TR_duplicate_recovery_row"] += 1
                        continue
                    recovered.add(key)
                else:
                    try:
                        keep = (number(row["n_change_support_reads"]) >= 5
                                and number(row["n_change_reads"]) >= 5
                                and number(row["n_baseline_reads"]) >= 10)
                    except ValueError:
                        audit["TR_invalid_numeric"] += 1
                        continue
                    if not keep:
                        audit["TR_below_support_threshold"] += 1
                        continue
                by_class[("TR", cell)] += 1
                audit["TR_kept"] += 1

        shared = read_shared(source("shared_sv_vcf", run["shared_sv_vcf"], sample, coverage))
        rows = read_table(source("sv_assignment", run["sv_assignment"], sample, coverage),
                          ["id", "variant_class", "has_hard_conflict", "linked_celltypes",
                           "primary_celltype", "n_overlapped", "sv_chr", "sv_pos", "sv_len", "sv_type"])
        for row in rows:
            audit["SV_assignment_input_rows"] += 1
            if row["variant_class"].upper() != "SV":
                audit["SV_non_sv_row"] += 1
                continue
            if boolean(row["has_hard_conflict"]):
                audit["SV_hard_conflict"] += 1
                continue
            a, b = GROUPS[tissue]
            linked = cell_set(row["linked_celltypes"])
            if linked & a and linked & b:
                audit["SV_both_cell_groups"] += 1
                continue
            try:
                cov = number(row["n_overlapped"])
                number(row["sv_pos"])
                number(row["sv_len"])
            except ValueError:
                audit["SV_invalid_numeric"] += 1
                continue
            if cov < sitecov:
                audit["SV_below_site_coverage"] += 1
                continue
            if shared_match(row, shared):
                audit["SV_shared_germline_match"] += 1
                continue
            cell = category(tissue, row["primary_celltype"] or row["linked_celltypes"])
            by_class[("SV", cell)] += 1
            audit["SV_kept"] += 1
        cells = ["Neuron", "Oligodendrocyte"] if tissue == "brain" else ["Lymphocyte group", "Monocyte"]
        for variant in ["TR", "SV"]:
            for cell in cells:
                counts.append(dict(tissue=tissue, sample=sample, coverage=coverage,
                                   variant_class=variant, cell_type=cell, count=by_class[(variant, cell)]))
        qc.extend(dict(sample=sample, coverage=coverage, stage=stage, count=count)
                  for stage, count in sorted(audit.items()))

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    write_table(outdir / "candidate_counts_by_sample.tsv",
                ["tissue", "sample", "coverage", "variant_class", "cell_type", "count"], counts)
    write_table(outdir / "filter_audit.tsv", ["sample", "coverage", "stage", "count"], qc)
    write_table(outdir / "brain_adjudication_audit.tsv", ["sample", "decision", "reason", "count"],
                [dict(sample=s, decision=d, reason=r, count=n) for (s, d, r), n in sorted(decisions.items())])
    write_table(outdir / "input_provenance.tsv", ["role", "sample", "coverage", "bytes", "sha256"], provenance)
    totals = Counter()
    for row in counts:
        if row["coverage"] == "full":
            totals[(row["tissue"], row["cell_type"])] += row["count"]
    aggregate = []
    for (tissue, cell), count in sorted(totals.items()):
        total = sum(n for (t, _), n in totals.items() if t == tissue)
        aggregate.append(dict(source="figure13_full_depth_sitecov15_v2_tr_matched" if sitecov == 15 else "configured_sitecov_v2_tr_matched",
                              tissue=tissue, display_tissue=tissue.capitalize() if tissue == "brain" else "PBMC",
                              category=cell, count=count, total=total,
                              percent=round(100 * count / total, 1) if total else 0))
    write_table(outdir / "cell_type_candidate_counts.tsv",
                ["source", "tissue", "display_tissue", "category", "count", "total", "percent"], aggregate)
    (outdir / "parameters.json").write_text(json.dumps(dict(
        sv_site_coverage=sitecov, pbmc_tr_support=5, pbmc_tr_change_reads=5, pbmc_tr_baseline_reads=10,
        brain_full="recorded v2 adjudication", brain_subsets="exact full-depth locus and cell-type recovery",
        sv_matching="ID or same chromosome/type, distance <= max(1000, 0.5*max_length), length difference <= 50%"), indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--sv-site-coverage", type=int, default=15)
    args = parser.parse_args()
    if args.sv_site_coverage < 0:
        parser.error("Site coverage must be nonnegative")
    rebuild(args.manifest, args.outdir, args.sv_site_coverage)
