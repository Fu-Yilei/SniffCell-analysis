#!/usr/bin/env python3
"""Measure the original curated-callset FANS CIGAR proxy from indexed BAMs.

This automated proxy does not recreate or overwrite manual validation labels.
Only aggregate statuses and read counts leave this script, never call records.
"""
import argparse
from collections import Counter
from contextlib import ExitStack
import json
import re
import hashlib
from pathlib import Path
from statistics import median

import pysam

from rebuild_candidate_counts import cell_set, read_table, write_table, fingerprint


def spanning_reads(bam, chrom, lo, hi):
    if chrom not in bam.references:
        raise ValueError("Candidate contig missing from BAM; verify genome build")
    for read in bam.fetch(chrom, max(0, lo - 150), hi + 150):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.reference_start <= lo and read.reference_end >= hi:
            yield read


def indels(read):
    pos = read.reference_start
    for op, length in read.cigartuples or []:
        if op in (0, 7, 8, 3):
            pos += length
        elif op == 2:
            yield pos, -length
            pos += length
        elif op == 1:
            yield pos, length


def tr_values(bam, chrom, start, end, pad):
    lo, hi = max(0, start - pad), end + pad
    return [sum(size for pos, size in indels(read) if lo <= pos <= hi)
            for read in spanning_reads(bam, chrom, lo, hi)]


def sv_values(bam, chrom, start, size, window):
    reads = list(spanning_reads(bam, chrom, max(0, start - 20), start + 20))
    min_size = max(20, abs(size) * 0.5) if size else 30
    alt = sum(any(abs(pos-start) < window and abs(length) >= min_size
                  for pos, length in indels(read)) for read in reads)
    return len(reads), alt


def measure(candidates, manifest, outdir, pad=150, min_reads=5, min_delta=5):
    rows = read_table(candidates, ["sample", "variant_type", "cell_type", "chrom", "start", "end", "direction", "size"])
    manifest = Path(manifest)
    config = json.loads(manifest.read_text())
    status_counts, read_counts = Counter(), Counter()
    provenance = []
    with ExitStack() as stack:
        handles = {}
        for sample, record in config["samples"].items():
            if record["tissue"] not in {"brain", "pbmc"}:
                raise ValueError("Unknown tissue")
            if len(set(record["bams"].values())) != len(record["bams"]):
                raise ValueError("The same BAM cannot represent multiple sorted populations")
            handles[sample] = {}
            for cell, filename in record["bams"].items():
                path = Path(filename)
                if not path.is_absolute():
                    path = manifest.parent / path
                bam = stack.enter_context(pysam.AlignmentFile(path, "rb"))
                if not bam.has_index():
                    raise ValueError("FANS BAM index missing")
                handles[sample][cell] = bam
                provenance.append(dict(sample=sample, cell_type=cell, bytes=path.stat().st_size,
                                       header_sha256=hashlib.sha256(str(bam.header).encode()).hexdigest()))
        for row in rows:
            sample = row["sample"]
            # Missing donor matches are not negative validation evidence.
            tissue = config["samples"].get(sample, {}).get("tissue", "unmapped")
            kind = "TR" if "TR" in row["variant_type"] and any(
                x in row["direction"].lower() for x in ["expansion", "contraction"]) else "SV"
            status = "unscored"
            if sample not in handles:
                status = "no_matched_FANS"
            elif "recomb" in row["direction"].lower():
                status = "requires_structural_review"
            else:
                bams = handles[sample]
                aliases = config.get("cell_type_aliases", {})
                changed = {aliases.get(c, c) for c in cell_set(row["cell_type"])}
                contrast = set(bams) - changed
                if not changed or changed - set(bams) or not contrast:
                    status_counts[(tissue, sample, kind, "unmapped_cell_type")] += 1
                    continue
                start, end = int(row["start"]), int(row["end"])
                if start < 0 or end < start:
                    raise ValueError("Invalid candidate coordinates")
                if kind == "TR":
                    direction = row["direction"].lower()
                    if not any(x in direction for x in ["expansion", "contraction"]):
                        raise ValueError("TR direction must specify expansion or contraction")
                    a = [x for cell in sorted(changed) for x in tr_values(bams[cell], row["chrom"], start, end, pad)]
                    b = [x for cell in sorted(contrast) for x in tr_values(bams[cell], row["chrom"], start, end, pad)]
                    nc, nb = len(a), len(b)
                    if min(nc, nb) < min_reads:
                        status = "low_coverage"
                    else:
                        delta = median(a) - median(b)
                        ok = delta >= min_delta if "expansion" in direction else delta <= -min_delta
                        status = "proxy_concordant" if ok else "proxy_not_concordant"
                else:
                    # Legacy curated table sometimes stores a size range.
                    size_match = re.match(r"[+-]?\d+", row["size"].strip())
                    size = int(size_match[0]) if size_match else end-start
                    a = [sv_values(bams[c], row["chrom"], start, size, pad) for c in sorted(changed)]
                    b = [sv_values(bams[c], row["chrom"], start, size, pad) for c in sorted(contrast)]
                    nc, nb = sum(x[0] for x in a), sum(x[0] for x in b)
                    if min(nc, nb) < min_reads:
                        status = "low_coverage"
                    else:
                        ok = sum(x[1] for x in a)/nc >= 0.15 and sum(x[1] for x in b)/nb <= 0.05
                        status = "proxy_concordant" if ok else "proxy_not_concordant"
                read_counts[(tissue, sample, kind, "changed")] += nc
                read_counts[(tissue, sample, kind, "contrast")] += nb
            status_counts[(tissue, sample, kind, status)] += 1
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    write_table(outdir / "bam_input_provenance.tsv", ["sample", "cell_type", "bytes", "header_sha256"], provenance)
    write_table(outdir / "candidate_input_provenance.tsv", ["role", "sha256", "rows"],
                [dict(role="curated_bam_validation_candidates", sha256=fingerprint(candidates), rows=len(rows))])
    write_table(outdir / "automated_fans_proxy_counts.tsv", ["tissue", "sample", "variant_class", "status", "count"],
                [dict(tissue=t, sample=s, variant_class=k, status=v, count=n)
                 for (t, s, k, v), n in sorted(status_counts.items())])
    write_table(outdir / "automated_fans_read_counts.tsv", ["tissue", "sample", "variant_class", "group", "read_locus_observations"],
                [dict(tissue=t, sample=s, variant_class=k, group=g, read_locus_observations=n)
                 for (t, s, k, g), n in sorted(read_counts.items())])
    (outdir / "automated_fans_parameters.json").write_text(json.dumps(dict(
        pad_bp=pad, minimum_reads_per_group=min_reads, minimum_tr_delta_bp=min_delta,
        sv_changed_vaf_min=0.15, sv_contrast_vaf_max=0.05,
        interpretation="CIGAR proxy only; not manual validation or sequence-resolved genotyping"), indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--pad", type=int, default=150)
    parser.add_argument("--min-reads", type=int, default=5)
    parser.add_argument("--min-delta", type=int, default=5)
    args = parser.parse_args()
    if min(args.pad, args.min_reads, args.min_delta) < 1:
        parser.error("Thresholds must be positive")
    measure(args.candidates, args.manifest, args.outdir, args.pad, args.min_reads, args.min_delta)
