#!/usr/bin/env python3
"""Calculate row-level and unique-locus manual FANS validation denominators."""
import argparse
from collections import Counter, defaultdict
from pathlib import Path

from rebuild_candidate_counts import read_table, write_table, fingerprint


def summarize(source, outdir):
    rows = read_table(source, ["variant_type", "tissue", "donor", "chrom", "start", "end",
                               "observed_in_paired_FANS"])
    seen = set()
    for row in rows:
        if row["variant_type"] not in {"TR", "SV"} or row["tissue"] not in {"brain", "pbmc"}:
            raise ValueError("Unexpected variant type or tissue")
        start, end = int(row["start"]), int(row["end"])
        if start < 0 or end < start:
            raise ValueError("Invalid coordinates")
        # The source defines unique loci by chromosome/start/end within class.
        key = (row["variant_type"], row["tissue"], row["donor"], row["chrom"], start, end)
        if key in seen:
            raise ValueError("Duplicate donor/locus record; resolve before counting")
        seen.add(key)
        status = row["observed_in_paired_FANS"].strip()
        row["status"] = status if status in {"Yes", "No"} else "Unscored"

    summary = []
    for tissue in ["all", "brain", "pbmc"]:
        for kind in ["TR", "SV"]:
            subset = [r for r in rows if r["variant_type"] == kind and (tissue == "all" or r["tissue"] == tissue)]
            loci = defaultdict(list)
            for r in subset:
                loci[(r["chrom"], int(r["start"]), int(r["end"]))].append(r["status"])
            for unit in ["donor_locus_rows", "unique_loci"]:
                statuses = [r["status"] for r in subset] if unit == "donor_locus_rows" else [
                    "Yes" if "Yes" in values else ("No" if "No" in values else "Unscored")
                    for values in loci.values()]
                c = Counter(statuses)
                eligible = c["Yes"] + c["No"]
                summary.append(dict(tissue=tissue, variant_class=kind, counting_unit=unit,
                                    observed=c["Yes"], not_observed=c["No"], unscored=c["Unscored"],
                                    eligible=eligible, total=len(statuses),
                                    percent_observed=round(100*c["Yes"]/eligible, 3) if eligible else "NA",
                                    loci_with_mixed_yes_no=sum("Yes" in v and "No" in v for v in loci.values())))
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    write_table(outdir / "fans_validation_denominators.tsv", list(summary[0]), summary)
    legacy = Counter((r["variant_type"], r["status"]) for r in rows)
    write_table(outdir / "fans_validation_counts.tsv", ["variant_class", "validation_status", "count"],
                [dict(variant_class=k, validation_status=s, count=n) for (k, s), n in sorted(legacy.items())])
    write_table(outdir / "validation_input_provenance.tsv", ["role", "sha256", "rows"],
                [dict(role="manual_fans_annotations", sha256=fingerprint(source), rows=len(rows))])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotations", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    summarize(args.annotations, args.outdir)
