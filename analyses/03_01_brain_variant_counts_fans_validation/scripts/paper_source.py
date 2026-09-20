"""Strict per-candidate table reader shared by the paper analyses."""
import csv
import hashlib
import json
from pathlib import Path

REQUIRED = {"variant_id", "variant_type", "tissue", "donor", "chrom", "start", "end",
            "cell_type", "observed_in_paired_FANS"}
CELLS = {"brain": {"Neuron", "Oligodendrocyte"},
         "pbmc": {"T-cell+NK-cell+B-cell", "Monocyte"}}


def records_sha256(rows):
    fields = sorted(REQUIRED | {"coverages_detected", "motif"})
    records = [{k: r.get(k, "") for k in fields} for r in sorted(rows, key=lambda r: r["variant_id"])]
    return hashlib.sha256(json.dumps(records, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def read_candidates(path):
    path = Path(path)
    if path.suffix == ".json":
        obj = json.loads(path.read_text())
        values = obj.get("structuredContent", obj)["values"]
    else:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            values = list(csv.reader(handle, delimiter="\t" if path.suffix == ".tsv" else ","))
    header_at = next((i for i, row in enumerate(values[:3]) if REQUIRED <= set(row)), None)
    if header_at is None:
        raise ValueError("Missing required per-candidate headers; aggregate tables are not accepted")
    headers = values[header_at]
    if len(headers) != len(set(headers)):
        raise ValueError("Duplicate column headers")
    rows, ids, keys = [], set(), set()
    for number, values_row in enumerate(values[header_at + 1:], header_at + 2):
        if not any(str(x).strip() for x in values_row):
            continue
        row = {k: str(values_row[i]).strip() if i < len(values_row) else ""
               for i, k in enumerate(headers)}
        if not row["variant_id"] or row["variant_id"] in ids:
            raise ValueError(f"Missing/duplicate candidate ID at source row {number}")
        if row["variant_type"] not in {"TR", "SV"} or row["tissue"] not in CELLS:
            raise ValueError(f"Invalid variant class/tissue at source row {number}")
        if row["cell_type"] not in CELLS[row["tissue"]]:
            raise ValueError(f"Invalid cell type at source row {number}")
        if row["observed_in_paired_FANS"] not in {"Yes", "No"}:
            raise ValueError(f"Unscored FANS annotation at source row {number}; do not treat it as No")
        if "coverages_detected" not in row:
            raise ValueError("Missing coverages_detected column")
        depths = row["coverages_detected"].split(";")
        if not depths or len(depths) != len(set(depths)) or set(depths) - {"full", "30x", "60x"}:
            raise ValueError(f"Invalid coverage membership at source row {number}")
        row["coverages_detected"] = ";".join(c for c in ("30x", "60x", "full") if c in depths)
        for k in ("start", "end"):
            value = float(row[k])
            if value != int(value) or value < 0:
                raise ValueError(f"Invalid coordinate at source row {number}")
            row[k] = str(int(value))
        if int(row["end"]) <= int(row["start"]):
            raise ValueError(f"Invalid interval at source row {number}")
        key = tuple(row[k] for k in ("variant_type", "tissue", "donor", "chrom", "start", "end", "cell_type"))
        if key in keys:
            raise ValueError(f"Duplicate donor-locus-cell record at source row {number}")
        keys.add(key)
        ids.add(row["variant_id"])
        row["source_row"] = number
        rows.append(row)
    if not rows:
        raise ValueError("No candidate records")
    return rows


def depth_counts(rows):
    samples = [("pbmc1", "pbmc"), ("pbmc2", "pbmc"), ("bcontrol1", "brain"),
               ("bcontrol2", "brain"), ("bmsa", "brain"), ("bparkinsons", "brain")]
    if {(r["donor"], r["tissue"]) for r in rows} != set(samples):
        raise ValueError("Expected exactly the six Table 5 donors and their tissues")
    out = []
    for sample, tissue in samples:
        for coverage in ("30x", "60x", "full"):
            group = [r for r in rows if r["donor"] == sample
                     and coverage in r["coverages_detected"].split(";")]
            out.append(dict(tissue=tissue, sample=sample, coverage=coverage,
                            SV=sum(r["variant_type"] == "SV" for r in group),
                            TR=sum(r["variant_type"] == "TR" for r in group), total=len(group)))
    return out


def summarize(rows, pbmc_composition=None):
    validation = []
    for kind in ("TR", "SV"):
        group = [r for r in rows if r["variant_type"] == kind]
        if not group:
            raise ValueError(f"Missing {kind} candidates")
        for status in ("Yes", "No"):
            n = sum(r["observed_in_paired_FANS"] == status for r in group)
            validation.append(dict(variant_class=kind, status=status, count=n,
                                   total=len(group), percent=100 * n / len(group)))
    composition = []
    for tissue, cells in (("pbmc", ("T-cell+NK-cell+B-cell", "Monocyte")),
                          ("brain", ("Neuron", "Oligodendrocyte"))):
        source = pbmc_composition if tissue == "pbmc" and pbmc_composition is not None else rows
        group = [r for r in source if r["tissue"] == tissue]
        if not group:
            raise ValueError(f"Missing {tissue} candidates")
        for cell in cells:
            n = sum(r["cell_type"] == cell for r in group)
            composition.append(dict(tissue=tissue, cell_type=cell, count=n,
                                    total=len(group), percent=100 * n / len(group)))
    # Unlike Figure 2f, the Results' brain validation sentence counts loci once
    # across donors and cell types. Never sum per-cell unique-locus denominators.
    brain = [r for r in rows if r["tissue"] == "brain" and r["variant_type"] == "TR"]
    loci = {}
    for r in brain:
        key = r["chrom"], r["start"], r["end"]
        loci[key] = loci.get(key, False) or r["observed_in_paired_FANS"] == "Yes"
    unique = dict(scope="brain_TR", unit="unique_genomic_locus_any_donor",
                  observed=sum(loci.values()), total=len(loci))
    return validation, composition, unique
