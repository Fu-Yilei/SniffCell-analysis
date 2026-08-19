# Shared SniffCell read assignments for analyses 1.7-1.8

This folder describes the original per-read SniffCell assignment tables for
the six bulk samples and sixteen purified FACS/FANS samples used by the
read-level paper benchmarks.

The directory-preserving Zenodo archive contains:

- `bulk/`: four brain and two PBMC `deconv_read_summary.tsv` tables;
- `fans/`: eight brain and eight PBMC FACS/FANS tables;
- a path-free sample manifest and README.

These are original SniffCell outputs. They contain read identifiers and
methylation-based cell-type assignments, but no server paths, BAMs, or variant
calls.

## Download

After Zenodo record [10.5281/zenodo.22003085](https://doi.org/10.5281/zenodo.22003085)
is published, run:

```bash
python scripts/download_zenodo_data.py
```

The script downloads `sniffcell_read_assignments.tar.zst`, verifies its size
and SHA-256 checksum, validates the single archive root, and extracts
`sniffcell_read_assignments/` under `data/`.

The GitHub folder retains paper-section numbering. The Zenodo filename and
internal root intentionally do not use section numbers.
