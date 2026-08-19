# Shared modkit CpG pileups for analyses 1.4-1.6

This folder defines the common methylation input package used by paper
analyses 1.4, 1.5, and 1.6. The package contains 28 path-free, gzip-compressed
CpG pileup tables: 12 SniffCell-assigned bulk group tables and 16 matched FANS
tables from six benchmark samples.

The original workflow ran modkit v0.6.1 with CpG-aware, strand-combined
pileups. It streamed the bedMethyl output into the four columns used by these
analyses, so full bedMethyl files were not retained. Each public table contains
no header and exactly:

1. chromosome
2. zero-based CpG position
3. valid coverage
4. modified count

These are methylation summaries, not variant calls. They contain neither read
identifiers nor server paths.

## Download

After Zenodo record 10.5281/zenodo.22003085 is published, run:

```bash
python scripts/download_zenodo_data.py
```

The command downloads one directory-preserving tar archive, verifies the
archive checksum, extracts its named root folder under `data/`, and verifies
all 28 files against `file_manifest.tsv`.

The GitHub directory keeps the paper-section numbering for navigation. The
Zenodo object intentionally does not: its filename is
`sniffcell_modkit_cpg_pileups.tar`, and its internal root is
`sniffcell_modkit_cpg_pileups/`.
