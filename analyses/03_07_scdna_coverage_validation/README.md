# 3.7 scDNA coverage validation

This section supports Extended Data Figure 2a-b. It compares single-cell DNA
coverage outside and inside SniffCell-called brain variant intervals across
four samples and two sorted cell types.

The public analysis contains no BAM, VCF, variant-level call table, genomic
interval list, read identifier, or server path. The exact source panel is
retained as `results/extended_data_figure_2ab_source.png`.

## Public data

Download `scdna_coverage_summaries.tar.zst` from
[Zenodo record 22003085](https://doi.org/10.5281/zenodo.22003085), then extract
it:

```bash
tar --use-compress-program=unzstd -xf scdna_coverage_summaries.tar.zst
```

The archive contains sanitized per-run and per-group coverage summaries plus
the 528 selected `.mosdepth.summary.txt` files. Coordinate-bearing mosdepth
region and threshold BED files are intentionally excluded.

## Run

```bash
python scripts/plot_scdna_coverage.py \
  --data-dir /path/to/scdna_coverage_summaries
```

The script checks the expected 528-run composition, validates every raw
mosdepth summary against the sanitized per-run table, regenerates the panel,
and writes the group summary to `results/`.

The inside metric is the length-weighted mean depth across 105 deduplicated
SniffCell intervals (71,154 bp). The outside metric subtracts coverage in those
intervals from whole-reference coverage and uses the remaining GRCh38 span as
the denominator. No interval coordinates are released here.
