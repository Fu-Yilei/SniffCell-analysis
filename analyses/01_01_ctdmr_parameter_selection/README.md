# 1.1 ctDMR parameter selection

This analysis reproduces Supplementary Figure 1 of the SniffCell paper. The
figure shows the relationship between SniffCell `find` parameters, methylation
similarity to matched FACS/FANS references, and the number of detected ctDMRs.

## Paper panels

- **A:** minimum CpGs per ctDMR (`2`, `3`, `4`, and `5`)
- **B:** maximum gap between CpGs (`250`, `500`, `1000`, and `2000` bp)

The bundled table contains the eight aggregate rows used by the manuscript
figure. It contains no variant calls, read identifiers, or private paths.

## Reproduce the figure

From this analysis directory:

```bash
python scripts/plot_supplementary_figure_1.py
```

The command writes:

- `results/supplementary_figure_1.png`
- `results/supplementary_figure_1.pdf`

Custom locations can be supplied with `--input` and `--output-dir`.

## Data columns

- `parameter`, `value`: parameter represented on the x-axis.
- `brain_similarity`: mean brain benchmark similarity.
- `pbmc_similarity`: mean PBMC benchmark similarity.
- `brain_ctdmr_count`, `pbmc_ctdmr_count`: mean detected-region counts.
- `total_ctdmr_count`: sum of the brain and PBMC mean counts.

## Interpretation note

The published trend points are grouped means across the tested parameter
configurations sharing a given `min_cpgs` or `max_gap_bp` value. Other
parameters were not fixed across every point, so the figure describes the
observed parameter sweep rather than a strict one-variable-at-a-time
experiment.

## Provenance

The values were taken from the finalized SniffCell `find` benchmarking summary
used for the manuscript. The original workflow ran SniffCell `find`, compared
the resulting ctDMRs with matched FACS/FANS methylation references, and then
aggregated the benchmark metrics for plotting. Only the path-free, figure-ready
aggregate table is included here.
