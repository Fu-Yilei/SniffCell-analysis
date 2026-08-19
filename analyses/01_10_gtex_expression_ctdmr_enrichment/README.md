# 1.10 GTEx expression and ctDMR enrichment

This analysis reproduces the claims and panels in Supplementary Figure 9:
expressed genes are more likely to contain a gene-body ctDMR, ctDMR-density
distributions are shown among genes with ctDMRs, and ctDMR density is compared
with GTEx expression in brain and whole blood.

The committed gene-level table contains GTEx expression and ctDMR-density
summaries for 54,626 genes in each tissue comparison. It contains no sample
paths, read identifiers, BAMs, or variant calls.

Run:

```bash
python scripts/plot_supplementary_figure_9.py
```

Outputs are written to `results/supplementary_figure_9.png` and `.pdf`.
