# 3.6 FGF14 flanking variant and interruption analysis

This section supports Extended Data Figure 1 and the FGF14 sequence-context
analysis. It contains:

- the exact original CFV-stratified figure as
  `results/extended_data_figure_1_source.png`;
- an aggregate table for the four benchmark brain donors, with read-length
  ranges and counts split by 5'-CFV and non-CFV allele; and
- consensus/FANS interruption summaries with no allele sequence, read name,
  VCF record, or server path.

Run:

```bash
python scripts/summarize_fgf14_sequence_context.py
```

The script verifies the exact figure asset, renders an independently runnable
range summary, and writes a compact interruption table. The 5'-CFV is treated
as a flanking polymorphism, not as a GAA-tract interruption.
