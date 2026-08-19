# 3.2 TR motif composition

This section reproduces the motif-length and sequence-composition panels from
the aggregate Figure 13 TR table. The input contains category counts only, not
TR loci or variant calls.

The exact manuscript-source PNG is retained as
`results/tr_motif_composition_source.png`. Run:

```bash
python scripts/plot_tr_motif_composition.py
```

The source totals are 10 PBMC and 105 brain TR candidates. The script validates
that every dimension sums to those totals before rendering the PNG and PDF.
