# 2.4 CALCOCO2 T-cell expansion

This analysis reproduces Figure 3c from a public, de-identified per-read table
containing only donor, source, cell label, and flank-corrected repeat length.
SniffCell T/B/NK labels are grouped as `Lymphocytes`; FACS fractions remain
separate so the expansion can be localized to T cells.

Run:

```bash
python scripts/plot_calcoco2_expansion.py
```

The exact manuscript-source PNG is retained as
`results/figure_3c_calcoco2_expansion_source.png`. The script writes its
rebuild, PDF, and series-level summary table to `results/`. There are no read
names, input paths, BAMs, or variant calls.
