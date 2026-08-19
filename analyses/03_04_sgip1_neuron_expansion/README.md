# 3.4 SGIP1 neuron-associated expansion

This section reproduces Supplementary Figure 15 for the SGIP1 repeat at
`chr1:66630776-66631127`. The public input retains only donor, source, cell
type, and flank-corrected read length; read identifiers and server paths were
removed.

The exact manuscript-source PNG is retained as
`results/supplementary_figure_15_sgip1_source.png`. Run:

```bash
python scripts/plot_sgip1_expansion.py
```

The output compares SniffCell-assigned and independently sorted FANS neuron
and oligodendrocyte reads across four donors and writes a numerical summary.
