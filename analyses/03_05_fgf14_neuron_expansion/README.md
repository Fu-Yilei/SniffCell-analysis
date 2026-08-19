# 3.5 FGF14 neuron-associated expansion

This section reproduces Figure 3d using the current SniffCell v0.9.4 TRGT
output. The public per-read table excludes read names and server paths while
retaining the numeric allele-length and motif estimates needed for the plot.

The exact manuscript-source PNG is retained as
`results/figure_3d_fgf14_neuron_expansion_source.png`. Run:

```bash
python scripts/plot_fgf14_neuron_expansion.py
```

Points are distribution-quantile matches between neuron and oligodendrocyte
reads within donor and allele; they are not matched biological reads. The
estimated repeat count is anchored to each assigned TRGT consensus allele:

`estimated GAA = consensus motif count + (read bp - consensus bp) / 3`.

No VCF, consensus sequence, or read identifier is included in this section.
