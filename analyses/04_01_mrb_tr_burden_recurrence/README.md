# Multi-region brain TR burden and recurrence

This analysis reproduces the aggregate results underlying Figure 4a. The
included tables contain counts only, not variant calls or per-sample calls.

Key values are 56 neuron-expanded and 74 oligodendrocyte-expanded directional
events. Of the 25 neuron-exclusive and 43 oligodendrocyte-exclusive loci, 10
and 17 are recurrent, respectively. Thus the current oligodendrocyte value is
17/43 (39.53%); the 32.55% value in an earlier manuscript draft is stale.

Run:

```bash
python scripts/plot_burden_recurrence.py
```

`results/figure_4a_source.png` is the exact retained source panel.
