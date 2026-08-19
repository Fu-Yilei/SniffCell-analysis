# 1.9 genomic-region coverage

This analysis reproduces the genomic-region coverage comparisons in Figure 2d
and Supplementary Figure 8. Coverage is the percentage of bases in each region
class covered at `>=5x` by SniffCell-assigned reads, averaged across samples and
cell-type groups within brain or PBMC.

Figure 2d uses all qualifying GTEx tissue-specific genes, not two example
genes per tissue:

- brain: 2,809 GTEx tissue-specific genes;
- blood/PBMC: 720 GTEx tissue-specific genes.

The GTEx tissue-specific bar is intentionally first (leftmost), followed by
GENCODE genes, medically relevant genes (MRG), cCREs, and the remaining genome.
Each panel prints its tissue-specific gene count directly under the first bar.

Run:

```bash
python scripts/plot_coverage_figures.py
```

The script uses only path-free aggregate coverage tables and writes Figure 2d
and Supplementary Figure 8 as PNG and PDF files under `results/`.
