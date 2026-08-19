# 1.4 ctDMR atlas benchmarking against FANS

This analysis reproduces Supplementary Figure 3, benchmarking the directional
methylation effects in the SniffCell ctDMR atlases against matched FANS
samples. The six panels are PBMC1, PBMC2, bcontrol1, bcontrol2, bmsa, and
bparkinsons.

The exact inputs are shared with analyses 1.5 and 1.6:

- FANS modkit-derived CpG pileups from
  `../shared_data/modkit_cpg_pileups_01_04_to_01_06/`;
- the PBMC `combined_dt04` atlas downloaded by analysis 1.2;
- `data/brain_cereb.maxgap3000_dt04.ctdmr.tsv`, with an identical copy stored
  in Zenodo record 10.5281/zenodo.22003085.

The committed `data/supplementary_figure_3_panel_data.npz` is a compact,
path-free regional cache for exact figure rendering. It contains the catalog
and FANS directional effects, not variant calls or read identifiers.

Render the manuscript figure:

```bash
python scripts/plot_supplementary_figure_3.py
```

The full pileup-to-panel cache builder is `scripts/build_panel_cache.py`. It
requires `bedtools`, the shared Zenodo pileups, the analysis-1.2 PBMC atlas,
and the committed brain benchmark atlas. `data/atlas_provenance.tsv` records
the exact row counts and checksums used by the figure.
