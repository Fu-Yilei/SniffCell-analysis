# 5hmC atlas update

The [v0.9.8a1 alpha](https://github.com/Fu-Yilei/SniffCell/releases/tag/v0.9.8a1) is now published on GitHub and PyPI; install it explicitly with `pip install sniffcell==0.9.8a1`.

The [held-out validation blog](analyses/06_06_heldout_validation/BLOG_UPDATE.md)
compares the original paper backbone with a channel-aware integrated atlas
on identical bcontrol1 FANS reads. It includes bar plots, aggregate tables,
and the detailed implementation record. Precision improves from 87.01% to
88.64%, and recall from 18.33% to 31.03%, across the same 46,373,872 controls.
This is one held-out donor, with other folds and bulk benchmarking deferred.

The [four-donor atlas download](https://github.com/Fu-Yilei/SniffCell/tree/5hmC_compatitible_atlas/atlases/5hmc)
is a separate distributable build that includes bcontrol1; its training
membership differs from the benchmark atlas.

The [new visualization blog](analyses/06_07_assay_visualization/BLOG_UPDATE.md)
shows the bmsa FGF14 example with proportional m/h ctDMR boxes, native haplotype
filtering and an explicitly disclosed one-read display exclusion. It documents
the development-branch `viz`/`report` update; this is not a new benchmark.

# SniffCell paper analyses

This repository contains figure-reproduction analyses for the SniffCell paper.

Each analysis is organized under `analyses/` and includes its figure-ready
data, scripts, generated tables or figures, and provenance. Large inputs used
by more than one section live once under `analyses/shared_data/`.

Public repository rules:

- no server-specific file paths;
- no variant call files or variant-level call tables;
- no sequencing read identifiers;
- manuscript figures must be reproducible from the included public data.

Generated figures use a Nature-style layout: figure-level titles and bottom
caption/footnote text are omitted, while panel labels, facet headings, legends,
axes, and data annotations are retained. Run
`python scripts/audit_figure_style.py` to enforce this rule. Immutable
`*_source` assets are retained only for provenance and checksum verification.

Exact retained source-figure assets can be checked with
`sha256sum -c SOURCE_FIGURE_SHA256SUMS`.

Analyses are added one paper section at a time after review.

## Included analyses

- `analyses/01_01_ctdmr_parameter_selection/`: ctDMR parameter selection and
  Supplementary Figure 1.
- `analyses/01_02_eight_tissue_ctdmr_atlases/`: eight-tissue SniffCell `find`
  inputs, generated ctDMR atlases, and Supplementary Figure 2.
- `analyses/01_03_loyfer_atlas_overlap/`: overlap and directional concordance
  between the public Loyfer atlas and the eight SniffCell tissue atlases.
- `analyses/01_04_ctdmr_fans_benchmark/`: ctDMR atlas directional-effect
  benchmarking against matched FANS samples (Supplementary Figure 3).
- `analyses/01_05_per_cpg_methylation_concordance/`: pooled per-CpG concordance
  between SniffCell-assigned bulk reads and FANS (Supplementary Figure 4).
- `analyses/01_06_read_level_ctdmr_disagreement/`: read-level disagreement with
  FANS at ctDMRs (Supplementary Figure 5).
- `analyses/01_07_conflicting_brain_reads_monocyte_signal/`: Monocyte-like
  ctDMR signal among conflicting brain reads (Supplementary Figure 6).
- `analyses/01_08_fans_sniffcell_ctdmr_pca/`: PCA of matched FACS/FANS and
  SniffCell-assigned ctDMR methylation profiles (Supplementary Figure 7).
- `analyses/01_09_genomic_region_coverage/`: genomic-region coverage, including
  all GTEx tissue-specific genes (Figure 2d and Supplementary Figure 8).
- `analyses/01_10_gtex_expression_ctdmr_enrichment/`: GTEx expression and
  gene-body ctDMR enrichment (Supplementary Figure 9).
- `analyses/02_01_pbmc_aggregate_variant_counts/`: aggregate PBMC germline and
  cell-type-associated variant counts for Supplementary Tables 1-2 and the
  PBMC component of Figure 2f.
- `analyses/02_02_vdj_rearrangement/`: Figure 3a-b V(D)J deletion evidence,
  paired with sanitized raw event tables on Zenodo.
- `analyses/02_03_pbmc_structural_variant_example/`: checksum-verified
  Supplementary Figure 14 IGV asset and its provenance discrepancy.
- `analyses/02_04_calcoco2_t_cell_expansion/`: Figure 3c CALCOCO2 per-read
  repeat-length analysis in SniffCell and sorted PBMC fractions.
- `analyses/03_01_brain_variant_counts_fans_validation/`: brain candidate
  counts, paired-FANS validation totals, and manuscript reconciliation.
- `analyses/03_02_brain_tr_motif_composition/`: aggregate brain/PBMC TR motif
  length and sequence-composition figure.
- `analyses/03_03_ont_low_quality_snv_artifact/`: checksum-verified
  Supplementary Figure 12 FBXL7 low-base-quality IGV asset.
- `analyses/03_04_sgip1_neuron_expansion/`: Supplementary Figure 15 SGIP1
  per-read expansion analysis across four brain donors.
- `analyses/03_05_fgf14_neuron_expansion/`: Figure 3d FGF14 TRGT
  neuron-versus-oligodendrocyte distribution analysis.
- `analyses/03_06_fgf14_flanking_variant_interruptions/`: Extended Data
  Figure 1, 5'-CFV allele ranges, and interruption summaries.
- `analyses/03_07_scdna_coverage_validation/`: Extended Data Figure 2a-b scDNA
  coverage validation, paired with sanitized summaries and 528 mosdepth
  summary files on Zenodo.
- `analyses/04_01_mrb_tr_burden_recurrence/`: Figure 4a multi-region brain TR
  burden, recurrence, and AAG-like motif enrichment.
- `analyses/04_02_mrb_regional_prevalence/`: Figure 4b-c regional prevalence
  and donor-stratified permutation results.
- `analyses/04_03_mrb_regional_overlap/`: Figure 4d and Supplementary Figure 16
  pairwise brain-region overlap.
- `analyses/04_04_mrb_donor_region_distribution/`: Figure 4e donor-by-region
  distribution of globally cell-type-specific loci.
- `analyses/04_05_pja2_read_level_sizing/`: Figure 4f PJA2 read-level sizing,
  paired with sanitized point data on Zenodo.
- `analyses/04_06_matched_blood_tr_expansion/`: Figure 4g-h matched-blood
  comparison, paired with sanitized figure-source rows on Zenodo.
- `analyses/05_01_card_gaa_recurrence_replication/`: Figure 5 CARD GAA-family
  recurrence and directional cross-cohort replication from aggregate tables.
- `analyses/05_02_read_length_n50/`: Supplementary Figure 17 read-length N50
  comparison across the FANS, multi-region, and large-scale brain cohorts.
- `analyses/05_03_card_gaa_permutation/`: age-independent component of
  Supplementary Figure 18, using precomputed within-donor permutation results.
- `analyses/shared_data/modkit_cpg_pileups_01_04_to_01_06/`: one Zenodo-backed
  set of modkit-derived CpG pileups shared by analyses 1.4-1.6.
- `analyses/shared_data/read_assignments_01_07_to_01_08/`: the Zenodo-backed
  original per-read SniffCell assignment tables for six bulk and sixteen
  purified FACS/FANS samples.

Third-party data retain their upstream terms. In particular, the Loyfer UXM
atlas under analysis 1.3 is subject to the Software Research License included
in that analysis folder.
