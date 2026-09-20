# 3.1 Per-locus candidates, coverage membership, and FANS validation

This analysis generates Figure 2e-g from the individual candidate records in
[Supplementary Table 5](https://docs.google.com/spreadsheets/d/1o31UEyJkIofEj7v_beWsKueVgT5R-KwZPNxyhfV1ubw/edit#gid=1888934510).
Source snapshot: 2026-09-20; 123 candidate/donor records.

## Figures

- [Figure 2e: counts by donor and coverage](results/table5/figure_2e_variant_counts.svg)
- [Figure 2f: paired-FANS validation](results/table5/figure_2f_fans_validation.svg)
- [Figure 2g: cell-type composition](results/table5/figure_2g_celltype_composition.svg)

PDF and PNG versions are beside each SVG. These retain the paper's six-facet
stacked count panel and horizontal percentage-bar panels, colors, and ordering.
There is no global title, footer, or panel letter; sample facet labels identify
the six donors. Axis labels are not bold. SVG text is editable and uses
DejaVu Sans, matching the original panel PDFs; install this font before editing
the SVGs to avoid font substitution. The PDFs embed the fonts.

## Reproduce from individual locus records

Run from this analysis folder after installing `requirements.txt`:

```bash
python scripts/fetch_table5.py \
  --out private_inputs/supplementary_table5.csv \
  --expected-record-sha256 5eefc604313077adf3be7d4685d019d09fa1c3431f9534203755181e4ea1a0c8

python scripts/run_analysis.py \
  --candidates private_inputs/supplementary_table5.csv \
  --outdir results/table5
```

The downloader reads the public supplement; it does not upload it. It refuses to
overwrite an existing snapshot. The record checksum is independent of CSV/JSON
format and row order. If the source changes, review the new per-locus records
before updating the snapshot checksum.

To run the GAA analysis in section 3.2 using the same candidate records:

```bash
python scripts/run_analysis.py \
  --candidates private_inputs/supplementary_table5.csv \
  --outdir results/table5 \
  --motif-outdir ../03_02_brain_tr_motif_composition/results/table5
```

Section 3.2 uses the motif and interval supplied in Table 5 by default.
Optionally add `--catalog private_inputs/adotto.v2.trgt.bed` to match motifs
from the Adotto reference, which is supplied separately.

The execution path is:

1. Read every candidate row, normalize whitespace and integer coordinates, and
   validate IDs, donor/locus/cell keys, cell labels, coverage membership, and
   explicit Yes/No annotations. Empty spreadsheet rows are ignored. Unscored
   annotations, duplicate records, and unknown labels are errors, not zeros.
2. Write a private normalized candidate table and a per-record panel-membership
   ledger. Every bar can be traced back to candidate IDs and source rows.
3. Calculate coverage counts, FANS numerators/denominators, cell-type totals,
   and the separate union of unique validated brain TR loci.
4. Plot the derived aggregate tables.
5. Run the GAA analysis from the same per-locus snapshot when `--motif-outdir`
   is supplied, optionally rematching intervals to the Adotto catalog.

## Counting definitions and results

| Output | Counting unit and selection | Result |
| --- | --- | --- |
| Figure 2e | Candidate/donor records with the requested depth in `coverages_detected`; stacked SV/TR | Full-depth PBMC 13; brain 110 |
| Figure 2f | All listed candidate/donor rows, grouped by variant class; manual FANS Yes/No | TR 98/115 observed; SV 4/8 observed |
| Figure 2g | All listed TR and SV rows, grouped by tissue and cell type | PBMC 11 lymphocyte / 2 monocyte; brain 59 neuron / 51 oligodendrocyte |
| Brain TR validation sentence | Unique chromosome/start/end intervals across donors and cell types; any donor Yes validates that locus | 83/92 |
| Brain TR classes | Candidate/donor rows, TR only | 54 neuron / 51 oligodendrocyte |
| Brain SV classes | Candidate/donor rows, SV only | 5 neuron / 0 oligodendrocyte |

Figure 2e shows depth membership of the Table 5 candidate set. An absent depth
token means that the listed candidate is not recorded as detected there.
These counts summarize the listed candidates rather than independent discovery
totals at each depth.

FANS Yes/No are the supplement's manual, matched-donor validation decisions.
For TRs the paper describes checking the same cell-type direction in per-read
length distributions. Figure 2f counts donor records; the 83/92 result counts
unique brain loci.

## Detailed inputs, outputs, and privacy

Required candidate columns:
`variant_id, variant_type, tissue, donor, chrom, start, end, cell_type, coverages_detected, observed_in_paired_FANS`.
A TSV/CSV may have one title row before its header; native spreadsheet connector
JSON with a `values` matrix is also accepted. Intervals must be valid integer
coordinates, and catalog analysis uses the source's GRCh38 half-open intervals.

`results/table5/` contains public aggregate tables, figures, and path-free
source URLs/checksums. `results/table5/private_inputs/` contains detailed records
and membership ledgers and is git-ignored. Do not publish these files or VCFs.
The downloader's `private_inputs/` is also ignored. No age metadata or private
server paths are included in public output.

The public supplement permits both table-to-figure and GAA analysis. A matching
Adotto reference is needed only for optional catalog matching. Manual FANS
annotations are read from the supplement, not inferred automatically from BAMs.

## Additional candidate-filtering analysis

`scripts/derive_brain_tr_candidates.py` implements sequence-artifact filtering,
cross-donor and depth reproducibility, and explicitly specified curated rescues
from per-locus sequence evidence. This optional upstream analysis uses the
private inputs described in its manifest, separately from the Table 5 workflow.

Copy `config/brain_filter.example.json` to `local.brain_filter.json` in this
analysis folder and populate the private input locations; paths resolve relative
to the manifest. Run:

```bash
python scripts/derive_brain_tr_candidates.py \
  --manifest local.brain_filter.json --outdir results/historical_filter
```

Its private manifest supplies `brain_sequence_metrics`, twelve
`brain_reproducibility_runs` (four donors at full/60x/30x), and explicitly
reasoned `curated_rescues`. It computes low-complexity scores from GC, entropy,
top-4-mer fraction and dinucleotide purity; support/effect-size selection;
cross-donor and depth reproducibility; residual low-complexity exclusion; and
the rescue. Detailed decisions are written only below `private_inputs/`.

The supplementary `rebuild_candidate_counts.py`, `summarize_fans_validation.py`,
`measure_fans_bam_evidence.py`, bundled `data/`, and results outside
`results/table5/` provide archived upstream analyses. Use the Table 5
workflow above to generate Figure 2e-g.

Run tests with `python -m unittest discover -s tests -v`.
