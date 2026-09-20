# 3.1 Candidate selection and paired-FANS validation

This folder derives counts from private SniffCell outputs and measures evidence
from indexed FANS BAMs. Only aggregate results are exported. Calls, BAMs, real
server paths, and local manifests are not public repository inputs.

## Run sections 3.1 and 3.2 together

Install `requirements.txt`. Create a private manifest using
`config/inputs.example.json`; paths resolve relative to the manifest. Supply
all six samples to reproduce the archived full-depth analysis. From this folder:

```bash
python scripts/run_analysis.py \
  --manifest private_inputs/local.inputs.json \
  --fans-annotations private_inputs/manual_fans_annotations.tsv \
  --motif-candidates private_inputs/full_depth_annotated_candidates.tsv \
  --outdir results/rebuilt
```

This executes candidate selection, manual FANS denominators, motif analysis in
section 3.2, and both plotting scripts. It starts from SniffCell outputs rather
than rerunning discovery. The checksummed versions reproduced here are archived
paper sources, not an automatically synchronized current manuscript spreadsheet.

## Candidate-selection rules

`scripts/rebuild_candidate_counts.py` adapts the original
`build_selfcontained_sv_germline_filtered_figures.py` and
`build_figure2E_top_panel_fig13_matched.py`:

- PBMC TRs require `n_change_support_reads >= 5`, `n_change_reads >= 5`, and
  `n_baseline_reads >= 10` in each run's `tr_changes.bed.tsv`.
- Full-depth brain TRs use recorded v2 artifact adjudication: 184 input rows,
  104 main keeps, 1 curated rescue, and 79 removals. Decisions are keyed by
  sample, TR locus, and cell type.
- Optional 30x/60x brain runs count exact recovery of those full-depth
  locus/cell-type keys, not all new subset discoveries.
- SVs require class SV, no hard conflict, assignment to only one comparison
  group, and at least 15 overlapping reads. Shared/germline SVs are excluded by
  ID or same chromosome/type, distance <= `max(1000 bp, 0.5 * larger length)`,
  and length difference <= 50%. Missing shared-SV files are errors.
- Output retains sample, coverage, variant class, and cell type. Invalid
  numeric values are excluded with audit counts. Missing columns, unknown cell
  labels, and duplicate adjudication keys fail explicitly.

Archived full-depth results are 105 brain TR rows (53 neuron, 52 oligodendrocyte)
and 4 SV rows (3 neuron, 1 oligodendrocyte), giving 56/53 combined candidates.
Candidate status does not imply experimental validation.

**Upstream gap:** the script deriving the v2 sequence metrics and decisions
was not located. This implementation reads the recorded adjudication; it does
not claim to recreate that investigation from BAMs. Keep/removal reasons and
the curated rescue remain explicit in the audit.

## Manual FANS validation and denominators

`scripts/summarize_fans_validation.py` calculates donor-locus rows and unique
genomic loci from per-row annotations. Unique loci use chromosome/start/end
within variant class and tissue scope. A locus is observed if any donor is Yes;
mixed Yes/No loci are counted explicitly. Unknown annotations are Unscored,
not No, and excluded from the eligible denominator.

The archived `Supplementary_Table_9_fig13_FANS_validation.tsv` contains brain
and PBMC, including subset-detected candidates. It is not the full-depth motif
table's candidate set.

| Scope | Observed / donor-locus rows | Observed / unique TR loci |
| --- | --- | --- |
| Brain and PBMC | 47 / 73 | 44 / 65 |
| Brain | 42 / 63 | 40 / 58 |
| PBMC | 5 / 10 | 5 / 9 |

47/73 versus 44/65 is a counting-unit difference, not a manuscript discrepancy.
The previous README incorrectly described a mismatch and the previous panel
labeled combined-cohort counts as brain. The regenerated panel uses brain-only
unique loci. `data/manuscript_reported_counts.tsv` is retained as a historical
prose snapshot, not analysis input or an assertion about the current paper.

## Direct BAM evidence

```bash
python scripts/measure_fans_bam_evidence.py \
  --candidates private_inputs/curated_candidates.tsv \
  --manifest private_inputs/local.fans.json \
  --outdir results/fans_bam_evidence
```

Use `config/fans.example.json` and explicitly matched donor IDs. The script
adapts the original `validate_curated_in_fans.py`:

- TR evidence is changed-minus-contrast median net CIGAR insertion length
  across reads spanning the locus plus 150-bp flanks. Each group needs 5 reads;
  the directional difference must be at least 5 bp.
- SV evidence counts indels within 150 bp of the breakpoint and length at least
  `max(20 bp, 0.5 * candidate length)`. Each group needs 5 spanning reads;
  changed VAF must be >= 0.15 and contrast VAF <= 0.05.
- Only primary mapped alignments are used. The historical proxy has no MAPQ
  threshold and does not require indel-type agreement. It is not a
  sequence-resolved genotype or a statistical significance test.
- Low coverage, missing donor matches, unmapped cell types, and rearrangements
  needing structural review are separate outcomes. Aggregated read counts
  represent read-locus observations, not globally unique reads.

This automated proxy cannot recreate manual inspection or medaka allele
concordance and never overwrites manual Yes/No labels. Only aggregate evidence
is exported.

## Required private input columns

| Input | Columns |
| --- | --- |
| Brain decisions | `sample, cell_type, trid, recommended_v2_keep_with_rescue, recommended_v2_reason` |
| TR changes | `trid, change_group, n_change_support_reads, n_change_reads, n_baseline_reads` |
| SV assignments | `id, variant_class, has_hard_conflict, linked_celltypes, primary_celltype, n_overlapped, sv_chr, sv_pos, sv_len, sv_type` |
| Shared SVs | VCF or VCF.gz, with `SVTYPE` and `SVLEN` INFO fields |
| Manual FANS annotations | `variant_type, tissue, donor, chrom, start, end, observed_in_paired_FANS` |
| BAM candidates | `sample, variant_type, cell_type, chrom, start, end, direction, size` |

BAM-validation coordinates are zero-based, half-open. Confirm the source
coordinate convention before supplying candidates. BAMs must be indexed and
use the same reference build and contig names. Section 3.2 documents motif input.

## Outputs and reproducibility limits

`results/rebuilt/counts_and_validation/` contains sample/class counts, sequential
filter exclusions, adjudication reasons, both FANS denominators, SHA256 input
checksums, parameters, and plots. Provenance uses logical roles, not paths.

Replot the bundled snapshot with `python scripts/plot_brain_counts_and_validation.py`.
That command alone is figure reproduction, not the upstream analysis. The
archived source PNG is retained. Run methodology tests with
`python -m unittest discover -s tests -v`.

The current public Zenodo deposition does not supply these private callsets or
adjudication inputs. Public readers can inspect the analysis and regenerate
the bundled plots, but cannot derive all candidates from the public deposition
alone. No variant calls are added to GitHub or Zenodo by these scripts.
