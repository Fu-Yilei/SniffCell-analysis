import pathlib,json,shutil,sys
import pandas as pd
R=pathlib.Path(sys.argv[1]).resolve()
P=pathlib.Path(__file__).resolve().parents[1];D=P/'data'
text='''# A first held-out test of the 5hmC-aware SniffCell atlas

This quick validation holds out one control donor, bcontrol1. Both its neuronal and oligodendrocyte FANS samples are excluded from atlas discovery; bmsa, bparkinsons, and bcontrol2 provide the training data. This is a single-donor proof of concept, not a completed four-donor validation.

The sole baseline is the paper’s 152,961-region modifiedC catalog, `brain_cereb.maxgap3000_dt04.ctdmr.tsv`. The integrated atlas keeps this backbone and adds eligible training-only ONT 5mC and 5hmC markers. Earlier comparisons using a different backbone are superseded.

Sorted FANS sample identity provides the reference cell-type label. TP is a read assigned to its source population; FP is a read from the other population assigned to that type. FN includes wrong assignments, ambiguous reads, and reads with no evidence. Precision is TP/(TP+FP); recall is TP/(TP+FN). These labels inherit the purity limitations of population sorting. Actual bulk reads do not have per-read truth here; matched-FANS methylation concordance is a separate assessment, deferred in this quick run.

'''
complete=True
for title,folder,prefix in [('FGF14',R/'fixed_fgf14/summary','fixed_fgf14'),('Genome-wide',R/'summary/genome_controls','genomewide_control')]:
 text+='## '+title+'\n\n'
 if title=='FGF14':text+='Both atlases receive the same pooled FANS reads spanning the full interval chr13:102161539–102161869. Selection uses primary alignments only, without methylation or assignment filters. All catalog markers overlapping those read spans are supplied, with no atlas-dependent nearest-marker window.\n\n'
 else:text+='Both atlases use the same FANS BAMs and frozen unique primary-read identities. Reads without assignment evidence remain in the recall denominator.\n\n'
 if not (folder/'DONE').exists():text+='Results pending.\n\n';complete=False;continue
 for name in ['precision_recall.tsv','per_cell_type.tsv','confusion.tsv','selected_markers.tsv','legacy_parity.tsv','control_identity_checks.tsv']:
  p=folder/name
  if p.exists():shutil.copy2(p,D/(prefix+'_'+name))
 f=pd.read_csv(folder/'precision_recall.tsv',sep='\t');f=f[(f.donor=='ALL')&f['mode'].isin(['legacy','dual'])]
 text+='| Atlas | Control reads | Correct / assigned | Micro-precision | Micro-recall | No evidence | Ambiguous |\n|---|---:|---:|---:|---:|---:|---:|\n'
 for x in f.itertuples():text+=f'| {"Paper backbone" if x.mode=="legacy" else "Integrated dual"} | {x.control_reads} | {x.correct_reads} / {x.assigned_reads} | {100*x.micro_precision:.2f}% | {100*x.micro_recall:.2f}% | {x.no_evidence_reads} | {x.ambiguous_reads} |\n'
 text+='\n'
par=R/'fixed_fgf14/summary/legacy_parity.tsv'
text+='## Backward compatibility and remaining scope\n\n'
if par.exists():
 f=pd.read_csv(par,sep='\t');text+=f'The paper backbone passed exact comparison of all historical columns in {len(f)} FGF14 result tables against pre-5hmC code. Added provenance columns are excluded from that comparison. This regional check does not establish genome-wide parity.\n\n'
else:text+='Exact historical-column compatibility on the FGF14 control is pending.\n\n'
text+='Other donor folds, single-channel ablations, and bulk/FANS checks are deferred to keep this first test small. One held-out donor cannot establish performance across donors. The [v0.9.8a1 alpha](https://github.com/Fu-Yilei/SniffCell/releases/tag/v0.9.8a1) is now published on GitHub and PyPI; install it explicitly with `pip install sniffcell==0.9.8a1`.\n'
text+='\n## Benchmark plots\n\n![FGF14 precision and recall](results/fixed_fgf14_precision_recall.png)\n\n![Genome-wide precision and recall](results/genomewide_control_precision_recall.png)\n\n![Genome-wide read outcomes](results/genomewide_read_outcomes.png)\n\n'
text+='\n'+(P/'IMPLEMENTATION.md').read_text()+'\n'
(P/'BLOG_UPDATE.md').write_text(text.rstrip()+'\n')
(D/'status.json').write_text(json.dumps({'status':'quick_validation_complete' if complete else 'in_progress','heldout_donor':'bcontrol1','backbone':'brain_cereb.maxgap3000_dt04.ctdmr.tsv','deferred':['other donor folds','bulk benchmarking','single-channel ablations']},indent=2))

print(P/'BLOG_UPDATE.md')
