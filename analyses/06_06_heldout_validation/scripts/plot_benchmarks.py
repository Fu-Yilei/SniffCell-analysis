"""Plot verified fixed-control aggregates; no individual read data required."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import PercentFormatter
P=Path(__file__).resolve().parents[1]; D=P/'data'; O=P/'results';O.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.dpi':220})
colors=['#687D90','#007F79'];names=['Paper backbone','Integrated dual']
def save(fig,name):
 for ext in ['png','pdf']:fig.savefig(O/f'{name}.{ext}',bbox_inches='tight',facecolor='white')
 plt.close(fig)
for prefix,title in [('fixed_fgf14','FGF14'),('genomewide_control','Genome-wide')]:
 f=pd.read_csv(D/f'{prefix}_precision_recall.tsv',sep='\t');f=f[f.donor=='ALL'].set_index('mode').loc[['legacy','dual']]
 assert f.control_reads.nunique()==1
 fig,axs=plt.subplots(1,2,figsize=(9,4.8))
 for ax,key,label in zip(axs,['micro_precision','micro_recall'],['Precision · correct / assigned','Recall · correct / all controls']):
  vals=100*f[key];bars=ax.bar(names,vals,color=colors,width=.58)
  ax.bar_label(bars,labels=[f'{v:.2f}%' for v in vals],padding=5,fontweight='bold');ax.set_ylim(0,105);ax.yaxis.set_major_formatter(PercentFormatter());ax.set_title(label,fontsize=12,pad=16);ax.set_axisbelow(True);ax.grid(axis='y',alpha=.18)
 fig.suptitle(f'{title}: one held-out donor (bcontrol1)',fontsize=16,fontweight='bold',y=1.03)
 fig.text(.5,.015,f'Identical {int(f.control_reads.iloc[0]):,} FANS control reads · both sorted populations\nNo-evidence and ambiguous reads remain in the recall denominator.',ha='center',fontsize=10,color='#45515B')
 fig.tight_layout(rect=(0,.13,1,.96));save(fig,prefix+'_precision_recall')
f=pd.read_csv(D/'genomewide_control_precision_recall.tsv',sep='\t');f=f[f.donor=='ALL'].set_index('mode').loc[['legacy','dual']]
fig,ax=plt.subplots(figsize=(10,4.7));left=[0.,0.]
for key,label,color in [('correct_reads','Correct','#007F79'),('wrong_reads','Wrong','#C45748'),('ambiguous_reads','Ambiguous','#DBAE55'),('no_evidence_reads','No evidence','#DDE3E8')]:
 vals=100*f[key]/f.control_reads;bars=ax.barh(names,vals,left=left,color=color,height=.5,label=label)
 for i,(v,n) in enumerate(zip(vals,f[key])):
  if v>=5:ax.text(left[i]+v/2,i,f'{v:.1f}%\n{n/1e6:.2f}M',ha='center',va='center',color='white' if key=='correct_reads' else '#26343C',fontsize=10)
 left=[a+b for a,b in zip(left,vals)]
ax.set_xlim(0,100);ax.xaxis.set_major_formatter(PercentFormatter());ax.invert_yaxis();ax.set_xlabel('Share of the same 46,373,872 FANS control reads');ax.set_title('Genome-wide assignment outcomes',loc='left',fontsize=16,fontweight='bold',pad=18);ax.legend(ncol=4,loc='upper center',bbox_to_anchor=(.5,-.2),frameon=False);fig.tight_layout();save(fig,'genomewide_read_outcomes')
print('Wrote three PNG/PDF figure pairs')
