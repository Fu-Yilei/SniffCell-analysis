"""Reproduce the selected illustration using authorized local inputs.

The exclusion file contains the single borderline support read discussed in
the blog. It is supplied privately; identifiers are not distributed here.
Input BAM: original full alignments of all 98 primary repeat-spanning reads.
"""
import argparse
from pathlib import Path
import subprocess
import sys

parser=argparse.ArgumentParser(description=__doc__)
for name in ['bam','variants','reference','atlas','read-assignment','read-summary','exclude-read-names']:
    parser.add_argument('--'+name,required=True)
parser.add_argument('--output-dir',type=Path,default=Path('results'))
a=parser.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
for ext in ['png','pdf']:
    subprocess.run([sys.executable,'-m','sniffcell.main','viz',
        '-i',a.bam,'-v',a.variants,'-r',a.reference,'-b',a.atlas,
        '--sv_id','chr13_102161539_102161869','--read_assignment',a.read_assignment,
        '--read_summary',a.read_summary,'--window','8000','--exact_window',
        '--max_reads','1000','--linked_ctdmr_mode','strict',
        '--exclude_read_names',a.exclude_read_names,'--sample_label','bmsa bulk · FGF14',
        '--dpi','180','--format',ext,'--export_tables',
        '-o',str(a.output_dir/f'FGF14_bmsa_5hmc_viz.{ext}')],check=True)
# Generated TSV exports contain private identifiers; do not publish them.
