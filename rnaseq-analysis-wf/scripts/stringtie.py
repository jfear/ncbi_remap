"""Script to generate stringtie command."""
from pathlib import Path

from snakemake.shell import shell

from ncbi_remap.snakemake import get_flag

pth = Path(f'../output/prealn-wf/samples/{snakemake.wildcards.srx}')
glob = list(pth.glob('*/STRAND'))
strand = get_flag(glob[0])

if strand == 'samestrand':
    param = '--rf'
elif strand == 'oppositestrand':
    param = '--fr'
else:
    param = ''


cmd = f"""stringtie \
    {snakemake.input.bam} \
    -G {snakemake.input.gtf} \
    {param} \
    -p {snakemake.threads} \
    -o {snakemake.output[0]}
"""

print(cmd)
shell(cmd)
