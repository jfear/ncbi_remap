"""Script to run samtools merge on a large number of BAMs."""
import os
from tempfile import NamedTemporaryFile

from snakemake.shell import shell

TMPDIR = os.path.join('/lscratch', os.getenv('SLURM_JOBID'))
shell.prefix("set -euo pipefail; export TMPDIR={};".format(TMPDIR))

bams = list(snakemake.input)
tmp = NamedTemporaryFile(mode='w', dir=TMPDIR)
tmp.write('\n'.join(bams))

cmd = f"""samtools merge \
    -b {tmp.name} \
    {snakemake.output[0]}
"""

print(cmd)
shell(cmd)

tmp.close()
