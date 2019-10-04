"""Script to run samtools merge on a large number of BAMs."""
import os
from time import sleep
from tempfile import NamedTemporaryFile

from snakemake.shell import shell

if os.getenv("SLURM_JOBID", False):
    TMPDIR = os.path.join('/lscratch', os.getenv('SLURM_JOBID'))
else:
    TMPDIR = os.getenv('TMPDIR', "/tmp")
shell.prefix("set -euo pipefail; export TMPDIR={};".format(TMPDIR))

threads = snakemake.threads

bams = list(snakemake.input)

tmp = NamedTemporaryFile(mode='w', dir=TMPDIR, delete=False)
tmp.write('\n'.join(bams))
tmp.close()

cmd = f"""samtools merge \
    -b {tmp.name} \
    -@ {threads} \
    {snakemake.output[0]}
"""

print(cmd)
shell(cmd)

os.unlink(tmp.name)

