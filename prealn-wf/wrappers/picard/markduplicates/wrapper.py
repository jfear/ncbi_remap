from snakemake.shell import shell
log = snakemake.log_fmt_shell()
extra = snakemake.params.get('extra', '')
java_args = snakemake.params.get('java_args', '')
mem = int(snakemake.resources.get('mem_gb', ''))

if mem >= 10:
    mem = '-Xmx{}g'.format(mem - 2)
else:
    mem = '-Xmx{}g'.format(2)

shell(
    'picard '
    '{mem} '
    '{java_args} '
    'MarkDuplicates '
    'INPUT={snakemake.input.bam} '
    'OUTPUT={snakemake.output.bam} '
    'METRICS_FILE={snakemake.output.metrics} '
    '{extra} '
    '{log} '
)
