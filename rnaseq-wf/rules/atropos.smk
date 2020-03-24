"""Run Atropos to Trim Reads"""
THREADS = 8
GROUP = "atropos"
RESOURCES = dict(
    mem_gb=lambda wildcards, attempt: attempt * 8,
    time_hr=lambda wildcards, attempt: attempt * 12
)

rule atropos:
    """Filter reads that are less than 25bp."""
    input:
        R1=rules.fastq_dump.output.r1,
        R2=rules.fastq_dump.output.r2,
        layout=rules.fastq_check.output.layout
    output:
        R1=temp("../output/rnaseq-wf/samples/{srx}/{srr}/{srr}_1.trim.clean.fastq"),
        R2=temp("../output/rnaseq-wf/samples/{srx}/{srr}/{srr}_2.trim.clean.fastq")
    params:
        extra_pe='-a file:../data/adapters.fa -A file:../data/adapters.fa -q 20 --minimum-length 25',
        extra_se='-a file:../data/adapters.fa -q 20 --minimum-length 25',
    log: "../output/rnaseq-wf/samples/{srx}/{srr}/{srr}_1.trim.clean.fastq.gz.log"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../../scripts/atropos.py"


rule atropos_check:
    input: 
        R1=lambda wildcards: queue.expand(rules.atropos.output.R1, wildcards.srr),
        R2=lambda wildcards: queue.expand(rules.atropos.output.R2, wildcards.srr),
        log=lambda wildcards: queue.expand(rules.atropos.log[0], wildcards.srr)
    output: "../output/rnaseq-wf/atropos/{srr}.parquet"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../../scripts/atropos_check.py"
