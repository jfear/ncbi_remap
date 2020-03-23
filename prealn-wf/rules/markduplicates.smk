"""Run Markduplicates."""
THREADS = 2
GROUP = "markduplicates"
RESOURCES = dict(
    mem_gb=lambda wildcards, attempt: attempt * 20,
    time_hr=lambda wildcards, attempt: attempt * 12
)

rule markduplicates:
    input:
        bam=rules.hisat2.output.bam,
        _=rules.hisat2_check.output[0]
    output:
        bam=temp("../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam.picard.markduplicates.bam"),
        metrics=temp("../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam.picard.markduplicates.metrics"),
    log: "../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam.picard.markduplicates.metrics.log"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../scripts/markduplicates.py"


rule markduplicates_summary:
    input: 
        lambda wildcards: queue.expand(rules.markduplicates.output.metrics, wildcards.srr)
    output: "../output/prealn-wf/markduplicates/{srr}.parquet"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../scripts/markduplicates_summary.py"
