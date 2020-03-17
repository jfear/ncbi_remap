"""Run Feature Counts"""
THREADS = 4
GROUP = "feature_counts"
RESOURCES = dict(
    mem_gb=lambda wildcards, attempt: attempt * 4,
    time_hr=lambda wildcards, attempt: attempt * 4
)

rule feature_counts:
    input:
        annotation=config['references']['dmel']['gtf'],
        bam=rules.hisat2.output.bam,
        layout=rules.fastq_check.output.layout,
        strand=rules.collectrnaseqmetrics_summary.output.flag,
    output:
        counts=temp("../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam.feature_counts.counts"),
        jcounts=temp("../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam.feature_counts.counts.jcounts"),
        summary=temp("../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam.feature_counts.counts.summary"),
    log: "../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam.feature_counts.counts.log"
    params:
        extra_pe='-p -P -C -J ',
        extra_se='-J '
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../../scripts/featurecounts.py"


rule featurecounts_summary:
    input:
        counts=lambda wildcards: queue.expand(rules.feature_counts.output.counts, wildcards.srr),
        jcounts=lambda wildcards: queue.expand(rules.feature_counts.output.jcounts, wildcards.srr)
    output: "../output/prealn-wf/count_summary/{srr}.parqeut",
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../scripts/featurecounts_summary.py"
    