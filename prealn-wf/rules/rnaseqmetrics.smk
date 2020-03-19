"""Run Picard CollectRNASeqMetrics to determine strandedness."""
THREADS = 2
GROUP = "CollectRNASeqMetrics"
RESOURCES = dict(
    mem_gb=lambda wildcards, attempt: attempt * 12,
    time_hr=lambda wildcards, attempt: attempt * 12
)

rule collectrnaseqmetrics_unstrand:
    input:
        bam=rules.hisat2.output.bam,
        refflat=config['references']['dmel']['refflat'],
        _=rules.hisat2_check.output[0]
    output: temp("../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam.NONE.picard.collectrnaseqmetrics")
    params:
        extra='STRAND=NONE',
    log: "../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam.NONE.picard.collectrnaseqmetrics.log"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../scripts/picard_rnaseqmetrics.py"


rule collectrnaseqmetrics_first:
    input:
        bam=rules.hisat2.output.bam,
        refflat=config['references']['dmel']['refflat'],
        _=rules.hisat2_check.output[0]
    output: temp("../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam.FIRST_READ_TRANSCRIPTION_STRAND.picard.collectrnaseqmetrics")
    params:
        extra='STRAND=FIRST_READ_TRANSCRIPTION_STRAND',
    log: "../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam.FIRST_READ_TRANSCRIPTION_STRAND.picard.collectrnaseqmetrics.log"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../scripts/picard_rnaseqmetrics.py"


rule collectrnaseqmetrics_second:
    input:
        bam=rules.hisat2.output.bam,
        refflat=config['references']['dmel']['refflat'],
        _=rules.hisat2_check.output[0]
    output: temp("../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam.SECOND_READ_TRANSCRIPTION_STRAND.picard.collectrnaseqmetrics")
    params:
        extra='STRAND=SECOND_READ_TRANSCRIPTION_STRAND'
    log: "../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam.SECOND_READ_TRANSCRIPTION_STRAND.picard.collectrnaseqmetrics.log"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../scripts/picard_rnaseqmetrics.py"


rule collectrnaseqmetrics_summary:
    input:
        unstranded=lambda wildcards: queue.expand(rules.collectrnaseqmetrics_unstrand.output[0], wildcards.srr),
        first=lambda wildcards: queue.expand(rules.collectrnaseqmetrics_first.output[0], wildcards.srr),
        second=lambda wildcards: queue.expand(rules.collectrnaseqmetrics_second.output[0], wildcards.srr)
    output: 
        flag="../output/prealn-wf/strand/{srr}",
        table="../output/prealn-wf/rnaseqmetrics/{srr}.parquet",
        gene_coverage="../output/prealn-wf/gene_coverage/{srr}.parquet"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../scripts/rnaseqmetrics_summary.py"
