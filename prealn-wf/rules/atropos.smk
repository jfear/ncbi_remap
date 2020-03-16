"""Run Atropos to Trim Reads"""
rule atropos:
    """Filter reads that are less than 25bp."""
    input:
        R1=rules.fastq_dump.output.r1,
        R2=rules.fastq_dump.output.r2,
        layout=rules.fastq_check.output.layout
    output:
        R1=temp("../output/prealn-wf/samples/{srx}/{srr}/{srr}_1.trim.clean.fastq"),
        R2=temp("../output/prealn-wf/samples/{srx}/{srr}/{srr}_2.trim.clean.fastq"),
    params:
        extra_pe='-U 0 --minimum-length 25',
        extra_se='--minimum-length 25',
    log: temp("../output/prealn-wf/samples/{srx}/{srr}/{srr}_1.trim.clean.fastq.gz.log")
    threads: 8
    group: "atropos"
    resources:
        mem_gb=lambda wildcards, attempt: attempt * 8,
        time_hr=lambda wildcards, attempt: attempt * 12
    script: "../../scripts/atropos.py"


rule atropos_summary:
    input: 
        _=lambda wildcards: queue.expand(rules.atropos.output.R1, wildcards.srr),
        log=lambda wildcards: queue.expand(rules.atropos.log[0], wildcards.srr)
    output: "../output/prealn-wf/atropos/{srr}.parquet"
    threads: 8
    group: "atropos"
    resources:
        mem_gb=lambda wildcards, attempt: attempt * 8,
        time_hr=lambda wildcards, attempt: attempt * 12
    script: "../../scripts/atropos_check.py"
