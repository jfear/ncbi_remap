"""Rules to download FASTQ from the SRA."""
THREADS = 8
GROUP = "fastq-dump"
RESOURCES = dict(
    mem_gb=lambda wildcards, attempt: attempt * 4,
    time_hr=lambda wildcards, attempt: attempt * 12
)

rule sra_prefetch:
    """Downloads container file from the SRA"""
    output: "../output/fastq-wf/sra_cache/{srr}.sra"
    log: "../output/fastq-wf/sra_download_logs/{srr}.log"
    resources:
        mem_gb=lambda wildcards, attempt: attempt * 2,
        time_hr=lambda wildcards, attempt: attempt * 8
    shell: "prefetch --output-file {output[0]} --max-size 40G {wildcards.srr} > {log} 2>&1"


rule fastq_dump:
    input: rules.sra_prefetch.output[0]
    output: 
        r1=temp("../output/fastq-wf/fastqs/{srr}_1.fastq"),
        r2=temp("../output/fastq-wf/fastqs/{srr}_2.fastq"),
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../scripts/fastq_extract.py"


rule fastq_check:
    input: 
        sra=rules.sra_prefetch.output[0],
        r1=rules.fastq_dump.output.r1,
        r2=rules.fastq_dump.output.r2,
        download_log=rules.sra_prefetch.log[0]
    output: 
        layout="../output/fastq-wf/fastq_info/{srr}/LAYOUT",
        summary="../output/fastq-wf/fastq_info/{srr}/summary.tsv"
    params:
        download_bad="../output/fastq-wf/download_bad",
        abi_solid="../output/fastq-wf/abi_solid",
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../scripts/fastq_check.py"
