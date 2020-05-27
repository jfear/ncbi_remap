"""Rules to download FASTQ from the SRA."""
rule sra_prefetch:
    """Downloads container file from the SRA"""
    output: "../output/fastq-wf/sra_cache/{srr}.sra"
    log: "../output/fastq-wf/sra_download_logs/{srr}.log"
    resources:
        mem_gb=lambda wildcards, attempt: attempt * 2,
        time_hr=lambda wildcards, attempt: attempt * 8
    shell: "prefetch --output-file {output[0]} --max-size 40G {wildcards.srr} > {log} 2>&1"


rule fastq_dump:
    input: 
        sra=rules.sra_prefetch.output[0],
        log=rules.sra_prefetch.log[0]
    output: 
        r1=temp("../output/fastq-wf/fastqs/{srr}_1.fastq.gz"),
        r2=temp("../output/fastq-wf/fastqs/{srr}_2.fastq.gz"),
        layout="../output/fastq-wf/layout/{srr}.parquet",
        summary="../output/fastq-wf/libsize/{srr}.parquet"
    params:
        download_bad="../output/fastq-wf/download_bad",
        abi_solid="../output/fastq-wf/abi_solid",
    threads: 8
    resources:
        mem_gb=lambda wildcards, attempt: attempt * 4,
        time_hr=lambda wildcards, attempt: attempt * 12
    script: "../scripts/fastq_extract.py"
