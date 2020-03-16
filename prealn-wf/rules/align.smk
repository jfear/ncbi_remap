
rule hisat2:
    """Basic alignment."""
    input:
        index=config['references']['dmel']['hisat2'],
        splice_sites=config['references']['dmel']['known_splice_sites'],
        R1=rules.atropos.output.R1,
        R2=rules.atropos.output.R2,
        layout=rules.fastq_check.output.layout,
    output:
        bam=temp("../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam"),
        bai=temp("../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam.bai"),
    log: temp("../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam.log")
    params:
        hisat2_extra='--max-intronlen 300000 --known-splicesite-infile {splice} '.format(splice=config['references']['dmel']['known_splice_sites']),
        samtools_sort_extra='--threads 4 -l 9 -m 3G -T $TMPDIR/samtools_sort'
    threads: 8
    group: "hisat2"
    resources:
        mem_gb=lambda wildcards, attempt: attempt * 16,
        time_hr=lambda wildcards, attempt: attempt * 8
    script: "../scripts/hisat2.py"


rule hisat2_summary:
    input: 
        _=lambda wildcards: queue.expand(rules.hisat2.output.bam, wildcards.srr),
        log=lambda wildcards: queue.expand(rules.hisat2.log[0], wildcards.srr)
    output: "../output/prealn-wf/hisat2/{srr}.parquet"
    threads: 8
    group: "hisat2"
    resources:
        mem_gb=lambda wildcards, attempt: attempt * 16,
        time_hr=lambda wildcards, attempt: attempt * 8
    script: "../../scripts/hisat2_check.py"


rule run_stats:
    input:
        bam=rules.hisat2.output.bam,
        bai=rules.hisat2.output.bai,
    output:
        samtools_stats="../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam.samtools.stats",
        samtools_idxstats="../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam.samtools.idxstats",
        bamtools_stats="../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam.bamtools.stats",
    threads: 8
    group: "hisat2"
    resources:
        mem_gb=lambda wildcards, attempt: attempt * 16,
        time_hr=lambda wildcards, attempt: attempt * 8
    shell:
        'BAM=$(mktemp --suffix=".bam") '
        '&& cp {input.bam} $BAM '
        '&& cp {input.bam}.bai $BAM.bai '
        '&& samtools stats $BAM > {output.samtools_stats} '
        '&& samtools idxstats $BAM > {output.samtools_idxstats} '
        '&& bamtools stats -in $BAM > {output.bamtools_stats} '
        '&& rm $BAM'

