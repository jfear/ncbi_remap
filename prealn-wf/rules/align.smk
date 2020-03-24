THREADS = 8
GROUP = "hisat2"
RESOURCES = dict(
    mem_gb=lambda wildcards, attempt: attempt * 16,
    time_hr=lambda wildcards, attempt: attempt * 8
)

rule hisat2:
    """Basic alignment."""
    input:
        index=config['references']['dmel']['hisat2'],
        splice_sites=config['references']['dmel']['known_splice_sites'],
        R1=rules.atropos.output.R1,
        R2=rules.atropos.output.R2,
        layout=rules.fastq_check.output.layout,
        _=rules.atropos_check.output[0]
    output:
        bam=temp("../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam"),
        bai=temp("../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam.bai"),
    log: "../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam.log"
    params:
        hisat2_extra='--max-intronlen 300000 --known-splicesite-infile {splice} '.format(splice=config['references']['dmel']['known_splice_sites']),
        samtools_sort_extra='--threads 4 -l 9 -m 3G -T $TMPDIR/samtools_sort'
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../../scripts/hisat2.py"


rule hisat2_check:
    input: 
        bam=lambda wildcards: queue.expand(rules.hisat2.output.bam, wildcards.srr),
        bai=lambda wildcards: queue.expand(rules.hisat2.output.bai, wildcards.srr),
        log=lambda wildcards: queue.expand(rules.hisat2.log[0], wildcards.srr)
    output: "../output/prealn-wf/hisat2/{srr}.parquet"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../../scripts/hisat2_check.py"


rule aln_stats:
    input:
        bam=rules.hisat2.output.bam,
        bai=rules.hisat2.output.bai,
        _=rules.hisat2_check.output[0]
    output:
        samtools_stats=temp("../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam.samtools.stats"),
        bamtools_stats=temp("../output/prealn-wf/samples/{srx}/{srr}/{srr}.hisat2.bam.bamtools.stats"),
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    shell:
        'BAM=$(mktemp --suffix=".bam") '
        '&& cp {input.bam} $BAM '
        '&& cp {input.bam}.bai $BAM.bai '
        '&& samtools stats $BAM > {output.samtools_stats} '
        '&& bamtools stats -in $BAM > {output.bamtools_stats} '
        '&& rm $BAM'


rule aln_stats_summary:
    input:
        samtools=lambda wildcards: queue.expand(rules.aln_stats.output.samtools_stats, wildcards.srr),
        bamtools=lambda wildcards: queue.expand(rules.aln_stats.output.bamtools_stats, wildcards.srr)
    output: "../output/prealn-wf/aln_stats/{srr}.parquet"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../../scripts/aln_stats_summary.py"
