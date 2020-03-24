THREADS = 8
GROUP = "merge"
RESOURCES = dict(
    mem_gb=lambda wildcards, attempt: attempt * 8,
    time_hr=lambda wildcards, attempt: attempt * 4
)


rule bamMerge:
    input: 
        bams=lambda wildcards: queue.expand(rules.hisat2.output.bam, wildcards.srx),
        _=lambda wildcards: queue.expand(rules.hisat2_check.output[0], wildcards.srx),
    output:
        bam="../output/rnaseq-wf/samples/{srx}/{srx}.bam",
        bai="../output/rnaseq-wf/samples/{srx}/{srx}.bam.bai"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    shell:
        'samtools merge '
        '-f '
        '-@ {threads} '
        '{output.bam} '
        '{input.bams} && '
        'samtools index {output.bam}'


rule aln_stats:
    input:
        bam=rules.bamMerge.output.bam,
        bai=rules.bamMerge.output.bai,
    output:
        samtools_stats=temp("../output/rnaseq-wf/samples/{srx}/{srx}.bam.samtools.stats"),
        bamtools_stats=temp("../output/rnaseq-wf/samples/{srx}/{srx}.bam.bamtools.stats"),
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
        samtools=rules.aln_stats.output.samtools_stats,
        bamtools=rules.aln_stats.output.bamtools_stats
    output: "../output/rnaseq-wf/aln_stats/{srx}.parquet"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../../scripts/aln_stats_summary.py"
