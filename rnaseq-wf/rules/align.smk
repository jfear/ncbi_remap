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
        strand="../../output/prealn-wf/strand/{srr}.parquet",
        _=rules.atropos_check.output[0]
    output:
        bam=temp("../output/rnaseq-wf/samples/{srx}/{srr}/{srr}.fq.bam"),
        bai=temp("../output/rnaseq-wf/samples/{srx}/{srr}/{srr}.fq.bam.bai")
    log: "../output/rnaseq-wf/samples/{srx}/{srr}/{srr}.fq.bam.log"
    params:
        hisat2_extra='--dta --max-intronlen 300000 --known-splicesite-infile {splice} '.format(splice=config['references']['dmel']['known_splice_sites']),
        samtools_view_extra="--threads 4 -q 20",
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
    output: "../output/rnaseq-wf/hisat2/{srr}.parquet"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../../scripts/hisat2_check.py"

