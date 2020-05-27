"""Run the FASTQ-Screen Utility"""
THREADS = 8
GROUP = "fastq_screen"
RESOURCES = dict(
    mem_gb=lambda wildcards, attempt: attempt * 4,
    time_hr=lambda wildcards, attempt: attempt * 2
)

rule fastq_screen:
    """Check for contamination."""
    input:
        fastq=rules.fastq_dump.output.r1,
        layout=rules.fastq_dump.output.layout,
        dm6=config['references']['dmel']['bowtie2'],
        hg19=config['references']['human']['bowtie2'],
        wolbachia=config['references']['wolbachia']['bowtie2'],
        ecoli=config['references']['ecoli']['bowtie2'],
        yeast=config['references']['yeast']['bowtie2'],
        rRNA=config['references']['rRNA']['bowtie2'],
        phix=config['references']['phix']['bowtie2'],
        ercc=config['references']['ercc']['bowtie2'],
        adapters=config['references']['adapters']['bowtie2']
    output: txt=temp("../output/prealn-wf/samples/{srx}/{srr}/{srr}_1.fastq_screen.txt")
    log: "../output/prealn-wf/samples/{srx}/{srr}/{srr}_1.fastq_screen.txt.log"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../scripts/fastq_screen.py"


rule fastq_screen_summary:
    input: lambda wildcards: queue.expand(rules.fastq_screen.output[0], wildcards.srr)
    output: "../output/prealn-wf/fastq_screen/{srr}.parquet"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../scripts/fastq_screen_summary.py"