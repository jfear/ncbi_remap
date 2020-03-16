"""Run the FASTQ-Screen Utility"""

rule fastq_screen:
    """Check for contamination."""
    input:
        fastq=rules.fastq_dump.output.r1,
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
    log: temp("../output/prealn-wf/samples/{srx}/{srr}/{srr}_1.fastq_screen.txt.log")
    threads: 8
    resources:
        mem_gb=lambda wildcards, attempt: attempt * 3,
        time_hr=lambda wildcards, attempt: attempt * 1
    script: "../scripts/fastq_screen.py"


rule parse_fastq_screen:
    input: lambda wildcards: queue.expand(rules.fastq_screen.output[0], wildcards.srr)
    output: "../output/prealn-wf/fastq_screen/{srr}.parquet"
    script: "../scripts/parse_fastq_screen.py"