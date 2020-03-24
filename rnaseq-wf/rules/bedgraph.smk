THREADS = 8
GROUP = "bedgraph"
RESOURCES = dict(
    mem_gb=lambda wildcards, attempt: attempt * 16,
    time_hr=lambda wildcards, attempt: attempt * 8
)

bamCoverage_options = (
    '--outFileFormat bedgraph '
    '--binSize 1 '
    '--effectiveGenomeSize 129000000 '
    '--normalizeUsing RPGC '
    '--ignoreForNormalization chrX '
)

rule bamCoverage_first:
    input:
        bam=rules.expMerge.output.bam,
        bai=rules.expMerge.output.bai,
    output: temp('../output/rnaseq-wf/samples/{srx}/{srx}.first.bedgraph')
    params:
        extra=bamCoverage_options + '--filterRNAstrand forward'
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../scripts/bamcoverage.py"


rule bamCoverage_second:
    input:
        bam=rules.expMerge.output.bam,
        bai=rules.expMerge.output.bai,
    output: temp('../output/rnaseq-wf/samples/{srx}/{srx}.second.bedgraph')
    params:
        extra=bamCoverage_options + '--filterRNAstrand reverse'
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../scripts/bamcoverage.py"


rule convertToFlybase_first:
    input: rules.bamCoverage_first.output[0]
    output: temp("../output/rnaseq-wf/samples/{srx}/{srx}.flybase.first.bedgraph")
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../scripts/chrom_convert.py"


rule convertToFlybase_second:
    input: rules.bamCoverage_second.output[0]
    output: temp("../output/rnaseq-wf/samples/{srx}/{srx}.flybase.second.bedgraph")
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../scripts/chrom_convert.py"
