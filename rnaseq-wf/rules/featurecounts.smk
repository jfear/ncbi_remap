"""Run Feature Counts"""
THREADS = 4
GROUP = "feature_counts"
RESOURCES = dict(
    mem_gb=lambda wildcards, attempt: attempt * 4,
    time_hr=lambda wildcards, attempt: attempt * 4
)

##################################################################################
# Gene Counts
##################################################################################
rule feature_counts:
    input:
        annotation=config['references']['dmel']['gtf'],
        bam=rules.bamMerge.output.bam,
        layout=lambda wildcards: queue.expand(rules.fastq_check.output.layout, wildcards.srx)[0],
        strand=lambda wildcards: queue.expand("../output/prealn-wf/strand/{srr}.parquet", wildcards.srx)[0]
    output:
        counts=temp("../output/rnaseq-wf/samples/{srx}/{srx}.bam.counts"),
        jcounts=temp("../output/rnaseq-wf/samples/{srx}/{srx}.bam.counts.jcounts"),
        summary=temp("../output/rnaseq-wf/samples/{srx}/{srx}.bam.counts.summary"),
    log: "../output/rnaseq-wf/samples/{srx}/{srx}.bam.counts.log"
    params:
        extra_pe = '-p -P -C -J ',
        extra_se = '-J '
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../../scripts/featurecounts.py"


rule parse_gene_count:
    input: rules.feature_counts.output.counts
    output: "../output/rnaseq-wf/gene_counts/{srx}.parquet"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../scripts/parse_feature_counts.py"


rule parse_junction_count:
    input: rules.feature_counts.output.jcounts
    output: "../output/rnaseq-wf/junction_counts/{srx}.parquet"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../scripts/parse_junction_counts.py"

##################################################################################
# Intergenic Counts
##################################################################################
rule feature_counts_intergenic:
    input:
        annotation=config['references']['dmel']['intergenic'],
        bam=rules.bamMerge.output.bam,
        layout=lambda wildcards: queue.expand(rules.fastq_check.output.layout, wildcards.srx)[0],
        strand=lambda wildcards: queue.expand("../output/prealn-wf/strand/{srr}.parquet", wildcards.srx)[0]
    output:
        counts=temp("../output/rnaseq-wf/samples/{srx}/{srx}.bam.intergenic.counts"),
        jcounts=temp("../output/rnaseq-wf/samples/{srx}/{srx}.bam.intergenic.counts.jcounts"),
        summary=temp("../output/rnaseq-wf/samples/{srx}/{srx}.bam.intergenic.counts.summary"),
    log: "../output/rnaseq-wf/samples/{srx}/{srx}.bam.intergenic.counts.log"
    params:
        extra_pe = '-p -P -C -J -t gene ',
        extra_se = '-J -t gene '
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../../scripts/featurecounts.py"


rule parse_intergenic_count:
    input: rules.feature_counts_intergenic.output.counts
    output: "../output/rnaseq-wf/intergenic_counts/{srx}.parquet"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../scripts/parse_feature_counts.py"


##################################################################################
# Segment Counts
##################################################################################
rule feature_counts_segments:
    input:
        annotation=config['references']['dmel']['nonstranded_segments'],
        bam=rules.bamMerge.output.bam,
        layout=lambda wildcards: queue.expand(rules.fastq_check.output.layout, wildcards.srx)[0],
        strand=lambda wildcards: queue.expand("../output/prealn-wf/strand/{srr}.parquet", wildcards.srx)[0]
    output:
        counts=temp("../output/rnaseq-wf/samples/{srx}/{srx}.bam.exon_segments.counts"),
        jcounts=temp("../output/rnaseq-wf/samples/{srx}/{srx}.bam.exon_segments.counts.jcounts"),
        summary=temp("../output/rnaseq-wf/samples/{srx}/{srx}.bam.exon_segments.counts.summary"),
    log: "../output/rnaseq-wf/samples/{srx}/{srx}.bam.exon_segments.counts.log"
    params:
        extra_pe = '-p -P -C -J -t segment -g ID ',
        extra_se = '-J -t segment -g ID '
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../../scripts/featurecounts.py"


rule parse_segment_count:
    input: rules.feature_counts_segments.output.counts
    output: "../output/rnaseq-wf/segment_counts/{srx}.parquet"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../scripts/parse_feature_counts.py"


##################################################################################
# Fusion Counts
##################################################################################
rule feature_counts_fusions:
    input:
        annotation=config['references']['dmel']['nonstranded_fusions'],
        bam=rules.bamMerge.output.bam,
        layout=lambda wildcards: queue.expand(rules.fastq_check.output.layout, wildcards.srx)[0],
        strand=lambda wildcards: queue.expand("../output/prealn-wf/strand/{srr}.parquet", wildcards.srx)[0]
    output:
        counts=temp("../output/rnaseq-wf/samples/{srx}/{srx}.bam.exon_fusions.counts"),
        jcounts=temp("../output/rnaseq-wf/samples/{srx}/{srx}.bam.exon_fusions.counts.jcounts"),
        summary=temp("../output/rnaseq-wf/samples/{srx}/{srx}.bam.exon_fusions.counts.summary"),
    log: "../output/rnaseq-wf/samples/{srx}/{srx}.bam.exon_fusions.counts.log"
   params:
        extra_pe = '-p -P -C -J -t fusion -g ID ',
        extra_se = '-J -t fusion -g ID '
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../../scripts/featurecounts.py"


rule parse_fusion_count:
    input: rules.feature_counts_fusions.output.counts
    output: "../output/rnaseq-wf/fusion_counts/{srx}.parquet"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    script: "../scripts/parse_feature_counts.py"
