THREADS = 2
GROUP = "bigwig"
RESOURCES = dict(
    mem_gb=lambda wildcards, attempt: attempt * 12,
    time_hr=lambda wildcards, attempt: attempt * 4
)


rule bigwig_first:
    input:
        bedgraph=rules.bamCoverage_first.output[0],
        chromSizes=config['references']['dmel']['chromSizes']
    output: "../output/rnaseq-wf/ucsc_bigwigs/{srx}.first.bw"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    shell:
        'tmpBg=`mktemp --suffix=.bedgraph` '
        '&& bedSort {input.bedgraph} $tmpBg '
        '&& bedGraphToBigWig $tmpBg {input.chromSizes} {output[0]} '
        '&& rm $tmpBg '


rule bigwig_second:
    input:
        bedgraph=rules.bamCoverage_second.output[0],
        chromSizes=config['references']['dmel']['chromSizes']
    output: "../output/rnaseq-wf/uscs_bigwigs/{srx}.second.bw"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    shell:
        'tmpBg=`mktemp --suffix=.bedgraph` '
        '&& bedSort {input.bedgraph} $tmpBg '
        '&& bedGraphToBigWig $tmpBg {input.chromSizes} {output[0]} '
        '&& rm $tmpBg '


rule flybase_bigwig_first:
    input:
        bedgraph=rules.convertToFlybase_first.output[0],
        chromSizes=config['references']['dmel']['fb_chromSizes']
    output: "../output/rnaseq-wf/flybase_bigwigs/{srx}.flybase.first.bw"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    shell:
        'tmpBg=`mktemp --suffix=.bedgraph` '
        '&& bedSort {input.bedgraph} $tmpBg '
        '&& bedGraphToBigWig $tmpBg {input.chromSizes} {output[0]} '
        '&& rm $tmpBg '


rule flybase_bigwig_second:
    input:
        bedgraph=rules.convertToFlybase_second.output[0],
        chromSizes=config['references']['dmel']['fb_chromSizes']
    output: "../output/rnaseq-wf/flybase_bigwigs/{srx}.flybase.second.bw"
    group: GROUP
    threads: THREADS
    resources: **RESOURCES
    shell:
        'tmpBg=`mktemp --suffix=.bedgraph` '
        '&& bedSort {input.bedgraph} $tmpBg '
        '&& bedGraphToBigWig $tmpBg {input.chromSizes} {output[0]} '
        '&& rm $tmpBg '
