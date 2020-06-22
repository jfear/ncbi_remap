import builtins

from ncbi_remap.mock import MockSnake


builtins.snakemake = MockSnake(
    input=dict(
        srx2srr="../../output/srx2srr.csv",
        prealn_done="../../output/prealn-wf/done.txt",
    ),
    threads=8,
    params=dict(
        libsize="../../output/fastq-wf/libsize",
        layout="../../output/fastq-wf/layout",
        fastq_screen="../../output/prealn-wf/fastq_screen",
        atropos="../../output/prealn-wf/atropos",
        hisat2="../../output/prealn-wf/hisat2",
        aln_stats="../../output/prealn-wf/aln_stats",
        rnaseqmetrics="../../output/prealn-wf/rnaseqmetrics",
        genebody_coverage="../../output/prealn-wf/genebody_coverage",
        strand="../../output/prealn-wf/strand",
        markduplicates="../../output/prealn-wf/markduplicates",
        count_summary="../../output/prealn-wf/count_summary",
    ),
)
