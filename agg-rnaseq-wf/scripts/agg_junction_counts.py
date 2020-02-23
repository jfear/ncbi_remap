import os
from pathlib import Path

from ncbi_remap.parser import parse_featureCounts_jcounts
from ncbi_remap.aggregation import RnaSeqAggregator


def main():
    RnaSeqAggregator(
        count_files=snakemake.input,
        data_store=snakemake.params[0],
        parser=parse_featureCounts_jcounts,
        sample_type="srx",
        threads=snakemake.threads,
    )

    Path(snakemake.output[0]).touch()


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG"):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=[
                "../../output/rnaseq-wf/samples/SRX014459/SRX014459.bam.counts.jcounts",
                "../../output/rnaseq-wf/samples/SRX014472/SRX014472.bam.counts.jcounts",
            ],
            params="../../output/agg-rnaseq-wf/junction_counts.tsv",
            threads=2,
        )

    main()
