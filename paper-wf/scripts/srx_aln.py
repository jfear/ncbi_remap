import os
import pandas as pd


def main():
    srxs = (
        pd.read_parquet(snakemake.input.library_strategy)
        .query("Fear_et_al_library_strategy == 'RNA-Seq'")
        .index
    )

    df = pd.concat(
        [hisat2(srxs), samtools_stats(srxs), bamtools_stats(srxs), agg_feature_counts(srxs)],
        axis=1,
        sort=True,
    )

    df.to_csv(snakemake.output[0], sep="\t")


def hisat2(srxs):
    cols = ["num_reads", "num_unaligned", "num_uniquely_aligned", "num_multimappers"]
    return (
        STORE.select("/aln/workflow/hisat2/", "srx == srxs", columns=cols)
        .groupby("srx")
        .sum()
        .assign(pct_uniquely_alignmed=lambda x: x.num_uniquely_aligned / x.num_reads * 100)
    )


def samtools_stats(srxs):
    cols = ["reads_properly_paired", "reads_mapped_and_paired"]
    return STORE.select("/aln/workflow/samtools_stats", "index == srxs", columns=cols)


def bamtools_stats(srxs):
    cols = ["Both pairs mapped", "Percent Forward", "Percent Reverse"]
    return STORE.select("/aln/workflow/bamtools_stats", "index == srxs", columns=cols)


def agg_feature_counts(srxs):
    return pd.concat((read_feature_count(srx) for srx in srxs), axis=1, sort=True).T


def read_feature_count(srx):
    return (
        pd.read_table(
            snakemake.params.genic_pattern.format(srx=srx), comment="#", index_col=0
        )
        .iloc[:, -1]
        .squeeze()
        .rename(srx)
        .rename_axis("FBgn")
    )


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                store="../../output/sra.h5",
                library_strategy="../../output/library_strategy-wf/summarized_metadata.parquet",
            ),
            params=dict(genic_pattern="../../output/aln-wf/samples/{srx}/{srx}.bam.counts"),
        )

    try:
        STORE = pd.HDFStore(snakemake.input.store, mode="r")
        main()
    finally:
        STORE.close()
