import os
import pandas as pd


def main():
    srxs = pd.read_table(snakemake.input.rnaseq_srxs, header=None).squeeze().values

    df = pd.concat(
        [metadata(srxs), hisat2(srxs), samtools_stats(srxs), bamtools_stats(srxs)],
        axis=1,
        sort=True,
    ).rename_axis("srx")

    df.to_csv(snakemake.output[0], sep="\t")


def metadata(srxs):
    return (
        pd.read_table(snakemake.input.srx_metadata, index_col=0)
        .date_created.squeeze()
        .reindex(srxs)
    )


def hisat2(srxs):
    cols = ["num_reads", "num_unaligned", "num_uniquely_aligned", "num_multimappers"]
    return (
        STORE.select("/aln/workflow/hisat2/", "srx == srxs", columns=cols)
        .groupby("srx")
        .sum()
        .assign(pct_uniquely_aligned=lambda x: x.num_uniquely_aligned / x.num_reads * 100)
    )


def samtools_stats(srxs):
    cols = ["reads_properly_paired", "reads_mapped_and_paired"]
    return STORE.select("/aln/workflow/samtools_stats", "index == srxs", columns=cols)


def bamtools_stats(srxs):
    cols = ["Both pairs mapped", "Percent Forward", "Percent Reverse"]
    return STORE.select("/aln/workflow/bamtools_stats", "index == srxs", columns=cols)


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                store="../../output/sra.h5",
                rnaseq_srxs="../../output/paper-wf/rnaseq_srxs.txt",
                srx_metadata="../../output/paper-wf/srx_metadata.tsv",
            )
        )

    try:
        STORE = pd.HDFStore(snakemake.input.store, mode="r")
        main()
    finally:
        STORE.close()
