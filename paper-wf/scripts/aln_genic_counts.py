import os
import pandas as pd


def main():
    srxs = pd.read_table(snakemake.input[0], header=None).squeeze().values

    pd.concat((read_feature_count(srx) for srx in srxs), axis=1, sort=True).T.to_csv(
        snakemake.output[0], sep="\t"
    )


def read_feature_count(srx):
    try:
        return (
            pd.read_table(snakemake.params.genic_pattern.format(srx=srx), comment="#", index_col=0)
            .iloc[:, -1]
            .squeeze()
            .rename(srx)
            .rename_axis("FBgn")
        )
    except FileNotFoundError:
        return None


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input="../../output/paper-wf/rnaseq_srxs.txt",
            params=dict(genic_pattern="../../output/aln-wf/samples/{srx}/{srx}.bam.counts"),
        )

    main()
