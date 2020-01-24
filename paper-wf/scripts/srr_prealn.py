import os

import pandas as pd


def main():
    df = pd.concat(
        [
            DATA_STORE.select("layout", "srr == SRRS"),
            DATA_STORE.select("strand", "srr == SRRS"),
            fastq(),
            hisat2(),
            samtools_stats(),
            bamtools_stats(),
            collectrnaseqmetrics(),
            markduplicates(),
            fastq_screen(),
        ],
        axis=1,
        sort=True,
    )

    df.to_csv(snakemake.output[0], sep="\t")


def fastq():
    return DATA_STORE.select("/prealn/workflow/fastq", "srr == SRRS")


def hisat2():
    cols = ["num_unaligned", "num_uniquely_aligned", "num_multimappers", "per_alignment"]
    return DATA_STORE.select("/prealn/workflow/hisat2/", "srr == SRRS", columns=cols)


def samtools_stats():
    cols = ["reads_properly_paired", "reads_mapped_and_paired"]
    return DATA_STORE.select("/prealn/workflow/samtools_stats", "srr == SRRS", columns=cols)


def bamtools_stats():
    cols = ["Both pairs mapped", "Percent Forward", "Percent Reverse"]
    return DATA_STORE.select("/prealn/workflow/bamtools_stats", "srr == SRRS", columns=cols)


def collectrnaseqmetrics():
    cols = [
        "PCT_CODING_BASES",
        "PCT_INTERGENIC_BASES",
        "PCT_INTRONIC_BASES",
        "PCT_MRNA_BASES",
        "PCT_UTR_BASES",
        "MEDIAN_3PRIME_BIAS",
        "MEDIAN_5PRIME_BIAS",
        "MEDIAN_5PRIME_TO_3PRIME_BIAS",
        "MEDIAN_CV_COVERAGE",
        "CORRECT_STRAND_READS",
        "INCORRECT_STRAND_READS",
        "PCT_CORRECT_STRAND_READS",
    ]

    mapper = {
        "CORRECT_STRAND_READS": "first_strand_reads",
        "INCORRECT_STRAND_READS": "second_strand_reads",
        "PCT_CORRECT_STRAND_READS": "pct_first_strand_reads",
    }

    df = DATA_STORE.select(
        "/prealn/workflow/collectrnaseqmetrics/first", "srr == SRRS", columns=cols
    ).rename(columns=lambda x: mapper.get(x, x.lower()))

    df["pct_second_strand_reads"] = df["second_strand_reads"] / (
        df["first_strand_reads"] + df["second_strand_reads"]
    )

    for column in df.columns:
        if "pct_" in column:
            df[column] = df[column] * 100

    return df


def markduplicates():
    return (
        DATA_STORE.select(
            "/prealn/workflow/markduplicates", "srr == SRRS", columns=["PERCENT_DUPLICATION"]
        )
        .squeeze()
        .rename("pct_duplication")
        * 100
    )


def fastq_screen():
    cols = [
        "reference",
        "reads_processed_count",
        "one_hit_one_library_count",
        "multiple_hits_one_library_count",
    ]
    df = DATA_STORE.select("/prealn/workflow/fastq_screen", "srr == SRRS", columns=cols).set_index(
        "reference", append=True
    )
    df["pct_specific_mapping"] = (
        (df["one_hit_one_library_count"] + df["multiple_hits_one_library_count"])
        / df["reads_processed_count"]
        * 100
    )

    return df.pct_specific_mapping.unstack().rename(columns=lambda x: f"pct_specific_mapping_{x}")


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG"):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(input="../../output/sra.h5")

    DATA_STORE = pd.HDFStore(snakemake.input[0], mode="r")
    SRRS = DATA_STORE["prealn/complete"].srr.unique()
    main()
    DATA_STORE.close()
