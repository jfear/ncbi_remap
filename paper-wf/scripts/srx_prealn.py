import os
import pandas as pd


def main():
    df = pd.concat([get_metadata(), aggregate_prealn()], axis=1, sort=False).rename_axis("srx")
    df.to_csv(snakemake.output[0], sep="\t")


def get_metadata():
    cols = ["date_created", "Fear_et_al_library_strategy", "Fear_et_al_library_selection"]
    return pd.read_table(snakemake.input.srx_metadata, index_col="srx").reindex(columns=cols)


def aggregate_prealn():
    aggregator = {
        "layout": "first",
        "strand": "first",
        "libsize_R1": "sum",
        "avgLen_R1": "mean",
        "libsize_R2": "sum",
        "avgLen_R2": "mean",
        "num_unaligned": "sum",
        "num_uniquely_aligned": "sum",
        "num_multimappers": "sum",
        "per_alignment": "mean",
        "reads_properly_paired": "sum",
        "reads_mapped_and_paired": "sum",
        "Both pairs mapped": "sum",
        "Percent Forward": "mean",
        "Percent Reverse": "mean",
        "pct_coding_bases": "mean",
        "pct_intergenic_bases": "mean",
        "pct_intronic_bases": "mean",
        "pct_mrna_bases": "mean",
        "pct_utr_bases": "mean",
        "median_3prime_bias": "mean",
        "median_5prime_bias": "mean",
        "median_5prime_to_3prime_bias": "mean",
        "median_cv_coverage": "mean",
        "first_strand_reads": "sum",
        "second_strand_reads": "sum",
        "pct_first_strand_reads": "mean",
        "pct_second_strand_reads": "mean",
        "pct_duplication": "mean",
        "pct_specific_mapping_adapters": "mean",
        "pct_specific_mapping_dm6": "mean",
        "pct_specific_mapping_ecoli": "mean",
        "pct_specific_mapping_ercc": "mean",
        "pct_specific_mapping_hg19": "mean",
        "pct_specific_mapping_phix": "mean",
        "pct_specific_mapping_rRNA": "mean",
        "pct_specific_mapping_wolbachia": "mean",
        "pct_specific_mapping_yeast": "mean",
    }
    return pd.read_table(snakemake.input.srr, index_col="srr").groupby("srx").agg(aggregator)


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG"):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                srr="../../output/paper-wf/srr_prealn.tsv",
                srx_metadata="../../output/paper-wf/srx_metadata.tsv",
            )
        )

    main()
