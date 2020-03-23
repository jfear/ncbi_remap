"""Create a feature set for machine learning.

Identify and munge a set of features from the files generated by the prealn-wf
and aln-wf.

Features include:
* CollectRNASeqMetrics
    * PCT_CODING_BASES
    * PCT_UTR_BASES
    * PCT_INTRONIC_BASES
    * PCT_INTERGENIC_BASES
    * PCT_MRNA_BASES
    * MEDIAN_CV_COVERAGE
    * MEDIAN_5PRIME_BIAS
    * MEDIAN_3PRIME_BIAS
* CollectRNASeqMetrics Gene Body Coverage
* Markduplicates
    * PERCENT_DUPLICATION
* Fastq Screen
    * Percent reads mapping to rRNA.
* FeatureCounts 
    * Number of reads mapping to junction
"""
import os
from pathlib import Path

import pandas as pd

# NOTE: features commented out are being dropped because they are repetitive or not important.
FEATURE_AGG = {
    "libsize_R1": "sum",
    "avgLen_R1": "mean",
    "libsize_R2": "sum",
    "avgLen_R2": "mean",
    "layout": "first",
    # "adapters_pct_reads_mapped": "mean",
    "dm6_pct_reads_mapped": "mean",
    # "ecoli_pct_reads_mapped": "mean",
    # "ercc_pct_reads_mapped": "mean",
    # "hg19_pct_reads_mapped": "mean",
    # "phix_pct_reads_mapped": "mean",
    "rRNA_pct_reads_mapped": "mean",
    # "wolbachia_pct_reads_mapped": "mean",
    # "yeast_pct_reads_mapped": "mean",
    # "total_processed": "sum",
    # "total_written": "sum",
    "too_short": "sum",
    # "num_reads": "sum",
    # "num_reads_paired": "sum",
    # "num_reads_unpaired": "sum",
    "num_concordant_reads_unaligned": "sum",
    "num_concordant_reads_uniquely_aligned": "sum",
    "num_concordant_multimappers": "sum",
    "num_discordant_reads_aligned": "sum",
    "num_unaligned": "sum",
    "num_uniquely_aligned": "sum",
    "num_multimappers": "sum",
    "per_alignment": "mean",
    "reads_MQ0": "sum",
    "average_quality": "mean",
    "insert_size_average": "mean",
    "insert_size_standard_deviation": "mean",
    "inward_oriented_pairs": "sum",
    "outward_oriented_pairs": "sum",
    "pairs_with_other_orientation": "sum",
    "pairs_on_different_chromosomes": "sum",
    "Percent Forward": "mean",
    "Percent Reverse": "mean",
    "percent_coding_bases": "mean",
    "percent_utr_bases": "mean",
    "percent_intronic_bases": "mean",
    "percent_intergenic_bases": "mean",
    "percent_mrna_bases": "mean",
    "median_cv_coverage": "mean",
    # "median_5prime_bias": "mean",
    # "median_3prime_bias": "mean",
    "median_5prime_to_3prime_bias": "mean",
    "pos_0": "mean",
    "pos_1": "mean",
    "pos_2": "mean",
    "pos_3": "mean",
    "pos_4": "mean",
    "pos_5": "mean",
    "pos_6": "mean",
    "pos_7": "mean",
    "pos_8": "mean",
    "pos_9": "mean",
    "pos_10": "mean",
    "pos_11": "mean",
    "pos_12": "mean",
    "pos_13": "mean",
    "pos_14": "mean",
    "pos_15": "mean",
    "pos_16": "mean",
    "pos_17": "mean",
    "pos_18": "mean",
    "pos_19": "mean",
    "pos_20": "mean",
    "pos_21": "mean",
    "pos_22": "mean",
    "pos_23": "mean",
    "pos_24": "mean",
    "pos_25": "mean",
    "pos_26": "mean",
    "pos_27": "mean",
    "pos_28": "mean",
    "pos_29": "mean",
    "pos_30": "mean",
    "pos_31": "mean",
    "pos_32": "mean",
    "pos_33": "mean",
    "pos_34": "mean",
    "pos_35": "mean",
    "pos_36": "mean",
    "pos_37": "mean",
    "pos_38": "mean",
    "pos_39": "mean",
    "pos_40": "mean",
    "pos_41": "mean",
    "pos_42": "mean",
    "pos_43": "mean",
    "pos_44": "mean",
    "pos_45": "mean",
    "pos_46": "mean",
    "pos_47": "mean",
    "pos_48": "mean",
    "pos_49": "mean",
    "pos_50": "mean",
    "pos_51": "mean",
    "pos_52": "mean",
    "pos_53": "mean",
    "pos_54": "mean",
    "pos_55": "mean",
    "pos_56": "mean",
    "pos_57": "mean",
    "pos_58": "mean",
    "pos_59": "mean",
    "pos_60": "mean",
    "pos_61": "mean",
    "pos_62": "mean",
    "pos_63": "mean",
    "pos_64": "mean",
    "pos_65": "mean",
    "pos_66": "mean",
    "pos_67": "mean",
    "pos_68": "mean",
    "pos_69": "mean",
    "pos_70": "mean",
    "pos_71": "mean",
    "pos_72": "mean",
    "pos_73": "mean",
    "pos_74": "mean",
    "pos_75": "mean",
    "pos_76": "mean",
    "pos_77": "mean",
    "pos_78": "mean",
    "pos_79": "mean",
    "pos_80": "mean",
    "pos_81": "mean",
    "pos_82": "mean",
    "pos_83": "mean",
    "pos_84": "mean",
    "pos_85": "mean",
    "pos_86": "mean",
    "pos_87": "mean",
    "pos_88": "mean",
    "pos_89": "mean",
    "pos_90": "mean",
    "pos_91": "mean",
    "pos_92": "mean",
    "pos_93": "mean",
    "pos_94": "mean",
    "pos_95": "mean",
    "pos_96": "mean",
    "pos_97": "mean",
    "pos_98": "mean",
    "pos_99": "mean",
    "pos_100": "mean",
    "strand": "first",
    "unpaired_reads_examined": "sum",
    "read_pairs_examined": "sum",
    "unpaired_read_duplicates": "sum",
    "read_pair_duplicates": "sum",
    "percent_duplication": "mean",
    "estimated_library_size": "sum",
    "number_genic_reads": "sum",
    "percent_genes_on": "mean",
    "number_junction_reads": "sum",
    "number_junctions_on": "sum",
}


def main():
    srx2srr = pd.read_csv(snakemake.input[0], index_col="srr")
    done = {pth.stem for pth in Path(snakemake.params.done).iterdir()}
    df = pd.concat(
        [
            pd.read_parquet(v, engine="pyarrow", use_threads=True)
            for k, v in snakemake.params.items()
            if k != "done"
        ],
        axis=1,
        sort=False,
    ).reindex(done)

    df.join(srx2srr).groupby("srx").agg(FEATURE_AGG).to_parquet(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        import sys

        sys.path.insert(0, "../../src")
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input="../../output/srx2srr_full.csv",
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
                done="../../output/prealn-wf/done",
            ),
        )

    main()
