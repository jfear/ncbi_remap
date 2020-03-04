"""Aggregates layout into a TSV"""
import os
from pathlib import Path
import csv

from ncbi_remap.aggregation import pull_processed_samples, touch_agg_file
from ncbi_remap.parser import FastqSummary


def main():
    header = ["srr", "libsize_R1", "avgLen_R1", "libsize_R2", "avgLen_R2"]
    touch_agg_file(snakemake.params[0], header)
    complete_srrs = pull_processed_samples(snakemake.params[0])

    with open(snakemake.params[0], "a", newline="") as csvfile:
        writer = csv.writer(csvfile)
        for file_name in snakemake.input:
            srr = Path(file_name).parent.stem
            if srr in complete_srrs:
                continue
            writer.writerow([srr, *FastqSummary(file_name).values])

    Path(snakemake.output[0]).touch()


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG"):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=[
                "../../output/fastq-wf/fastq_info/SRR031714/summary.tsv",
                "../../output/fastq-wf/fastq_info/SRR031715/summary.tsv",
            ],
            params="../../output/agg-prealn-wf/fastq_summary.tsv",
        )

    main()
