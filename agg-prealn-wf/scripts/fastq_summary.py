"""Aggregates layout into a TSV"""
import os
from pathlib import Path
import csv

HEADER = ["srx", "srr", "libsize_R1", "avgLen_R1", "libsize_R2", "avgLen_R2"]


def main():
    queue = snakemake.params.queue
    data_store = snakemake.params.data_store

    complete_srrs = check_srrs(data_store)

    with open(data_store, "a") as file_out:
        writer = csv.writer(file_out, delimiter="\t")
        for file_name in snakemake.input:
            srr = Path(file_name).parent.stem
            if srr in complete_srrs:
                continue
            writer.writerow([queue.get_srx(srr), srr, *parse_summary(file_name)])

    Path(snakemake.output[0]).touch()


def check_srrs(file_name: str) -> set:
    if Path(file_name).exists():
        with open(file_name, "r") as file_in:
            reader = csv.reader(file_in, delimiter="\t")
            _header = next(reader)
            return {row[1] for row in reader}

    with open(file_name, "w") as file_out:
        writer = csv.writer(file_out, delimiter="\t")
        writer.writerow(HEADER)
    return set()


def parse_summary(file_name: str) -> list:
    with open(file_name, "r") as file_in:
        reader = csv.reader(file_in, delimiter="\t")
        _header = next(reader)
        row = next(reader)
        ncols = len(row)
        if ncols == 4:
            return row
        elif ncols == 6:
            # Older version of summary table. Contains md5sums that I am now
            # ignoring.
            return [row[1], row[2], row[4], row[5]]


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG"):
        from ncbi_remap.debug import snakemake_debug
        from ncbi_remap.queue import Queue

        queue = Queue(targets="../../output/prealn-wf/done", srx2srr="../../output/srx2srr.csv")

        snakemake = snakemake_debug(
            input=[
                "../../output/fastq-wf/fastq_info/SRR031714/summary.tsv",
                "../../output/fastq-wf/fastq_info/SRR031715/summary.tsv",
            ],
            params=dict(queue=queue, data_store="../../output/agg-prealn-wf/fastq_summary.tsv"),
        )

    main()
