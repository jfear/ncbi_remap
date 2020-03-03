"""Aggregates layout into a TSV"""
import os
from pathlib import Path
import csv
from ncbi_remap.snakemake import get_flag

HEADER = ["srx", "srr", "layout"]


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
            writer.writerow([queue.get_srx(srr), srr, get_flag(file_name)])

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


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG"):
        from ncbi_remap.debug import snakemake_debug
        from ncbi_remap.queue import Queue

        queue = Queue(targets="../../output/prealn-wf/done", srx2srr="../../output/srx2srr.csv")

        snakemake = snakemake_debug(
            input=[
                "../../output/fastq-wf/fastq_info/SRR031714/LAYOUT",
                "../../output/fastq-wf/fastq_info/SRR031715/LAYOUT",
            ],
            params=dict(queue=queue, data_store="../../output/agg-prealn-wf/layout.tsv"),
        )

    main()
