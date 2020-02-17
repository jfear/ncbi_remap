"""Aggregate counts tables from the rnaseq-wf"""
import os
import multiprocessing
from multiprocessing import Pool
from pathlib import Path
import csv

import pandas as pd
from more_itertools import grouper

CPUS = int(os.getenv("SLURM_CPUS_PER_TASK", max(1, multiprocessing.cpu_count() - 2)))

INPUT_DIR = "../output/rnaseq-wf"
GENE_OUTPUT = "../output/agg-rnaseq-wf/gene_counts.tsv"
INTERGENIC_OUTPUT = "../output/agg-rnaseq-wf/intergenic_counts.tsv"
JUNCTION_OUTPUT = "../output/agg-rnaseq-wf/junction_counts.tsv"
SEGMENT_OUTPUT = "../output/agg-rnaseq-wf/segment_counts.tsv"
FUSION_OUTPUT = "../output/agg-rnaseq-wf/fusion_counts.tsv"


def main():
    pool = multiprocessing.Pool(CPUS)
    srxs = get_completed_srxs()
    aggregate_gene_counts(srxs, pool)
    aggregate_intergenic_counts(srxs, pool)
    aggregate_junction_counts(srxs, pool)
    aggregate_segment_counts(srxs, pool)
    aggregate_fusion_counts(srxs, pool)


def get_completed_srxs() -> set:
    return {pth.stem for pth in Path(INPUT_DIR, "done").iterdir()}


def aggregate_gene_counts(srxs, pool):
    parsed_srxs = check_already_parsed(GENE_OUTPUT)
    files = [Path(INPUT_DIR, "samples", srx, f"{srx}.bam.counts") for srx in srxs - parsed_srxs]
    for _files in grouper(files, 1000):
        df = load_data(_files, pool)
        write_output(GENE_OUTPUT, df)


def aggregate_intergenic_counts(srxs, pool):
    parsed_srxs = check_already_parsed(INTERGENIC_OUTPUT)
    files = [
        Path(INPUT_DIR, "samples", srx, f"{srx}.bam.intergenic.counts")
        for srx in srxs - parsed_srxs
    ]
    for _files in grouper(files, 1000):
        df = load_data(_files, pool)
        write_output(INTERGENIC_OUTPUT, df)


def aggregate_junction_counts(srxs, pool):
    parsed_srxs = check_already_parsed(JUNCTION_OUTPUT)
    files = [
        Path(INPUT_DIR, "samples", srx, f"{srx}.bam.counts.jcounts") for srx in srxs - parsed_srxs
    ]
    for _files in grouper(files, 1000):
        df = load_data(_files, pool, parse_junction_counts)
        write_output(JUNCTION_OUTPUT, df)


def aggregate_segment_counts(srxs, pool):
    parsed_srxs = check_already_parsed(SEGMENT_OUTPUT)
    files = [
        Path(INPUT_DIR, "samples", srx, f"{srx}.bam.exon_segments.counts")
        for srx in srxs - parsed_srxs
    ]
    for _files in grouper(files, 1000):
        df = load_data(_files, pool)
        write_output(SEGMENT_OUTPUT, df)


def aggregate_fusion_counts(srxs, pool):
    parsed_srxs = check_already_parsed(FUSION_OUTPUT)
    files = [
        Path(INPUT_DIR, "samples", srx, f"{srx}.bam.exon_fusions.counts")
        for srx in srxs - parsed_srxs
    ]
    for _files in grouper(files, 1000):
        df = load_data(_files, pool)
        write_output(FUSION_OUTPUT, df)


def check_already_parsed(file_name) -> set:
    if Path(file_name).exists():
        reader = csv.reader(open(file_name, "r"), delimiter="\t")
        return {row[0] for row in reader}
    return set()


def parse_feature_counts(file_name) -> pd.Series:
    if file_name is None:
        return

    srx = Path(file_name).parent.stem
    return (
        pd.read_csv(file_name, sep="\t", comment="#", index_col=0)
        .iloc[:, -1]
        .squeeze()
        .rename(srx)
        .to_frame()
        .T.rename_axis("srx")
    )


def parse_junction_counts(file_name) -> pd.Series:
    if file_name is None:
        return

    srx = Path(file_name).parent.stem
    df = pd.read_csv(file_name, sep="\t", comment="#").assign(srx=srx).set_index("srx")
    cols = df.columns.tolist()
    cols[-1] = "Count"
    df.columns = cols
    return df


def load_data(files, pool, func=parse_feature_counts) -> pd.DataFrame:
    try:
        return pd.concat(pool.map(func, files))
    except ValueError:
        pass


def write_output(file_name, df):
    if df is None:
        return

    if Path(file_name).exists():
        cols = pd.read_csv(file_name, index_col=0, nrows=0).columns
        df.reindex(columns=cols).to_csv(file_name, header=False, mode="a", sep="\t")
    else:
        df.to_csv(file_name, sep="\t")


if __name__ == "__main__":
    Path(GENE_OUTPUT).parent.mkdir(exist_ok=True)
    main()
