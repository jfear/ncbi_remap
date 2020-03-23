"""Aggregate junction and genic counts."""
import numpy as np
import pandas as pd


def main():
    gene_counts = get_counts(snakemake.input.counts)
    genic_reads = gene_counts.sum()
    percent_genes_on = (gene_counts > 0).mean() * 100

    junction_counts = get_counts(snakemake.input.jcounts)
    junction_reads = junction_counts.sum()
    number_junctions_on = junction_counts.shape[0]

    df = pd.DataFrame(
        [[genic_reads, percent_genes_on, junction_reads, number_junctions_on]],
        columns=[
            "number_genic_reads",
            "percent_genes_on",
            "number_junction_reads",
            "number_junctions_on",
        ],
        index=pd.Index([snakemake.wildcards.srr], name="srr"),
    )

    df.to_parquet(snakemake.output[0])


def get_counts(file_name):
    col_name = pd.read_table(file_name, comment="#", nrows=1).columns[-1]
    return pd.read_table(file_name, comment="#", usecols=[col_name], dtype=np.int64).values


if __name__ == "__main__":
    main()
