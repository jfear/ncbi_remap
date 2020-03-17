"""Aggregate junction and genic counts."""
import pandas as pd

# Gene Counts
gene_counts = pd.read_table(
    snakemake.input.counts, comment="#", index_col=0
).iloc[:, -1]
genic_reads = gene_counts.sum()
percent_genes_on = (gene_counts > 0).mean() * 100

# Junction Counts
junction_counts = pd.read_table(
    snakemake.input.jcounts, comment="#", index_col=0
).iloc[:, -1]
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
    index=[snakemake.wildcards.srr],
)

df.to_parquet(snakemake.output[0])
