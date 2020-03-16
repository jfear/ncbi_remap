"""Summarizes fastq screen results
Cacluates the number of reads mapping to each specific reference. Ignores reads the map to multiple references.

Goes from this:
| reference   |   multiple_hits_multiple_libraries_count |   multiple_hits_multiple_libraries_percent |   multiple_hits_one_library_count |   multiple_hits_one_library_percent |   one_hit_multiple_libraries_count |   one_hit_multiple_libraries_percent |   one_hit_one_library_count |   one_hit_one_library_percent |   reads_processed_count |   unmapped_count |   unmapped_percent |
|:------------|-----------------------------------------:|-------------------------------------------:|----------------------------------:|------------------------------------:|-----------------------------------:|-------------------------------------:|----------------------------:|------------------------------:|------------------------:|-----------------:|-------------------:|
| adapters    |                                       48 |                                       0.05 |                                 0 |                                0    |                                  0 |                                 0    |                           0 |                          0    |                   99973 |            99925 |              99.95 |
| dm6         |                                     1713 |                                       1.71 |                              6278 |                                6.28 |                                224 |                                 0.22 |                       88393 |                         88.42 |                   99973 |             3365 |               3.37 |
| ecoli       |                                        1 |                                       0    |                                 0 |                                0    |                                  0 |                                 0    |                           2 |                          0    |                   99973 |            99970 |             100    |
...

To this:
|            |   adapters_pct_reads_mapped |   dm6_pct_reads_mapped |   ecoli_pct_reads_mapped |   ercc_pct_reads_mapped |   hg19_pct_reads_mapped |   phix_pct_reads_mapped |   rRNA_pct_reads_mapped |   wolbachia_pct_reads_mapped |   yeast_pct_reads_mapped |
|:-----------|----------------------------:|-----------------------:|-------------------------:|------------------------:|------------------------:|------------------------:|------------------------:|-----------------------------:|-------------------------:|
| SRR0000001 |                           0 |                94.6966 |               0.00200054 |                       0 |               0.0160043 |                       0 |              0.00100027 |                            0 |               0.00500135 |

"""
import os
import sys

import pandas as pd

sys.path.insert(0, "../src")
from ncbi_remap.parser import parse_fastq_screen


df = parse_fastq_screen(snakemake.input[0]).set_index("reference")
summarized = (
    (
        (df.one_hit_one_library_count + df.multiple_hits_one_library_count)
        / df.reads_processed_count
        * 100
    )
    .rename(snakemake.wildcards.srr)
    .to_frame()
    .T
)

summarized.columns = [f"{col}_pct_reads_mapped" for col in summarized.columns]
summarized.to_parquet(snakemake.output[0])
