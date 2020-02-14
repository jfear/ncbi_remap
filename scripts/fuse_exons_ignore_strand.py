"""Builds segment GTF

This script segments exons based on overlap and strand. It does the
following:

1. Sorts exons by chrom, strand, start, end
2. Checks for overlapping exons
3. Segments those exons
4. Saves a GTF formatted file.

"""
import os

import gffutils

from ncbi_remap.gtf import FeatureAccumulator


def main():
    db = gffutils.FeatureDB(snakemake.input[0])
    exons = db.features_of_type("exon", order_by=["seqid", "start", "end"])
    with open(snakemake.output[0], "w") as file_out:
        accumulartor = FeatureAccumulator("fusion", stranded=False)
        for exon in exons:
            if accumulartor.add_feature(exon):
                continue

            file_out.write(accumulartor.merge())
            accumulartor.reset(exon)

        file_out.write(accumulartor.merge())


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="rnaseq-analysis-wf",
            input="../lcdb-references/dmel/r6-11/gtf/dmel_r6-11.gtf.db",
        )

    main()
