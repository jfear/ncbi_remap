import os

import pandas as pd

from lcdblib.utils.chrom_convert import import_conversion
from ncbi_remap.io import GffRow

HEADER = ["FBgn", "gene_symbol", "UCSC_chrom", "FB_chrom", "start", "end", "length", "strand"]


def main():
    df = pd.DataFrame(parse_gtf(), columns=HEADER)
    df.to_csv(snakemake.output[0], sep="\t", index=False)


def parse_gtf():
    fb_mapper = import_conversion("UCSC", "FlyBase")
    with open(snakemake.input[0]) as file_in:
        for row in file_in:
            if row.count("\t") < 5:
                continue

            if row.startswith("#"):
                continue

            gtf_row = GffRow(row)
            if gtf_row.is_gene:
                yield [
                    gtf_row["gene_id"],
                    gtf_row["gene_symbol"],
                    gtf_row.seqid,
                    fb_mapper[gtf_row.seqid],
                    gtf_row.start,
                    gtf_row.end,
                    str(int(gtf_row.end) - int(gtf_row.start)),
                    gtf_row.strand,
                ]



if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG"):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(input="../../lcdb-references/dmel/r6-11/gtf/dmel_r6-11.gtf")

    main()
