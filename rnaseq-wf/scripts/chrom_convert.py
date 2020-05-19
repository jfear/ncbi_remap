import pandas as pd
import pybedtools
import ssl

# I was getting an SSL error.
# From: https://stackoverflow.com/a/56230607/4605992
ssl._create_default_https_context = ssl._create_unverified_context

def main():
    conversion_table = import_conversion()
    pybedtools_convert(snakemake.input[0], snakemake.output[0], conversion_table)


def import_conversion():
    """Download conversion table from NCBI.

        col#   Header
        0      Sequence-Name
        1      Sequence-Role
        2      Assigned-Molecule
        3      Assigned-Molecule-Location/Type
        4      GenBank-Accn
        5      Relationship
        6      RefSeq-Accn
        7      Assembly-Unit
        8      Sequence-Length
        9      UCSC-style-name
    """
    url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_assembly_report.txt"

    # This table has a small error were FlyBase mitochondrion_genome is MT
    df = pd.read_table(url, comment="#", header=None).replace("MT", "mitochondrion_genome")
    mapping = {"FlyBase": 0, "UCSC": 9, "GenBank": 4, "RefSeq": 6}
    return {k: v for k, v in df[[mapping["UCSC"], mapping["FlyBase"]]].values}


def convertFeature(f, mapper):
    f.chrom = mapper[f.chrom]
    return f


def pybedtools_convert(input, output, mapper):
    """ Use pybedtools to convert chromosomes in BED, GTF, or GFF. """
    bt = pybedtools.BedTool(input)
    bt.each(convertFeature, mapper).saveas(output)


if __name__ == "__main__":
    main()
