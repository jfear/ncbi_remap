import os

import gffutils
from gffutils.pybedtools_integration import to_bedtool, featurefuncs


def main():
    db = gffutils.FeatureDB(snakemake.input[0])
    gene = to_bedtool(db.features_of_type('gene')).saveas()
    slopped = gene.slop(b=100, genome='dm6')
    merged = slopped.sort().merge()
    complement = merged.complement(genome='dm6').saveas()

    global cnt
    cnt = 1

    def interName(feature):
        global cnt
        feature = featurefuncs.extend_fields(feature, 4)
        feature.name = 'intergenic{}'.format(cnt)
        cnt += 1
        return feature

    def interGFF(feature):
        gff = featurefuncs.bed2gff(feature)
        gff[1] = 'bedtools'
        gff[2] = 'gene'
        gff.attrs['gene_id'] = gff.name
        return gff

    bed = complement.each(interName).saveas(snakemake.output.bed)
    bed.each(interGFF).saveas(snakemake.output.gtf)


if __name__ == '__main__':
    if os.getenv('SNAKE_DEBUG', False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input="../../lcdb-references/dmel/r6-11/gtf/dmel_r6-11.gtf.db"
        )

    main()
