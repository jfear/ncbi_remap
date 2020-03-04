from textwrap import dedent
from ncbi_remap import parser


def test_FastqScreen():
    data = dedent(
        """\
        #Fastq_screen version: 0.11.3   #Aligner: bowtie2       #Reads in subset: 100000
        Genome  #Reads_processed        #Unmapped       %Unmapped       #One_hit_one_genome     %One_hit_one_genome     #Multiple_hits_one_genome       %Multiple_hits_one_genome       #One_hit_multiple_genomes       %One_hit_multiple_genomes       Multiple_hits_multiple_genomes  %Multiple_hits_multiple_genomes
        dm6     100517  15153   15.08   73582   73.20   7788    7.75    122     0.12    3872    3.85
        hg19    100517  99683   99.17   10      0.01    14      0.01    109     0.11    701     0.70
        wolbachia       100517  100517  100.00  0       0.00    0       0.00    0       0.00    0       0.00

        %Hit_no_genomes: 15.05"""
    )

    print(data)

