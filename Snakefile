"""Get things started."""
workdir: '.'
configfile: 'config/reference_config.yaml'

rule targets:
    input:
        # Database Info
        "output/srx2srr.csv",
        "output/rnaseq_srxs.txt",
        # Additional References
        "output/known_splice_sites_r6-11.txt",
        "output/dmel_r6-11.intergenic.gtf",
        "output/nonstranded_exon_segments.gtf",
        "output/stranded_exon_segments.gtf",
        "output/nonstranded_exon_fusions.gtf",
        "output/stranded_exon_fusions.gtf",


rule sra2mongo:
    output: "output/db_download.date"
    shell: """
    sra2mongo --api $ENTREZ_API_KEY --query '"Drosophila melanogaster"[orgn]' && \
    date > {output[0]}
    """

rule srx2srr:
    input: rules.sra2mongo.output[0]
    output: "output/srx2srr.csv",
    script: "scripts/srx2srr.py"


rule rnaseq_srxs:
    """Pull out a list of SRXs annotated as RNA-Seq in SRA."""
    input: rules.sra2mongo.output[0]
    output: "output/rnaseq_srxs.txt"
    script: "scripts/rnaseq_srxs.py"


###############################################################################
# Build References
###############################################################################
rule hisat2_splice_site:
    """Generate splicesite information from known annotations."""
    input: gtf=config['references']['dmel']['gtf'].lstrip("../")
    output: "output/known_splice_sites_r6-11.txt"
    conda: "./prealn-wf/wrappers/hisat2/align/environment.yaml"
    shell: "hisat2_extract_splice_sites.py {input.gtf} > {output}"


rule intergenic:
    input: config['references']['dmel']['db'].lstrip("../")
    output:
        bed="output/dmel_r6-11.intergenic.bed",
        gtf="output/dmel_r6-11.intergenic.gtf"
    script: "scripts/intergenic_region.py"


rule segment_exons_ignore_strand:
    input: config['references']['dmel']['db'].lstrip("../")
    output: "output/nonstranded_exon_segments.gtf"
    script: "scripts/segment_exons_ignore_strand.py"


rule segment_exons_with_strand:
    input: config['references']['dmel']['db'].lstrip("../")
    output: "output/stranded_exon_segments.gtf"
    script: "scripts/segment_exons_with_strand.py"


rule fuse_exons_ignore_strand:
    input: config['references']['dmel']['db'].lstrip("../")
    output: "output/nonstranded_exon_fusions.gtf"
    script: "scripts/fuse_exons_ignore_strand.py"


rule fuse_exons_with_strand:
    input: config['references']['dmel']['db'].lstrip("../")
    output: "output/stranded_exon_fusions.gtf"
    script: "scripts/fuse_exons_with_strand.py"
