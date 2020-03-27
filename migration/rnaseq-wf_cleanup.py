import os
import re
import shutil
from pathlib import Path

CLEAN_UP = os.environ.get("CLEAN_UP", False)
SRR_PATTERN = re.compile(r"^[SED]RR\d+$")
TARGETS = [
    "../output/rnaseq-wf/aln_stats/{srx}.parquet",
    "../output/rnaseq-wf/gene_counts/{srx}.parquet",
    "../output/rnaseq-wf/junction_counts/{srx}.parquet",
    "../output/rnaseq-wf/intergenic_counts/{srx}.parquet",
    "../output/rnaseq-wf/segment_counts/{srx}.parquet",
    "../output/rnaseq-wf/fusion_counts/{srx}.parquet",
    "../output/rnaseq-wf/flybase_bigwigs/{srx}.flybase.first.bw",
    "../output/rnaseq-wf/flybase_bigwigs/{srx}.flybase.second.bw",
    "../output/rnaseq-wf/ucsc_bigwigs/{srx}.first.bw",
    "../output/rnaseq-wf/ucsc_bigwigs/{srx}.second.bw",
    "../output/rnaseq-wf/samples/{srx}/{srx}.bam",
    "../output/rnaseq-wf/samples/{srx}/{srx}.bam.bai",
]


def main():
    for srx_path in Path("../output/rnaseq-wf/samples").iterdir():
        srx = srx_path.name
        remove_temp(srx)

        if (
            Path(f"../output/rnaseq-wf/atropos_bad/{srx}").exists()
            or Path(f"../output/rnaseq-wf/alignment_bad/{srx}").exists()
        ):
            remove_srx_folder(srx)
            continue

        if all(check_target(target.format(srx=srx)) for target in TARGETS):
            Path(f"../output/rnaseq-wf/done/{srx}").touch()
            remove_srr_folders(srx)
            remove_processed_files(srx)
            remove_misc_files(srx)


def remove_temp(srx: str):
    for pth in Path(f"../output/rnaseq-wf/samples/{srx}").glob("*.tmp"):
        pth.unlink()


def remove_srx_folder(srx: str):
    pth = Path(f"../output/rnaseq-wf/samples/{srx}")
    if pth.exists() and CLEAN_UP:
        shutil.rmtree(pth)
    elif pth.exists():
        print("Removing SRX Folder:", pth, sep="\t")


def check_target(file_name: str):
    if Path(file_name).exists():
        return True
    print("Missing Target:", file_name, sep="\t")


def remove_srr_folders(srx: str):
    for pth in Path(f"../output/rnaseq-wf/samples/{srx}").iterdir():
        if pth.is_dir() and re.match(SRR_PATTERN, pth.name):
            if CLEAN_UP:
                shutil.rmtree(pth)
            else:
                print("Removing SRR Folder:", pth, sep="\t")


def remove_file(file_name: str):
    pth = Path(file_name)
    if pth.exists() and CLEAN_UP:
        pth.unlink()
    elif pth.exists():
        print("Removing File:", pth, sep="\t")


def remove_processed_files(srx: str):
    remove_file(f"../output/rnaseq-wf/samples/{srx}/{srx}.bam.samtools.stats")
    remove_file(f"../output/rnaseq-wf/samples/{srx}/{srx}.bam.bamtools.stats")

    remove_file(f"../output/rnaseq-wf/samples/{srx}/{srx}.bam.counts")
    remove_file(f"../output/rnaseq-wf/samples/{srx}/{srx}.bam.counts.jcounts")
    remove_file(f"../output/rnaseq-wf/samples/{srx}/{srx}.bam.intergenic.counts")
    remove_file(f"../output/rnaseq-wf/samples/{srx}/{srx}.bam.exon_fusions.counts")
    remove_file(f"../output/rnaseq-wf/samples/{srx}/{srx}.bam.exon_segments.counts")


def remove_misc_files(srx: str):
    remove_file(f"../output/rnaseq-wf/samples/{srx}/{srx}.trim.clean.tsv")
    remove_file(f"../output/rnaseq-wf/samples/{srx}/{srx}.hisat2.bam.tsv")
    remove_file(f"../output/rnaseq-wf/samples/{srx}/{srx}.bam.samtools.idxstats")
    remove_file(f"../output/rnaseq-wf/samples/{srx}/{srx}.counts")

    remove_file(f"../output/rnaseq-wf/samples/{srx}/{srx}.bam.counts.summary")
    remove_file(f"../output/rnaseq-wf/samples/{srx}/{srx}.bam.counts.log")

    remove_file(f"../output/rnaseq-wf/samples/{srx}/{srx}.bam.intergenic.counts.jcounts")
    remove_file(f"../output/rnaseq-wf/samples/{srx}/{srx}.bam.intergenic.counts.summary")
    remove_file(f"../output/rnaseq-wf/samples/{srx}/{srx}.bam.intergenic.counts.log")

    remove_file(f"../output/rnaseq-wf/samples/{srx}/{srx}.bam.exon_fusions.counts.jcounts")
    remove_file(f"../output/rnaseq-wf/samples/{srx}/{srx}.bam.exon_fusions.counts.summary")
    remove_file(f"../output/rnaseq-wf/samples/{srx}/{srx}.bam.exon_fusions.counts.log")

    remove_file(f"../output/rnaseq-wf/samples/{srx}/{srx}.bam.exon_segments.counts.jcounts")
    remove_file(f"../output/rnaseq-wf/samples/{srx}/{srx}.bam.exon_segments.counts.summary")
    remove_file(f"../output/rnaseq-wf/samples/{srx}/{srx}.bam.exon_segments.counts.log")


if __name__ == "__main__":
    main()
