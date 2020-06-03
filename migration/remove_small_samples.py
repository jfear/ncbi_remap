"""Remove samples with a small number of reads.

I just added a 100k read minimum for downloading a sample. Samples can fall
below this after trimming, but I want to go ahead and remove all samples that
I have now filtered out.

"""
import textwrap
from pathlib import Path
import shutil

import pandas as pd


PROJECT_DIR = Path(__file__).parents[1]
OUTPUT_DIR = PROJECT_DIR / "output"


def main():
    srx2srr = pd.read_csv(OUTPUT_DIR / "srx2srr.csv")
    ok_srxs = srx2srr.srx.unique()
    ok_srrs = srx2srr.srr.unique()

    clean_up_fastq_wf(ok_srrs)
    clean_up_prealn_wf(ok_srxs, ok_srrs)
    clean_up_rnaseq_wf(ok_srxs, ok_srrs)


def clean_up_fastq_wf(srrs):
    pth = OUTPUT_DIR / "fastq-wf"
    iterdir(Path(pth, "abi_solid"), srrs)
    iterdir(Path(pth, "download_bad"), srrs)
    iterdir(Path(pth, "fastqs"), srrs)
    iterdir(Path(pth, "layout"), srrs)
    iterdir(Path(pth, "libsize"), srrs)
    iterdir(Path(pth, "sra_cache"), srrs)
    iterdir(Path(pth, "sra_download_logs"), srrs)


def clean_up_prealn_wf(srxs, srrs):
    pth = OUTPUT_DIR / "prealn-wf"
    iterdir(Path(pth, "alignment_bad"), srrs)
    iterdir(Path(pth, "aln_stats"), srrs)
    iterdir(Path(pth, "atropos"), srrs)
    iterdir(Path(pth, "atropos_bad"), srrs)
    iterdir(Path(pth, "count_summary"), srrs)
    iterdir(Path(pth, "done"), srrs)
    iterdir(Path(pth, "fastq_screen"), srrs)
    iterdir(Path(pth, "genebody_coverage"), srrs)
    iterdir(Path(pth, "hisat2"), srrs)
    iterdir(Path(pth, "markduplicates"), srrs)
    iterdir(Path(pth, "revisit"), srrs)
    iterdir(Path(pth, "rnaseqmetrics"), srrs)
    iterdir(Path(pth, "samples"), srxs)
    iterdir_nested(Path(pth, "samples"), srrs)
    iterdir(Path(pth, "strand"), srrs)


def clean_up_rnaseq_wf(srxs, srrs):
    pth = OUTPUT_DIR / "rnaseq-wf"
    iterdir(Path(pth, "alignment_bad"), srrs)
    iterdir(Path(pth, "aln_stats"), srxs)
    iterdir(Path(pth, "atropos"), srrs)
    iterdir(Path(pth, "atropos_bad"), srrs)
    iterdir(Path(pth, "done"), srxs)
    iterdir(Path(pth, "flybase_bigwigs"), srxs)
    iterdir(Path(pth, "fusion_counts"), srxs)
    iterdir(Path(pth, "gene_counts"), srxs)
    iterdir(Path(pth, "hisat2"), srrs)
    iterdir(Path(pth, "intergenic_counts"), srxs)
    iterdir(Path(pth, "junction_counts"), srxs)
    iterdir(Path(pth, "samples"), srxs)
    iterdir_nested(Path(pth, "samples"), srrs)
    iterdir(Path(pth, "segment_counts"), srxs)
    iterdir(Path(pth, "ucsc_bigwigs"), srxs)


def iterdir(pth, ok_list):
    cnt = 0
    for pth_ in pth.iterdir():
        file_ = pth_.name.split(".")[0].split("_")[0]
        if file_ not in ok_list:
            cnt += 1
    try:
        print(f"{pth_.as_posix():<80}{file_}: {cnt:,}")
    except UnboundLocalError:
        pass


def iterdir_nested(pth, ok_list):
    cnt = 0
    for pth_ in pth.glob("*/*"):
        if pth_.is_file():
            continue
        folder_ = pth_.name.split(".")[0].split("_")[0]
        if folder_ not in ok_list:
            cnt += 1

    try:
        print(f"{pth_.as_posix():<80}{folder_}: {cnt:,}")
    except UnboundLocalError:
        pass


def remove_file(file_name):
    if Path(file_name).exists() & Path(file_name).is_file():
        Path(file_name).unlink()


def remove_folder(folder_name):
    if Path(folder_name).exists() & Path(folder_name).is_dir():
        shutil.rmtree(folder_name)


if __name__ == "__main__":
    main()
