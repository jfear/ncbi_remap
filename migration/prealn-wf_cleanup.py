from pathlib import Path
from collections import namedtuple
import shutil


PREALN_PATH = Path("../output/prealn-wf/samples")

Files = namedtuple(
    "Files",
    "srr path alignment_bad aln_stats atropos atropos_bad count_summary fastq_screen genebody_coverage hisat2 markduplciates rnaseqmetrics strand done",
)


def main():
    for srr_pth in PREALN_PATH.glob("**/SRR*"):
        if not srr_pth.is_dir():
            continue
        srr = srr_pth.name
        files = Files(
            srr,
            srr_pth,
            Path(f"../output/prealn-wf/alignment_bad/{srr}.parquet"),
            Path(f"../output/prealn-wf/aln_stats/{srr}.parquet"),
            Path(f"../output/prealn-wf/atropos/{srr}.parquet"),
            Path(f"../output/prealn-wf/atropos_bad/{srr}.parquet"),
            Path(f"../output/prealn-wf/count_summary/{srr}.parquet"),
            Path(f"../output/prealn-wf/fastq_screen/{srr}.parquet"),
            Path(f"../output/prealn-wf/genebody_coverage/{srr}.parquet"),
            Path(f"../output/prealn-wf/hisat2/{srr}.parquet"),
            Path(f"../output/prealn-wf/markduplicates/{srr}.parquet"),
            Path(f"../output/prealn-wf/rnaseqmetrics/{srr}.parquet"),
            Path(f"../output/prealn-wf/strand/{srr}.parquet"),
            Path(f"../output/prealn-wf/done/{srr}"),
        )

        if files.alignment_bad.exists() or files.atropos_bad.exists():
            shutil.rmtree(files.path)

        if all(
            pth.exists()
            for pth in [
                files.atropos,
                files.fastq_screen,
                files.hisat2,
                files.aln_stats,
                files.strand,
                files.rnaseqmetrics,
                files.genebody_coverage,
                files.markduplciates,
                files.count_summary,
            ]
        ):
            files.done.touch()
            shutil.rmtree(files.path)


if __name__ == "__main__":
    main()
