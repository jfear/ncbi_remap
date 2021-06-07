import gzip
import os
import shutil
from pathlib import Path
from subprocess import run
from typing import List

import pandas as pd


def main(
    genic: Path,
    intergenic: Path,
    fusion: Path,
    segment: Path,
    bams: List[Path],
    layouts: List[Path],
    strands: List[Path],
    done_flag: Path,
):
    fc = FeatureCountsRunner(genic, intergenic, fusion, segment)

    for bam, layout, strand in zip(bams, layouts, strands):
        fc.run(bam, pd.read_parquet(layout).layout[0], pd.read_parquet(strand).strand[0])

    done_flag.write_text("\n".join(bam.stem for bam in bams) + "\n")


class FeatureCountsException(Exception):
    """Fastq Screen Processing Exception"""


class FeatureCountsRunner:
    fc_options = {
        "genic": "",
        "intergenic": "-t gene",
        "fusion": "-t fusion -g ID",
        "segment": "-t segment -g ID",
        "PE": "-p -P -C",
        "same_strand": "-s 1",
        "opposite_strand": "-s 2",
        "unstranded": "-s 0",
    }

    file_patterns = {
        "genic": {
            "counts": "{prefix}/rnaseq-counts-wf/gene_counts/{srx}.bam.counts.txt.gz",
            "jcounts": "{prefix}/rnaseq-counts-wf/gene_counts/{srx}.bam.counts.jcounts.txt.gz",
        },
        "intergenic": {
            "counts": "{prefix}/rnaseq-counts-wf/intergenic_counts/{srx}.bam.intergenic.counts.txt.gz",
            "jcounts": "{prefix}/rnaseq-counts-wf/intergenic_counts/{srx}.bam.intergenic.counts.jcounts.txt.gz",
        },
        "fusion": {
            "counts": "{prefix}/rnaseq-counts-wf/fusion_counts/{srx}.bam.fusion.counts.txt.gz",
            "jcounts": "{prefix}/rnaseq-counts-wf/fusion_counts/{srx}.bam.fusion.counts.jcounts.txt.gz",
        },
        "segment": {
            "counts": "{prefix}/rnaseq-counts-wf/segment_counts/{srx}.bam.segment.counts.txt.gz",
            "jcounts": "{prefix}/rnaseq-counts-wf/segment_counts/{srx}.bam.segment.counts.jcounts.txt.gz",
        },
    }

    def __init__(self, genic, intergenic, fusion, segment):
        self.gtfs = {
            "genic": self.stage_file(genic),
            "intergenic": self.stage_file(intergenic),
            "fusion": self.stage_file(fusion),
            "segment": self.stage_file(segment),
        }

    def run(self, bam: Path, layout: str, strand: str):
        self.srx = bam.stem
        self.prefix = bam.resolve().parents[3]
        self.local_bam = self.stage_file(bam)
        self.layout = layout
        self.strand = strand

        for gtf in self.gtfs.keys():
            self.feature_counts(gtf)
            self.save_outputs(gtf)
            self.clean_up()

        self.local_bam.unlink()

    def feature_counts(self, gtf: str, attempt: int = 0):
        log = self.local_bam.with_suffix(".log")
        cmd = " ".join(
            [
                "featureCounts",
                f"-T {THREADS}",
                self.fc_options.get(self.layout, ""),
                "-J",
                self.fc_options[gtf],
                self.fc_options[self.strand],
                f"-a {self.gtfs[gtf]}",
                f"-o {self.local_bam.with_suffix('.counts')}",
                self.local_bam.as_posix(),
                f"> {log.as_posix()} 2>&1",
            ]
        )

        print(cmd)
        run(cmd, shell=True, check=True)
        try:
            self.check_log(log)
            log.unlink()
        except FeatureCountsException as err:
            if attempt < 1:
                self.feature_counts(gtf, attempt + 1)
            else:
                print(f"Cannot get {gtf} counts on: {self.srx}")
                raise err

    def save_outputs(self, gtf: str):
        local_path = self.local_bam.with_suffix(".counts")
        remote_path = Path(
            self.file_patterns[gtf]["counts"].format(prefix=self.prefix, srx=self.srx)
        )
        remote_path.parent.mkdir(exist_ok=True, parents=True)
        self.copy_and_zip(local_path, remote_path)

        local_path = self.local_bam.with_suffix(".counts.jcounts")
        remote_path = Path(
            self.file_patterns[gtf]["jcounts"].format(prefix=self.prefix, srx=self.srx)
        )
        remote_path.parent.mkdir(exist_ok=True, parents=True)
        self.copy_and_zip(local_path, remote_path)

    def clean_up(self):
        self.local_bam.with_suffix(".counts").unlink()
        self.local_bam.with_suffix(".counts.jcounts").unlink()
        self.local_bam.with_suffix(".counts.summary").unlink()

    @staticmethod
    def copy_and_zip(local: Path, target: Path):
        with local.open("rb") as fh, gzip.open(target, "wb") as fo:
            shutil.copyfileobj(fh, fo)

    @staticmethod
    def check_log(log: Path) -> None:
        log_text = log.read_text()
        if (
            "Read assignment finished" not in log_text
            or "Summary of counting results" not in log_text
        ):
            raise FeatureCountsException(f"{log.stem} not complete")

    @staticmethod
    def stage_file(file_name: Path) -> Path:
        tmp_path = TMPDIR / file_name.name
        shutil.copyfile(file_name, tmp_path)
        return tmp_path


def get_tmp_dir() -> Path:
    job_id = os.getenv("SLURM_JOBID", None)

    if job_id is None:
        return Path("/tmp")

    if Path("/lscratch", job_id).exists():
        return Path("/lscratch", job_id)

    if (Path("~/scratch/lscratch").expanduser() / job_id).exists():
        return Path("~/scratch/lscratch").expanduser() / job_id

    Path("/tmp", job_id).mkdir(exist_ok=True)
    return Path("/tmp", job_id)


if __name__ == "__main__" and "snakemake" in locals():
    TMPDIR = get_tmp_dir()
    THREADS = str(snakemake.threads)  # type: ignore # noqa

    main(
        genic=Path(snakemake.input.genic),  # type: ignore # noqa
        intergenic=Path(snakemake.input.intergenic),  # type: ignore # noqa
        fusion=Path(snakemake.input.fusion),  # type: ignore # noqa
        segment=Path(snakemake.input.segment),  # type: ignore # noqa
        bams=[Path(x) for x in snakemake.input.bams],  # type: ignore # noqa
        layouts=[Path(x) for x in snakemake.input.layouts],  # type: ignore # noqa
        strands=[Path(x) for x in snakemake.input.strands],  # type: ignore # noqa
        done_flag=Path(snakemake.output[0]),  # type: ignore # noqa
    )
