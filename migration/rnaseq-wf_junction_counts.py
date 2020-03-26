import os
import sys
from pathlib import Path
from collections import namedtuple
from multiprocessing import Pool

sys.path.insert(0, "../src")
from ncbi_remap.parser import parse_featureCounts_jcounts

CLEAN_UP = os.environ.get("CLEAN_UP", False)
THREADS = int(os.environ.get("SLURM_CPUS_PER_TASK", "2"))
Files = namedtuple("Files", "srx file_name output")

def main():
    workers = Pool(THREADS)
    work = []
    for srx_pth in Path("../output/rnaseq-wf/samples").iterdir():
        srx = srx_pth.name
        if not srx_pth.is_dir():
            continue

        files = Files(
            srx,
            srx_pth / f"{srx}.bam.counts.jcounts",
            f"../output/rnaseq-wf/junction_counts/{srx}.parquet"
        )

        if files.file_name.exists():
            work.append(files)

    workers.map(parse_srx, work)


def parse_srx(files):
    df = parse_featureCounts_jcounts(files.file_name, files.srx)
    df.to_parquet(files.output)
    if CLEAN_UP:
        files.file_name.unlink()


if __name__ == "__main__":
    main()
