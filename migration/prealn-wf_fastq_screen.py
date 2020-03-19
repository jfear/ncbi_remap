import sys
from pathlib import Path
from collections import namedtuple

import pandas as pd

sys.path.insert(0, "../src")
from ncbi_remap.parser import parse_fastq_screen


PREALN_PATH = Path("../output/prealn-wf/samples")

Files = namedtuple("Files", "srr idx table log output")


def main():
    for srr_pth in PREALN_PATH.glob("**/SRR*"):
        if not srr_pth.is_dir():
            continue
        srr = srr_pth.name
        files = Files(
            srr,
            pd.Index([srr], name="srr"),
            (srr_pth / f"{srr}_1.fastq_screen.txt"),
            (srr_pth / f"{srr}_1.fastq_screen.txt.log"),
            f"../output/prealn-wf/fastq_screen/{srr}.parquet",
        )

        if files.table.exists():
            convert_table(files)
            Path(files.table).unlink()
            Path(files.log).unlink()


def convert_table(files):
    try:
        df = parse_fastq_screen(files.table).set_index("reference").fillna(0)
        summarized = (
            (
                (df.one_hit_one_library_count + df.multiple_hits_one_library_count)
                / df.reads_processed_count
                * 100
            )
            .rename(files.srr)
            .rename_axis("")
            .to_frame()
            .T.rename_axis("srr")
        )
        summarized.columns = [f"{col}_pct_reads_mapped" for col in summarized.columns]
        summarized.to_parquet(files.output)
    except:
        pass


if __name__ == "__main__":
    main()
