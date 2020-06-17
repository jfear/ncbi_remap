import os
import re
import shutil
import sys
from pathlib import Path
from typing import Optional, Tuple

import pandas as pd
from snakemake.logging import logger
from snakemake.shell import shell

sys.path.insert(0, "../src")
from ncbi_remap.io import remove_file, remove_folder
from ncbi_remap.snakemake import StepLogger


LOG = StepLogger(str(snakemake.log))
SRR = snakemake.wildcards.srr
THREADS = snakemake.threads

TMPDIR = Path(os.getenv("TMPDIR", "/tmp"), f"{SRR}/atropos")
TMPDIR.mkdir(parents=True, exist_ok=True)


def main():
    r1, r2 = stage_data(snakemake.input.r1, snakemake.input.r2)
    log, r1_trim, r2_trim = atropos(snakemake.input.layout, r1, r2)
    summarize(log, snakemake.output.summary)
    save_trimmed_reads(r1_trim, r2_trim, snakemake.output.r1, snakemake.output.r2)


def stage_data(r1: str, r2: str) -> Tuple[Path, Path]:
    r1_local = TMPDIR / Path(r1).name
    shutil.copy2(r1, r1_local)

    r2_local = TMPDIR / Path(r2).name
    shutil.copy2(r2, r2_local)

    return r1_local, r2_local


def atropos(layout: str, r1: Path, r2: Path) -> Tuple[Path, Path, Path]:
    log = TMPDIR / "atropos.log"
    r1_trim = TMPDIR / f"{SRR}_1.trim.fastq.gz"
    r2_trim = TMPDIR / f"{SRR}_2.trim.fastq.gz"

    layout_ = pd.read_parquet(layout).layout[0]
    try:
        if layout_ == "PE":
            shell(
                f"atropos trim {snakemake.params.extra_pe} --threads {THREADS} "
                f"-pe1 {r1} -pe2 {r2} -o {r1_trim} -p {r2_trim} >{log} 2>&1"
            )
        elif layout_ == "keep_R2":
            r1_trim.touch()
            shell(
                f"atropos trim {snakemake.params.extra_se} --threads {THREADS} -se {r2} -o {r2_trim} >{log} 2>&1"
            )
        else:
            r2_trim.touch()
            shell(
                f"atropos trim {snakemake.params.extra_se} --threads {THREADS} -se {r1} -o {r1_trim} >{log} 2>&1"
            )
    finally:
        LOG.append("atropos", log)

    remove_file(r1)
    remove_file(r2)
    return log, r1_trim, r2_trim


def summarize(log: Path, output_file: str) -> None:
    df = pd.DataFrame(
        [parse_atropos_log(log)],
        index=pd.Index([SRR], name="srr"),
        columns=["total_processed", "total_written", "too_short"],
    )

    if df.total_written[0] < 1_000:
        raise AtroposException("<1,000 reads")

    df.to_parquet(output_file)


def save_trimmed_reads(r1: Path, r2: Path, r1_out: str, r2_out: str) -> None:
    shutil.copy2(r1, r1_out)
    shutil.copy2(r2, r2_out)


def parse_atropos_log(file_name: Path) -> Tuple[int, int, int]:
    """Example Atropos output
        # Example of SE
        Reads                                 records  fraction
        ----------------------------------- --------- ---------
        Total reads processed:                 85,319
        Reads with adapter:                     4,725      5.5%
        Reads that were too short:              1,879      2.2%
        Reads written (passing filters):       83,440     97.8%

        # Example of PE
        Pairs                                records fraction
        ----------------------------------- -------- --------
        Total read pairs processed:            1,354
        Read 1 with adapter:                   106     7.8%
        Read 2 with adapter:                   123     9.1%
        Pairs that were too short:                16     1.2%
        Pairs written (passing filters):       1,338    98.8%
    """
    with open(file_name) as fh:
        log_text = fh.read().replace(",", "")

    if "ERROR" in log_text:
        logger.warning(f"{SRR}: Atropos reported an error")
        LOG.append("Atropos Check", text="Atropos reported an error")

    try:
        tot_processed = int(re.findall(r"Total read.*processed:\s+(\d+)", log_text)[0])
        tot_written = int(
            re.findall(r"[Read|Pair]s written \(passing filters\):\s+(\d+)", log_text)[0]
        )
        too_short = int(re.findall(r"[Read|Pair]s that were too short:\s+(\d+)", log_text)[0])
        return tot_processed, tot_written, too_short
    except IndexError:
        raise AtroposException("Unable to parse log")


class AtroposException(Exception):
    """Basic Atropos Exception"""


if __name__ == "__main__":
    try:
        main()
    except AtroposException as error:
        logger.warning(f"Flagging {SRR} as Atropos Bad")
        LOG.append("Exception", text=str(error))

        # Add flag
        pth = Path(snakemake.params.atropos_bad)
        pth.mkdir(exist_ok=True)
        (pth / SRR).touch()

        # Remove outputs
        remove_folder(Path(snakemake.output.r1).parent)
        remove_file(snakemake.output.summary)
        
        raise SystemExit
    finally:
        remove_folder(TMPDIR)
