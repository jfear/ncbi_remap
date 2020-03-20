import sys
from pathlib import Path
from collections import namedtuple

import pandas as pd

sys.path.insert(0, "../src")
from ncbi_remap.parser import parse_picardCollect_summary, parse_picardCollect_hist


PREALN_PATH = Path("../output/prealn-wf/samples")

Files = namedtuple("Files", "srr idx unstranded first second strand table genebody_coverage")


class PicardException(Exception):
    """Problem with Picard file."""


def main():
    for srr_pth in PREALN_PATH.glob("**/SRR*"):
        if not srr_pth.is_dir():
            continue
        srr = srr_pth.name
        files = Files(
            srr,
            pd.Index([srr], name="srr"),
            (srr_pth / f"{srr}.hisat2.bam.NONE.picard.collectrnaseqmetrics"),
            (
                srr_pth
                / f"{srr}.hisat2.bam.FIRST_READ_TRANSCRIPTION_STRAND.picard.collectrnaseqmetrics"
            ),
            (
                srr_pth
                / f"{srr}.hisat2.bam.SECOND_READ_TRANSCRIPTION_STRAND.picard.collectrnaseqmetrics"
            ),
            f"../output/prealn-wf/strand/{srr}.parquet",
            f"../output/prealn-wf/rnaseqmetrics/{srr}.parquet",
            f"../output/prealn-wf/genebody_coverage/{srr}.parquet",
        )

        if not files.first.exists() & files.second.exists() & files.unstranded.exists():
            continue

        try:
            strand_flag(files)
            parse_table(files)
            genebody_coverage(files)
        except PicardException:
            remove_file(files.strand)
            remove_file(files.table)
            remove_file(files.genebody_coverage)
        
        remove_file(files.unstranded)
        remove_file(files.first)
        remove_file(files.second)


def genebody_coverage(files):
    try:
        # Parse genome coverage histogram
        df = parse_picardCollect_hist(files.unstranded)
        df.index = files.idx
        df.to_parquet(files.genebody_coverage)
    except AttributeError:
        raise PicardException("Missing genebody coverage:\t%s" % files.unstranded)
    except Exception as err:
        print(files.unstranded)
        raise err


def parse_stranded(file_name):
    return (parse_picardCollect_summary(file_name).PCT_CORRECT_STRAND_READS >= 0.75)[0]


def strand_flag(files):
    try:
        # Parse the flags
        if parse_stranded(files.first):
            strand = "same_strand"
        elif parse_stranded(files.second):
            strand = "opposite_strand"
        else:
            strand = "unstranded"

        df = pd.DataFrame([[strand]], index=files.idx, columns=["strand"])
        df.to_parquet(files.strand)
    except FileNotFoundError:
        raise PicardException("Missing files:\t%s, %s" % files.first, files.second)
    except Exception as err:
        print(files.first)
        raise err


def parse_table(files):
    try:
        df = parse_picardCollect_summary(files.unstranded)[
            [
                "PCT_CODING_BASES",
                "PCT_UTR_BASES",
                "PCT_INTRONIC_BASES",
                "PCT_INTERGENIC_BASES",
                "PCT_MRNA_BASES",
                "MEDIAN_CV_COVERAGE",
                "MEDIAN_5PRIME_BIAS",
                "MEDIAN_3PRIME_BIAS",
                "MEDIAN_5PRIME_TO_3PRIME_BIAS",
            ]
        ].fillna(0.0)
        df.columns = [x.lower() for x in df.columns]
        df.columns = [x.replace("pct_", "percent_") for x in df.columns]
        df.index = files.idx
        df.to_parquet(files.table)
    except FileNotFoundError:
        raise PicardException("Missing files: %s" % files.unstranded)
    except Exception as err:
        print(files.unstranded)
        raise err


def remove_file(file_name):
    if Path(file_name).exists():
        Path(file_name).unlink()


if __name__ == "__main__":
    main()
