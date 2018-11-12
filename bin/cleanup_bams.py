#!/usr/bin/env python
"""Remove BAMs from the alignment workflow that are no longer needed.

I have been keeping the BAMs from the alignment workflow b/c I wanted them for
RNA-Seq downstream analysis. Now I have solid calls about what is RNA-Seq/EST
and I want to go ahead and remove the other BAM files.

"""
from textwrap import dedent
from pathlib import Path
from logging import DEBUG as LDB
import pandas as pd

from ncbi_remap.logging import logger

DEBUG = False


class Tracker:
    def __init__(self, ftype=''):
        self.ftype = ftype
        self.count = 0
        self.size = 0

    def remove(self, pth):
        if not pth.exists():
            logger.debug('%s does not exists.', pth.as_posix())
            return

        # Track number of files and sizes
        self.count += 1
        self.size += pth.stat().st_size

        # Delete the file
        if DEBUG:
            logger.debug('This will delete: %s', pth.as_posix())
        else:
            pth.unlink()

    def __str__(self):
        return dedent(f"""
        There were {self.count:,} {self.ftype} files removed.
        For a total of {self.size / 1e9:,.4f} GB.
        """)


def main():
    # Get list of completed SRXs
    store = pd.HDFStore('../sra.h5', mode='r')
    complete = store['aln/complete'].srx.values
    srxs = not_rnaseq_srxs(complete)

    # Set up tracker
    tracker_bam = Tracker('BAM')
    tracker_bai = Tracker('BAI')

    # Remove files
    for srx in srxs:
        pth1 = Path(f'../output/aln-wf/samples/{srx}/{srx}.bam')
        tracker_bam.remove(pth1)

        pth2 = Path(f'../output/aln-wf/samples/{srx}/{srx}.bam.bai')
        tracker_bai.remove(pth2)

    # Report Summary
    logger.info(tracker_bam)
    logger.info(tracker_bai)


def not_rnaseq_srxs(complete):
    """Create a list of non RNA-Seq SRXs."""
    metadata = pd.read_parquet('../output/metadata-wf/select_library_strategy.parquet')
    metadata = metadata.reindex(complete)
    mask = metadata.Fear_et_al_library_strategy.str.contains('RNA-Seq') | \
        metadata.Fear_et_al_library_strategy.str.contains('EST')

    logger.debug('%d SRXs to be cleaned.', (~mask).sum())
    return metadata[~mask].index.tolist()


if __name__ == '__main__':
    if DEBUG:
        logger.setLevel(LDB)

    main()
