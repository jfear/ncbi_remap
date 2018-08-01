"""Calculate md5sums and copy files for RNA-Seq samples."""
from pathlib import Path

import pandas as pd

from ncbi_remap.logging import logger
from ncbi_remap.sample_lists import get_complete_rnaseq

from .md5sum import start_cluster, copy_files

DEBUG = False
OUTDIR = Path('output/justin.fear@nih.gov')
OUTDIR.mkdir(parents=True, exist_ok=True)

if DEBUG:
    import logging
    logger.setLevel(logging.DEBUG)


def main():
    rnaseq = get_complete_rnaseq()

    if DEBUG:
        rnaseq = rnaseq[:20]

    logger.info(f'Processing {len(rnaseq):,} SRXs')
    work = []
    for srx in rnaseq:
        files = {
            'firstFB': {
                'fname': Path(f'../aln-wf/output/samples/{srx}/{srx}.flybase.first.bw'),
                'ftype': 'BigWig',
            },
            'secondFB': {
                'fname': Path(f'../aln-wf/output/samples/{srx}/{srx}.flybase.second.bw'),
                'ftype': 'BigWig',
            },
            'gene': {
                'fname': Path(f'../aln-wf/output/samples/{srx}/{srx}.bam.counts'),
                'ftype': 'abundance measurements',
            },
            'geneJunc': {
                'fname': Path(f'../aln-wf/output/samples/{srx}/{srx}.bam.counts.jcounts'),
                'ftype': 'abundance measurements',
            },
            'inter': {
                'fname': Path(f'../aln-wf/output/samples/{srx}/{srx}.bam.intergenic.counts'),
                'ftype': 'abundance measurements',
            },
            'interJunc': {
                'fname': Path(f'../aln-wf/output/samples/{srx}/{srx}.bam.intergenic.counts.jcounts'),
                'ftype': 'abundance measurements',
            },
        }
        work.append(copy_files(files, OUTDIR, DEBUG))

    try:
        client = start_cluster()
        logger.info('Submitting to Mini Cluster')
        futures = client.compute(work)
        results = client.gather(futures)
    except KeyboardInterrupt as e:
        raise e
    finally:
        client.close()

    hashes = []
    for res in results:
        hashes.extend(res)

    df = pd.DataFrame(hashes, columns=['file name', 'file type', 'file checksum']).set_index('file name')

    if DEBUG:
        logger.debug("\n\n" + df.to_string() + "\n\n")
        return

    logger.info('Writing out hash table')
    df.to_csv(Path(OUTDIR, 'md5sum.tsv'), sep='\t')


if __name__ == '__main__':
    try:
        main()
        logger.info('Script complete')
    except KeyboardInterrupt:
        logger.error('Keyboard interrupted.')
