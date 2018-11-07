#!/usr/bin/env python
"""Constructs a geo submission."""
import os
import sys
from pathlib import Path
import shutil
import hashlib
from json import dumps

import pandas as pd
from pymongo import MongoClient
from dask import delayed
from dask.distributed import Client, LocalCluster

from ncbi_remap.logging import logger
from ncbi_remap.config import memory


DEBUG = False
OUTDIR = Path('../output/geo-wf/justin.fear@nih.gov')
OUTDIR.mkdir(parents=True, exist_ok=True)

if DEBUG:
    import logging
    logger.setLevel(logging.DEBUG)


def start_cluster():
    # Start up mini cluster
    cpus = int(os.environ['SLURM_CPUS_PER_TASK'])
    mem = int(os.environ['SLURM_MEM_PER_NODE']) / 1024 / cpus
    logger.info(f'Starting Mini Cluster: Cpu={cpus} Mem={mem}')

    cluster = LocalCluster(n_workers=cpus, threads_per_worker=1,
                           memory_limit=f'{mem}GB')

    client = Client(cluster)
    with open('dask_info.json', 'w') as fh:
        fh.write(dumps(client._scheduler_identity, indent=True))

    return client


def get_samples():
    logger.info('Connecting to store and getting list of complete samples')
    store = pd.HDFStore('../sra.h5', mode='r')
    samples = store['aln/complete'].srx.unique().tolist()
    store.close()
    return samples


def mongo_connect():
    try:
        with open('../output/.mongodb_host', 'r') as fh:
            host = fh.read().strip()
    except FileNotFoundError:
        host = 'localhost'
    client = MongoClient(host=host, port=27017)
    db = client['sra']
    return client, db['ncbi']


def get_labeled(samples):
    logger.info('Connecting to DB and getting list of RNA-Seq samples')
    client, ncbi = mongo_connect()
    labels = pd.DataFrame(list(ncbi.aggregate([
        {
            '$match': {
                '_id': {'$in': samples}
            }
        },
        {
            '$project': {
                '_id': 0,
                'srx': '$_id',
                'library_strategy': '$sra.experiment.library_strategy'
            }
        }

    ])))
    labels.set_index('srx', inplace=True)
    client.close()
    return labels


@memory.cache
def md5sum(filename, blocksize=65536):
    hash = hashlib.md5()
    with open(filename, "rb") as f:
        for block in iter(lambda: f.read(blocksize), b""):
            hash.update(block)
    return hash.hexdigest()


@delayed
def copy_files(files):
    hashes = []
    for k, v in files.items():
        fname = v['fname']
        ftype = v['ftype']

        # GEO wants us to use plus/minus instead of first/second
        if ftype == 'BigWig':
            _name = fname.name.replace('first', 'plus').replace('second', 'minus')
        else:
            _name = fname.name

        new = Path(OUTDIR, _name)

        # If the file exists then skip
        if new.exists():
            continue

        # Calculate hashes
        _hash = md5sum(fname)
        hashes.append((_name, ftype, _hash))

        # If DEBUG then don't copy files
        if DEBUG:
            continue

        shutil.copy(fname, new)
        hashNew = md5sum(new)

        try:
            assert _hash == hashNew
        except AssertionError:
            new.unlink()
            logger.warn('Error copying {fname.name}')

    return hashes


def main():
    samples = get_samples()
    labels = get_labeled(samples)
    rnaseq = labels.query('library_strategy == "RNA-Seq"').index.tolist()

    if DEBUG:
        rnaseq = rnaseq[:20]

    logger.info(f'Processing {len(rnaseq):,} SRXs')
    work = []
    for srx in rnaseq:
        files = {
            'firstFB': {
                'fname': Path(f'../output/aln-wf/samples/{srx}/{srx}.flybase.first.bw'),
                'ftype': 'BigWig',
            },
            'secondFB': {
                'fname': Path(f'../output/aln-wf/samples/{srx}/{srx}.flybase.second.bw'),
                'ftype': 'BigWig',
            },
            'gene': {
                'fname': Path(f'../output/aln-wf/samples/{srx}/{srx}.bam.counts'),
                'ftype': 'abundance measurements',
            },
            'geneJunc': {
                'fname': Path(f'../output/aln-wf/samples/{srx}/{srx}.bam.counts.jcounts'),
                'ftype': 'abundance measurements',
            },
            'inter': {
                'fname': Path(f'../output/aln-wf/samples/{srx}/{srx}.bam.intergenic.counts'),
                'ftype': 'abundance measurements',
            },
            'interJunc': {
                'fname': Path(f'../output/aln-wf/samples/{srx}/{srx}.bam.intergenic.counts.jcounts'),
                'ftype': 'abundance measurements',
            },
        }
        work.append(copy_files(files))

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
    df.to_csv(Path(OUTDIR, 'md5sum.tsv', mode='a'), sep='\t')


if __name__ == '__main__':
    try:
        main()
        logger.info('Script complete')
    except KeyboardInterrupt:
        logger.error('Keyboard interrupted.')
