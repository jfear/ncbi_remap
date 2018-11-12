#!/usr/bin/env python
"""Calculate the pairwise correlations among samples."""
import os
import sys
from pathlib import Path
from json import dumps

import pandas as pd

from dask import delayed
from dask.distributed import Client, LocalCluster

sys.path.insert(0, '../lib')
from ncbi_remap.logging import logger
from ncbi_remap.normalization import cpm


DEBUG = False
DEBUG_N = 1000


def start_cluster():
    # Start up mini cluster
    cpus = int(os.environ['SLURM_CPUS_PER_TASK'])
    mem = int(os.environ['SLURM_MEM_PER_NODE']) / 1024 / cpus
    logger.info(f'Starting Mini Cluster: Cpu={cpus} Mem={mem}')

    cluster = LocalCluster(n_workers=cpus, threads_per_worker=1,
                           memory_limit=f'{mem}GB')

    client = Client(cluster)
    with open('../output/downstream-analysis/dask_info.json', 'w') as fh:
        fh.write(dumps(client._scheduler_identity, indent=True))

    return client


def get_samples():
    # Connect to data store
    store = pd.HDFStore('../sra.h5', mode='r')
    samples = store['aln/complete'].srx.unique().tolist()
    store.close()
    logger.info(f'{len(samples):,} samples loaded')
    return samples[::-1]


@delayed
def get_norm_counts(srx):
    data = []
    # Gene level counts
    gene = pd.read_parquet(f'../output/aln-wf/gene_counts/{srx}.parquet')
    gene.reset_index(inplace=True)
    gene[['FBgn', 'srx', 'count']].copy()
    gene['var_type'] = 'gene'
    if gene.size > 0:
        data.append(gene)

    # junction level counts
    chroms = ['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrY']
    junc = pd.read_parquet(f'../output/aln-wf/junction_counts/{srx}.parquet')
    junc = junc.query(f'Site1_chr == {chroms} & Site1_chr == Site2_chr')
    junc = junc[['PrimaryGene', 'srx', 'count']].copy()
    junc.columns = ['FBgn', 'srx', 'count']
    junc.dropna(inplace=True)
    junc = junc.groupby(['FBgn', 'srx'])['count'].sum().reset_index()
    junc['var_type'] = 'junction'
    if junc.size > 0:
        data.append(junc)

    # intergenic counts
    inter_annot = pd.read_csv(
        '../output/dmel_r6-11.intergenic.bed',
        sep='\t',
        header=None,
        names=['chrom', 'start', 'end', 'FBgn'],
        index_col='FBgn'
    )
    inter_names = inter_annot.query(f'chrom == {chroms}').index.unique().tolist()
    inter = pd.read_parquet(f'../output/aln-wf/intergenic_counts/{srx}.parquet')
    inter = inter.query(f'FBgn == {inter_names}')
    inter.reset_index(inplace=True)
    inter = inter[['FBgn', 'srx', 'count']].copy()
    inter['var_type'] = 'intergenic'
    if inter.size > 0:
        data.append(inter)

    # combine, normalize, aggregate
    df = pd.concat(data, sort=True)
    norm = cpm(df.set_index(['FBgn', 'srx', 'var_type']), log='log10').reset_index()
    medians = norm.groupby(['srx', 'var_type']).median().unstack()
    medians.columns = medians.columns.droplevel(0)

    # If missing any type set to 0
    cols = medians.columns
    if 'gene' not in cols:
        medians['gene'] = 0.0
    elif 'junction' not in cols:
        medians['junction'] = 0.0
    elif 'intergenic' not in cols:
        medians['intergenic'] = 0.0

    return medians


def main():
    samples = get_samples()
    srxs = samples

    if DEBUG:
        logger.info(f'DEBUG ON: only running {DEBUG_N:,} samples')
        srxs = samples[:DEBUG_N]

    client = start_cluster()
    try:
        logger.info('Submitting jobs to Mini Cluster')
        futures = client.compute([get_norm_counts(srx) for srx in srxs])

        logger.info('Gathering results')
        df = pd.concat(client.gather(futures))
    except KeyboardInterrupt as e:
        raise e
    finally:
        logger.info('Shutting down Mini Cluster')
        client.close()

    logger.info('Saving results')
    df.to_parquet('../output/downstream-analysis/counts_norm_median.parquet')


if __name__ == '__main__':
    try:
        Path('../output/downstream-analysis').mkdir(parents=True, exist_ok=True)
        main()
    except KeyboardInterrupt:
        logger.info('Interrupted')
    else:
        logger.info('Script Complete')
