from pathlib import Path
import numpy as np
import pandas as pd

from dask import delayed
import dask.dataframe as dd
from dask.distributed import LocalCluster, Client

threads = int(snakemake.threads)
mem = int(snakemake.resources.mem_gb)
mem -= mem * .05
mem = int(np.floor(mem))

chroms = ['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3L', 'chr3R', 'chr4', 'chrY']

@delayed
def parse_juncs(srx):
    fname = f'../aln-wf/output/junction_counts/{srx}.parquet'
    df = pd.read_parquet(fname).query(f'Site1_chr == {chroms}')
    if df.shape[0] == 0:
        return srx, 0

    return srx, df['count'].sum()


def main(client):
    pth = Path('../aln-wf/output/junction_counts')
    srxs = [x.stem for x in pth.glob('*.parquet')]
    work = map(parse_juncs, srxs)
    futures = client.compute(work)
    data = client.gather(futures)
    df = pd.DataFrame(data, columns=['srx', 'count']).set_index('srx')
    df.to_parquet(snakemake.output[0])


if __name__ == '__main__':
    cluster = LocalCluster(n_workers=int(threads / 2), threads_per_worker=2, memory_limit=f'{mem}GB')
    client = Client(cluster)
    main(client)
