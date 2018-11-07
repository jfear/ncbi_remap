import numpy as np
import dask.dataframe as dd
from dask.distributed import LocalCluster, Client

threads = int(snakemake.threads)
mem = int(snakemake.resources.mem_gb)
mem -= mem * .05
mem = int(np.floor(mem))


def main(client):
    ddf = dd.read_parquet('../output/aln-wf/intergenic_counts/*.parquet').groupby('srx').sum()
    ddf.to_parquet(snakemake.output[0], write_index=True)


if __name__ == '__main__':
    cluster = LocalCluster(n_workers=int(threads / 2), threads_per_worker=2, memory_limit=f'{mem}GB')
    client = Client(cluster)
    main(client)
