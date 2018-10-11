import numpy as np
import pandas as pd
#from dask.distributed import LocalCluster, Client
#import dask.dataframe as dd

#threads = int(snakemake.threads)
#mem = int(snakemake.resources.mem_gb)
#mem -= mem * .05
#mem = int(np.floor(mem))


def main(client):
    with open('../output/rna_seq_srxs.txt') as fh:
        rnaseq_srx = fh.read().split('\n')

    cnts = pd.read_parquet('output/gene_counts_wide.parquet', columns=rnaseq_srx)
    corr = cnts.corr(method='spearman')
    corr.to_parquet(snakemake.output[0])
    #corr = cnts.corr().compute()
    #corr.to_parquet(snakemake.output[0])


if __name__ == '__main__':
#    cluster = LocalCluster(n_workers=int(threads / 2), threads_per_worker=2, memory_limit=f'{mem}GB')
#    client = Client(cluster)
    client = ''
    main(client)

