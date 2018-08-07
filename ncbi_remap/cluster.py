"""Functions necessary for handling compute clusters.

"""

import os
from tempfile import NamedTemporaryFile
from json import dumps

from dask.distributed import Client, LocalCluster

from ncbi_remap.logging import logger


def start_dask_cluster(cpus=None, mem=None):
    if cpus is None:
        _cpus = os.environ.get('SLURM_CPUS_PER_TASK', 4)
        cpus = int(_cpus)

    if mem is None:
        _mem = os.environ.get('SLURM_MEM_PER_NODE', 6180)
        mem = int(_mem) / 1024 / cpus

    logger.info(f'Starting Mini Cluster: Cpu={cpus} Mem={mem}')

    cluster = LocalCluster(n_workers=cpus, threads_per_worker=1,
                           memory_limit=f'{mem}GB')

    client = Client(cluster)

    log = NamedTemporaryFile('wt', prefix='dask_info.', suffix='.json', delete=False)
    logger.info(f'Writing dask JSON dump to {log.name}')
    log.write(dumps(client._scheduler_identity, indent=True))
    log.close()

    return client
