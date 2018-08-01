"""Calculate md5sums and copy files for RNA-Seq samples."""
import os
from pathlib import Path
import shutil
import hashlib
from json import dumps

from dask import delayed
from dask.distributed import Client, LocalCluster

from ncbi_remap.logging import logger
from ncbi_remap.config import memory


def start_cluster():
    # Start up mini cluster
    cpus = int(os.environ['SLURM_CPUS_PER_TASK'])
    mem = int(os.environ['SLURM_MEM_PER_NODE']) / 1024 / cpus
    logger.info(f'Starting Mini Cluster: Cpu={cpus} Mem={mem}')

    cluster = LocalCluster(n_workers=cpus, threads_per_worker=1,
                           memory_limit=f'{mem}GB')

    client = Client(cluster)
    with open('../output/dask_info.json', 'w') as fh:
        fh.write(dumps(client._scheduler_identity, indent=True))

    return client


@memory.cache
def md5sum(filename, blocksize=65536):
    hash = hashlib.md5()
    with open(filename, "rb") as f:
        for block in iter(lambda: f.read(blocksize), b""):
            hash.update(block)
    return hash.hexdigest()


@delayed
def copy_files(files, outdir, debug=False):
    hashes = []
    for k, v in files.items():
        fname = v['fname']
        ftype = v['ftype']
        # Calculate hashes
        _hash = md5sum(fname)

        # GEO wants us to use plus/minus instead of first/second
        if ftype == 'BigWig':
            _name = fname.name.replace('first', 'plus').replace('second', 'minus')
        else:
            _name = fname.name

        hashes.append((_name, ftype, _hash))
        new = Path(outdir, _name)

        # Copy files
        if debug:
            continue

        if new.exists():
            continue

        shutil.copy(fname, new)
        hashNew = md5sum(new)

        try:
            assert _hash == hashNew
        except AssertionError:
            new.unlink()
            logger.warn('Error copying {fname.name}')

    return hashes
