"""Helper functions to get various sample lists."""
import pandas as pd

from ncbi_remap.logging import logger
from ncbi_remap.mongo import mongo_connect


def get_complete():
    logger.info('Connecting to store and getting list of complete samples')
    store = pd.HDFStore('../output/sra.h5', mode='r')
    srxs = store['aln/complete'].srx.unique().tolist()
    store.close()
    return srxs


def get_complete_rnaseq():
    samples = get_complete()

    logger.info('Connecting to DB and getting list of RNA-Seq samples')
    client, db, ncbi = mongo_connect()
    srxs = [x['_id'] for x in ncbi.aggregate([
        {
            '$match': {
                '_id': {'$in': samples},
                'sra.experiment.library_strategy': 'RNA-Seq',
            }
        },
        {
            '$project': {
                '_id': 1,
            }
        }

    ])]

    client.close()

    return srxs


def get_complete_chipseq():
    samples = get_complete()

    logger.info('Connecting to DB and getting list of RNA-Seq samples')
    client, db, ncbi = mongo_connect()
    srxs = [x['_id'] for x in ncbi.aggregate([
        {
            '$match': {
                '_id': {'$in': samples},
                'sra.experiment.library_strategy': 'ChIP-Seq',
            }
        },
        {
            '$project': {
                '_id': 1,
            }
        }

    ])]

    client.close()

    return srxs
