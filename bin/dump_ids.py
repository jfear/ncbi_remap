#!/usr/bin/env python
import os
import sys
sys.path.insert(0, '../lib/python')

import pandas as pd
import argparse

from pymongo import MongoClient

from ncbi_remap.logging import logger
from ncib_remap.io import build_index


def getOptions():
    """Function to pull in command line arguments"""
    parser = argparse.ArgumentParser(description='Tool to pull ids out of sramongo database and make a table.')
    parser.add_argument('--host', dest='host', action='store', default='localhost', required=False, help='Hostname with mongo database. [localhost]')
    parser.add_argument('--port', dest='port', action='store', default=27022, type=int, required=False, help='Hostname with mongo database. [27022]')
    parser.add_argument('--db', dest='db', action='store', default='sra2', required=False, help='Name of the database. [sra2]')
    parser.add_argument('--collection', dest='collection', action='store', default='ncbi', required=False, help='Name of the collection. [ncbi]')
    parser.add_argument('--store', dest='store', action='store', default='../config/sra.h5', required=False,
                        help='Hostname with mongo database. [../output/sra.h5]')
    args = parser.parse_args()
    return args


def add_new_data(df, store, key, **kwargs):
    newDat = df[~df.srr.isin(store[key].srr)].copy()

    if newDat.shape[0] > 0:
        store.append(key, newDat, data_columns=True, index=False)
        build_index(store, 'ids', columns=['srx', 'srr'])
    else:
        logger.info('No New Data.')


if __name__ == '__main__':
    args = getOptions()

    # Connect to database
    logger.info('Connecting to {}@{}:{}'.format(args.db, args.host, args.port))
    client = MongoClient(host=args.host, port=args.port)
    db = client[args.db]
    ncbi = db[args.collection]

    logger.info('Querying Database.')
    df = pd.DataFrame(list(ncbi.aggregate([
        {'$unwind': '$sra.run'},
        {
            '$project': {
                '_id': 0,
                'srx': '$_id',
                'srr': '$sra.run.run_id'
            }
        }
    ])))

    df['flag_pre_aln_complete'] = False
    df['flag_aln_complete'] = False

    logger.info('Loading into HDF5 store: {}'.format(args.store))
    store = pd.HDFStore(args.store)
    add_new_data(df, store, 'ids')
    store.close()

    logger.info('Complete')
