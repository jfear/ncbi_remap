#!/usr/bin/env python
import os
import sys
import logging
import subprocess
from tempfile import NamedTemporaryFile

import numpy as np
from pymongo import MongoClient

sys.path.insert(0, '../lib/python')
from ncbi_remap.parser import parse_fastq_screen

DEBUG = False
logging.basicConfig(level=logging.INFO)

client = MongoClient(host='cn0330', port=27022)
db = client['sra2']
remap = db['remap']


def clean_types(dic):
    newDic = {}
    for k, v in dic.items():
        if isinstance(v, np.int64):
            newDic[k] = int(v)
        elif isinstance(v, np.float64):
            newDic[k] = float(v)
        elif v is np.nan:
            pass
        elif v is None:
            pass
        else:
            newDic[k] = v
    return newDic


def get_records():
    """Get SRX/SRR records.

    Get a list of SRX/SRR that have fastq_screen defined. Return a list of
    dicts.
    """
    return list(remap.aggregate([
        {'$unwind': '$runs'},
        {
            '$match': {
                'runs.pre_aln_workflow.fastq_screen': {'$exists': 1}
            }
        },
        {
            '$project': {
                'srx': '$_id',
                'srr': '$runs.srr',
                '_id': 0
            }
        },
    ]))


def parse_records(records):
    for record in records:
        srr = record['srr']
        fname = '/data/MiegNCBI/ncbi_remap/output/prealignment/raw/{srx}/{srr}/{srr}_1.fastq_screen.txt'.format(**record)

        if os.path.exists(fname):
            logging.info('Files exists for:  {srx}/{srr}'.format(**record))
            df = parse_fastq_screen(srr, fname)
            dd = {k[1]: v for k, v in df.to_dict('index').items()}

            if DEBUG:
                logging.info(dd)
            else:
                logging.info('Updating database.')
                remap.update_one(
                    {'runs.srr': srr},
                    {
                        '$set': {
                            'runs.$.pre_aln_workflow.fastq_screen': clean_types(dd)
                        }
                    }
                )


if __name__ == '__main__':
    # Pull out a list of srx/srr that have fastq screen
    records = get_records()
    logging.info('There were {:,} records'.format(len(records)))
    parse_records(records)
