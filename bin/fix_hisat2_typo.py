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
                'runs.pre_aln_workflow.hisat2.num_uniquely_algined': {'$exists': 1}
            }
        },
        {
            '$project': {
                'srx': '$_id',
                'srr': '$runs.srr',
                'curr': '$runs.pre_aln_workflow.hisat2.num_uniquely_algined',
                '_id': 0
            }
        },
    ]))


def parse_records(records):
    for record in records:
        logging.info('{srx}/{srr}'.format(**record))
        srr = record['srr']
        remap.update_one(
            {'runs.srr': srr},
            {
                '$set': {
                    'runs.$.pre_aln_workflow.hisat2.num_uniquely_aligned': record['curr']
                },
            }
        )

        remap.update_one(
            {'runs.srr': srr},
            {
                '$unset': {
                    'runs.$.pre_aln_workflow.hisat2.num_uniquely_aligned': ""
                },
            }
        )

if __name__ == '__main__':
    # Pull out a list of srx/srr that have fastq screen
    records = get_records()
    logging.info('There were {:,} records'.format(len(records)))
    parse_records(records)
