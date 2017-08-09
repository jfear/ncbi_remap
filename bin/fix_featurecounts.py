#!/usr/bin/env python
import os
import sys
import logging
import subprocess
from tempfile import NamedTemporaryFile

import numpy as np
from pymongo import MongoClient

sys.path.insert(0, '../lib/python')
from ncbi_remap.parser import parse_featureCounts_jcounts, parse_featureCounts_summary

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

    Get a list of SRX/SRR that have featurecounts defined. Return a list of
    dicts.
    """
    return list(remap.aggregate([
        {'$unwind': '$runs'},
        {
            '$match': {
                'runs.pre_aln_workflow.featurecounts': {'$exists': 1}
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
        fnameS = '/data/MiegNCBI/ncbi_remap/output/prealignment/raw/{srx}/{srr}/{srr}.hisat2.bam.feature_counts.counts.summary'.format(**record)
        fnameJ = '/data/MiegNCBI/ncbi_remap/output/prealignment/raw/{srx}/{srr}/{srr}.hisat2.bam.feature_counts.counts.jcounts'.format(**record)

        if os.path.exists(fnameS) & os.path.exists(fnameJ):
            logging.info('Files exists for:  {srx}/{srr}'.format(**record))
            dfJ = parse_featureCounts_jcounts(srr, fnameJ)

            if dfJ.shape[0] > 0:
                num_junction_reads = dfJ.loc[~dfJ.PrimaryGene.isnull(), 'count'].sum()
            else:
                num_junction_reads = 0

            dfS = parse_featureCounts_summary(srr, fnameS)
            dd = dfS.to_dict('index')[srr]
            dd['Assigned_Junction'] = num_junction_reads

            if DEBUG:
                logging.info(dd)
            else:
                logging.info('Updating database.')
                remap.update_one(
                    {'runs.srr': srr},
                    {
                        '$set': {
                            'runs.$.pre_aln_workflow.featurecounts': clean_types(dd)
                        }
                    }
                )


if __name__ == '__main__':
    # Pull out a list of srx/srr that have feature counts
    records = get_records()
    logging.info('There were {:,} records'.format(len(records)))
    parse_records(records)
