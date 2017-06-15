#!/usr/bin/env python

from pymongo import MongoClient
with open('/home/fearjm/Projects/ncbi_remap/output/.mongodb_host', 'r') as fh:
    host = fh.read().strip()

client = MongoClient(host=host, port=27022)
db = client['sra2']
remap = db['remap']


def remove_prealn(srr):
    remap.find_one_and_update({'runs.srr': srr}, {
        '$set': {
            'runs.$.pre_aln_flags': [],
            'runs.$.pre_aln_workflow': {},
            'runs.$.libsize': {},
            'runs.$.avgReadLen': {},
            'runs.$.md5': {},
        }})
