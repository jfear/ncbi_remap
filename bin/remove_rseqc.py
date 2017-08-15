#!/usr/bin/env python
"""Helper script to remove rseqc infer_experiment and bam_stat from the database."""
import os
import sys
import logging

from pymongo import MongoClient

client = MongoClient(host='cn0330', port=27022)
db = client['sra2']
remap = db['remap']

while True:
    ud = remap.update_many({'runs.pre_aln_workflow.infer_expeirment': {'$exists': 1}}, {'$unset': {'runs.$.pre_aln_workflow.infer_expeirment': ""}})
    print(ud.modified_count)
    if ud.modified_count == 0:
        break

while True:
    ud2 = remap.update_many({'runs.pre_aln_workflow.bam_stat': {'$exists': 1}}, {'$unset': {'runs.$.pre_aln_workflow.bam_stat': ""}})
    print(ud2.modified_count)
    if ud2.modified_count == 0:
        break
