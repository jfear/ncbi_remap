#!/usr/bin/env python
"""Populate fields of the Remap collection from the Ncbi collection.

This program initializes the Remap mongoDB collection to include fields from
the Ncbi collection. At this time I am only copying over identifiers. Other
metadata will eventually be copied over.
"""
import sys
import mongoengine as me
from sramongo.mongo_schema import Ncbi

sys.path.insert(0, '../lib/python')
from ncbi_remap.mongo_schema import Remap

client = me.connect(db='sra2', host='localhost', port=27022)

for ncbi in Ncbi.objects():
    sample = {
        'srs': ncbi.sra.sample.sample_id,
        'biosample': ncbi.sra.sample.BioSample,
        'gsm': ncbi.sra.sample.GEO,
    }

    runs = [{'srr': r.run_id} for r in ncbi.sra.run]

    contacts = []
    for s in  ncbi.biosample:
        for contact in s.contacts:
            contacts.append({
                'first_name': contact.first_name,
                'last_name': contact.last_name,
                'email': contact.email,
            })

    Remap.objects(srx=ncbi.srx).update_one(
        srx=ncbi.srx,
        sra=ncbi.sra.submission.submission_id,
        srp=ncbi.sra.study.study_id,
        bioproject=ncbi.sra.study.BioProject,
        sample=sample,
        add_to_set__runs=runs,
        add_to_set__contacts=contacts,
        add_to_set__papers=ncbi.pubmed,
        upsert=True
    )
