"""Construct the smaple section."""
import re

import numpy as np
import pandas as pd

from ncbi_remap.logging import logger
from ncbi_remap.sample_lists import get_complete_rnaseq
from ncbi_remap.mongo import mongo_connect


HEADER = ['title', 'organism', 'study', 'runs', 'GEO Experiment', 'GEO Sample', 'BioSample ID', 'BioProject',
          'pubmed', 'pubmed_title', 'pubmed_citation', 'pubmed_authors', 'contact',
          'sex', 'developmental stage', 'tissue', 'cell type', 'molecule', 'description', 'raw file', ]

OUTPUT = '../output/geo-wf/rnaseq_metadata.tsv'


def get_sample_info(samples):
    client, db, ncbi = mongo_connect()

    metadata = pd.DataFrame(list(ncbi.aggregate([
        {
            '$match': {
                '_id': {'$in': samples}
            }
        },
        {
            '$project': {
                '_id': 0,
                'sample_name': '$srx',
                'title': '$sra.sample.title',
                'organism': '$sra.sample.scientific_name',
                'study': '$sra.study.study_id',
                'runs': '$runs.srr',
                'GEO Experiment': '$sra.study.GEO',
                'GEO Sample': '$sra.sample.GEO',
                'BioSample ID': '$sra.sample.BioSample',
                'BioProject': '$sra.study.BioProject',
                'pubmed': '$pubmed.pubmed_id',
                'pubmed_title': '$pubmed.title',
                'pubmed_citation': '$pubmed.citation',
                'pubmed_authors': '$pubmed.authors',
                'contact': '$biosample.contacts',
                'molecule': '$sra.experiment.library_selection',
                'raw file': '$srx',
            }
        }

    ])))
    client.close()
    return metadata.set_index('sample_name')


def munge_contact(metadata):
    values = metadata['contact']

    res = []
    for v in values:
        if v in [np.nan, [], {}, [[]], [{}], [[{}]]]:
            res.append(np.nan)
            continue

        contact = v[0][0]
        if contact:
            first = contact.get('first_name', '')
            last = contact.get('last_name', '')
            email = contact.get('email', '')
            res.append(f'{first} {last} <{email}>')
            continue

    metadata['contact'] = res


def munge_authors(metadata):
    values = metadata['pubmed_authors']

    res = []
    for v in values:
        if v in [np.nan, [], {}, [[]], [{}], [[{}]]]:
            res.append(np.nan)
            continue

        authors = v[0]
        if authors:
            auts = []
            for author in authors:
                first = author.get('first_name', '')[0]
                last = author.get('last_name', '')
                auts.append(last + first)
            res.append('|'.join(auts))
            continue

    metadata['pubmed_authors'] = res


def munge_lists(metadata, col):
    values = metadata[col]

    res = []
    for v in values:
        if v in [np.nan, [], {}, [[]], [{}], [[{}]]]:
            res.append(np.nan)
            continue

        res.append('|'.join(v))

    metadata[col] = res


def clean_title(metadata):
    values = metadata['title']

    res = []
    for v in values:
        if re.search(r'currently private', str(v)):
            res.append('not provided')
        else:
            res.append(v)

    metadata['title'] = res


def main():
    samples = get_complete_rnaseq()
    metadata = get_sample_info(samples)
    munge_contact(metadata)
    munge_authors(metadata)
    munge_lists(metadata, 'runs')
    munge_lists(metadata, 'pubmed')
    munge_lists(metadata, 'pubmed_title')
    munge_lists(metadata, 'pubmed_citation')
    clean_title(metadata)

    # bag of words: keywords created using TfIDF
    bow = pd.read_parquet('../output/geo-wf/bow.parquet')
    bow.index.name = 'sample_name'
    bow.columns = ['description']

    # Hand cleaned characteristics
    characteristics = pd.read_csv('../output/geo-wf/rnaseq_characteristics.tsv', sep='\t', index_col=0)
    characteristics.index.name = 'sample_name'

    # TODO: Add QC metrics

    df = metadata.join(bow).join(characteristics)
    df[HEADER].to_csv(OUTPUT, sep='\t')


if __name__ == '__main__':
    main()
    logger.info('Script Complete')
