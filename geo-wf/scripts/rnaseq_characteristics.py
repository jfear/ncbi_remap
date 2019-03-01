"""Create a cleaner list of sample characteristics including:
    * sex
    * tissue
    * cell type
    * developmental stage
"""
import re
from yaml import load

import numpy as np
import pandas as pd

from ncbi_remap.sample_lists import get_complete_rnaseq
from ncbi_remap.mongo import mongo_connect

OUTPUT = '../output/geo-wf/rnaseq_characteristics.tsv'


def normalize_columns(attrs, replace_dict):
    regex = re.compile('|'.join(replace_dict.keys()))

    attrs['name'] = attrs['name'].str.lower().apply(lambda s: regex.sub(lambda x: replace_dict[x.group()], s))
    attrs['value'] = attrs['value'].str.lower().apply(lambda s: regex.sub(lambda x: replace_dict[x.group()], s))
    return attrs


def get_attributes(samples):
    client, db, ncbi = mongo_connect()

    attrs = pd.DataFrame(list(ncbi.aggregate([
        {
            '$match': {
                '_id': {'$in': samples}
            }
        },
        {
            '$unwind': {
                'path': '$sra.sample.attributes'
            }
        },
        {
            '$project': {
                '_id': 0,
                'srx': '$srx',
                'title': '$sra.sample.title',
                'name': '$sra.sample.attributes.name',
                'value': '$sra.sample.attributes.value'
            }
        },
    ])))

    client.close()

    replace_dict = {
        '_': ' ',
        '; ': ';',
        ' / ': '/',
        ' - ': '-',
        "''": '',
    }

    return normalize_columns(attrs, replace_dict)


def fix_sex(attrs):
    attrsp = attrs.pivot_table(values='value', index='srx', columns='name', aggfunc='first')
    dat = attrsp['sex'].copy()
    titles = attrs[['srx', 'title']].set_index('srx').title.drop_duplicates()

    # Fix based on hand mappings
    with open('config/sex.yaml') as fh:
        norm = load(fh)

    for k, v in norm.items():
        if v == 'None':
            norm[k] = np.nan

    dat.replace(norm, inplace=True)

    # If we are missing sex information try to pull it out from other attributes.
    bow = {}
    for srx, row in attrsp.iterrows():
        try:
            title = titles[srx].lower()
        except (KeyError, AttributeError):
            title = ''
        bow[srx] = ' '.join(row.dropna().values.tolist() + [title])

    missing = dat[dat.isna()].index.tolist()
    for srx in missing:
        try:
            string = bow[srx]
        except KeyError:
            continue

        mf = re.search(r'\bfemale\b', string)
        mm = re.search(r'\bmale\b', string)
        if (mf is not None) & (mm is not None):
            dat[srx] = 'mixed'
        elif mf:
            dat[srx] = 'female'
        elif mm:
            dat[srx] = 'male'

    return dat


def fix_dev(attrs):
    # First I want to combine multiple columns that are related to dev stage.
    _attrs = attrs.copy()
    replace_dict = {
        'dev stage': 'developmental stage',
        'dev-stage': 'developmental stage',
        'developemntal stage': 'developmental stage',
        'development stage': 'developmental stage',
        'developmental stage': 'developmental stage',
        'develpomental stage': 'developmental stage',
    }

    _attrs = normalize_columns(_attrs, replace_dict)

    attrsp = _attrs.pivot_table(values='value', index='srx', columns='name', aggfunc='first')
    dat = attrsp['developmental stage'].copy()

    # Fix based on hand mappings
    with open('config/dev_stage.yaml') as fh:
        norm = load(fh)

    for k, v in norm.items():
        if v == 'None':
            norm[k] = np.nan

    dat.replace(norm, inplace=True)
    # dev stage is too complicated to try to use other columns for info

    return dat


def fix_tissue(attrs):
    titles = attrs[['srx', 'title']].set_index('srx').title.drop_duplicates()
    # First I want to combine multiple columns that are related to tissue
    _attrs = attrs.copy()
    replace_dict = {
        'cell/tissue type': 'tissue',
        'tissue lib': 'tissue',
        'tissue source': 'tissue',
        'tissue type': 'tissue',
        'tissue/cell type': 'tissue',
    }

    _attrs = normalize_columns(_attrs, replace_dict)

    attrsp = _attrs.pivot_table(values='value', index='srx', columns='name', aggfunc='first')
    dat = attrsp['tissue'].copy()

    # Fix based on hand mappings
    with open('config/tissue.yaml') as fh:
        norm = load(fh)

    for k, v in norm.items():
        if v == 'None':
            norm[k] = np.nan

    dat.replace(norm, inplace=True)
    # If we are missing tissue information try to pull it out from other attributes.
    bow = {}
    for srx, row in attrsp.iterrows():
        try:
            title = titles[srx].lower()
        except (KeyError, AttributeError):
            title = ''

        bow[srx] = ' '.join(row.dropna().values.tolist() + [title])

    missing = dat[dat.isna()].index.tolist()
    for srx in missing:
        try:
            string = bow[srx]
        except KeyError:
            continue

        if re.search(r'\bwhole body\b', string):
            dat[srx] = 'whole body'
        elif re.search(r'\bhead\b', string):
            dat[srx] = 'head'
        elif re.search(r'\bovary\b', string):
            dat[srx] = 'ovary'
        elif re.search(r'\btestis\b', string):
            dat[srx] = 'testis'
        elif re.search(r'\bmushroom body\b', string):
            dat[srx] = 'mushroom body'

    return dat


def fix_cell(attrs):
    # First I want to combine multiple columns that are related to cell type
    _attrs = attrs.copy()
    replace_dict = {
        'cell/tissue type': 'cell type',
        'tissue/cell type': 'cell type',
        'cell class': 'cell type',
        'cell line': 'cell type',
        'cell line background': 'cell type',
        'cell subtype': 'cell type',
        'cell type': 'cell type',
        'cells': 'cell type',
        'cells derived from': 'cell type',
    }

    _attrs = normalize_columns(_attrs, replace_dict)

    attrsp = _attrs.pivot_table(values='value', index='srx', columns='name', aggfunc='first')
    dat = attrsp['cell type'].copy()

    # Fix based on hand mappings
    with open('config/cell_type.yaml') as fh:
        norm = load(fh)

    for k, v in norm.items():
        if v == 'None':
            norm[k] = np.nan

    dat.replace(norm, inplace=True)
    # cell type is too complicated to try to use other columns for info

    return dat


def main():
    samples = get_complete_rnaseq()
    sample_attributes = get_attributes(samples)
    sex = fix_sex(sample_attributes)
    dev_stage = fix_dev(sample_attributes)
    tissue = fix_tissue(sample_attributes)
    cell_type = fix_cell(sample_attributes)
    df = pd.concat([sex, dev_stage, tissue, cell_type], axis=1)
    df.to_csv(OUTPUT, sep='\t')


if __name__ == '__main__':
    main()
