#!/usr/bin/env python
""" Set of helper scripts for file handling """

import os
import numpy as np
import pandas as pd
from IPython.display import display


def build_index(store, key, columns=None):
    """Index a HDF5 table.

    Builds an index for an HDF5 table.

    Parameters
    ----------
    store : pd.io.pytables.HDFStore
        The data store to save to.
    key : str
        The path in the HDF5 store to save data to.
    columns : str or list, optional
        If 'all' then it will index using all data columns in the store. If a
        list of column headers is given then it will use those columns for
        indexing.

    """
    if columns is None:
        columns = []
    elif columns == 'all':
        columns = store[key].columns.tolist()

    store.create_table_index(key, columns=columns, optlevel=9, kind='full')


def add_table(store, key, data=None, force=None, **kwargs):
    """Create a new HDF5 table.

    Adds a dataframe to an HDF5 store and creates an index.

    Parameters
    ----------
    store : pd.io.pytables.HDFStore
        The data store to save to.
    key : str
        The path in the HDF5 store to save data to.
    data : pd.DataFrame
        The data to store.
    force : bool
        If True then delete the previous store if it exists.

    """
    defaults = {'columns': 'all'}
    defaults.update(kwargs)

    # If the store exists delete
    if store.__contains__(key) & (force is True):
            del store[key]
    elif store.__contains__(key):
        # Drop if duplicates
        if 'srr' in data.columns:
            data = data[~data.srr.isin(store[key].srr)].copy()
        elif 'srx' in data.columns:
            data = data[~data.srx.isin(store[key].srx)].copy()
        else:
            data = data[~data.isin(store[key].to_dict('list')).all(axis=1)].copy()

    if data.shape[0] > 0:
        store.append(key, data, data_columns=True, index=False)
        build_index(store, key, **defaults)


def add_id(store, key, **kwargs):
    """Adds an ID to the ids data store.

    Takes **kwargs and builds a table. Then tries to add the table with ncbi_remap.io.add_table.

    Parameters
    ----------
    store : pd.io.pytables.HDFStore
        The data store to save to.
    key : str
        The path in the HDF5 store to save data to.

    """
    table = pd.Series(kwargs).to_frame().T

    if store.__contains__(key):
        cols = store[key].columns
        table = table[cols]

    add_table(store, key, table)


def remove_id(store, key, **kwargs):
    """Removes an ID to the ids data store.

    Builds a query with the current kwargs, if the query matches then the record is removed.

    Parameters
    ----------
    store : pd.io.pytables.HDFStore
        The data store to save to.
    key : str
        The path in the HDF5 store to save data to.

    """
    query = ['{} == {}'.format(k, v) for k, v in kwargs.items()]
    store.remove(key, ' & '.join(query))


def remove_chunk(store, key, srrs, **kwargs):
    """Removes an ID to the ids data store.

    If the SRR is not in the current collection, then append the srx and srr.

    Parameters
    ----------
    store : pd.io.pytables.HDFStore
        The data store to save to.
    key : str
        The path in the HDF5 store to save data to.
    srrs : list
        A list of SRRs to remove.

    """
    defaults = {'columns': 'all'}
    defaults.update(kwargs)

    df = store[key]
    subset = df[~df.srr.isin(srrs)].copy()
    store.put(key, subset, format='table')
    build_index(store, key, **defaults)


def check_alignment(store, pattern, **kwargs):
    """Checks for ALIGNMENT_BAD file.

    If there is an ALIGNMENT_BAD file then remove from the queue, add to
    complete, and add to 'prealn/alignment_bad'.

    Parameters
    ----------
    store : pd.io.pytables.HDFStore
        The data store to save to.
    patter : str
        File naming pattern for the ALIGNEMNT_BAD file.
    **kwargs
        Keywords needed to fill the pattern.

    """
    ab = pattern.format(**kwargs)
    if os.path.exists(ab):
        remove_id(store, 'prealn/queue', **kwargs)
        add_id(store, 'prealn/alignment_bad', **kwargs)
        add_id(store, 'prealn/complete', **kwargs)
        return True


def check_layout(store, pattern, **kwargs):
    """Checks LAYOUT file.

    Parses layout file and adds to the corresponding hdf5 lists.

        * 'layout/SE'
        * 'layout/PE'
        * 'layout/keep_R1'
        * 'layout/keep_R2'

    Parameters
    ----------
    store : pd.io.pytables.HDFStore
        The data store to save to.
    patter : str
        File naming pattern for the ALIGNEMNT_BAD file.
    **kwargs
        Keywords needed to fill the pattern.

    """
    with open(pattern.format(**kwargs)) as fh:
        strand = fh.read().strip()
        if strand == 'SE':
            key = 'layout/SE'
        elif strand == 'PE':
            key = 'layout/PE'
        elif strand == 'keep_R1':
            key = 'layout/keep_R1'
        elif strand == 'keep_R2':
            key = 'layout/keep_R2'
        add_id(store, key, **kwargs)


def check_strand(store, pattern, **kwargs):
    """Checks STRAND file.

    Parses strand file and adds to the corresponding hdf5 lists.

        * 'strand/first'
        * 'strand/second'
        * 'strand/unstranded'

    Parameters
    ----------
    store : pd.io.pytables.HDFStore
        The data store to save to.
    patter : str
        File naming pattern for the ALIGNEMNT_BAD file.
    **kwargs
        Keywords needed to fill the pattern.

    """
    with open(pattern.format(**kwargs)) as fh:
        strand = fh.read().strip()
        if (strand == 'first_strand') | (strand == 'same_strand'):
            key = 'strand/first'
        elif (strand == 'second_strand') | (strand == 'opposite_strand'):
            key = 'strand/second'
        elif strand == 'unstranded':
            key = 'strand/unstranded'
        add_id(store, key, **kwargs)


class remapDesign(object):
    def __init__(self, sampleTable):
        """ Import sample table and make a set of sampleID lists for easy
        reference.

        Parameters
        ----------
        sampleTable: str
            Filename of sample table "design file"

        Attributes
        ----------
        sampleTable: pandas.DataFrame
            The imported sampleTable
        exp_unit: dict
            Dictionary for each experimental unit with a list of corresponding
            samples.

        """
        header = [
                  "Run in SRA (or group of runs from a single library)",
                  "Runs from sublibraries (one library ran in multiple lanes)",
                  "From SRP project",
                  "Project description and aim",
                  "reference (pm stands for PMID)",
                  "Authors of the study (not always precise enough in SRA)",
                  "Biosample ID (SAM*)",
                  "Summary of sample (information collected from SRA by Magic",
                  "RNA aim: total, polyA, nascent, small RNA, 3', 5'..(manual",
                  "Experiment detail (manual)",
                  "Targeted RNA, RIP CLIP, ribosome profiling, gene selection",
                  "Details on targeted experiment (manual)",
                  "Genomic DNA",
                  "Biological class, control (manual)",
                  "Control",
                  "Biological: treatment, disease, phenotype, gene modified ",
                  "Treatment details",
                  "Gene down regulated by RNAi or mutation (~hypomorph/amorph",
                  "Details on gene down",
                  "Gene overexpressed or ectopically expressed (~neomorph)",
                  "Details on ectopic or anomalous gene",
                  "Marker gene",
                  "Driver gene or promoter",
                  "Sex (NLM)",
                  "Sex (D+J)",
                  "Sex (Zhenxia) ",
                  "Sex (Brian)",
                  "NLM stage of development",
                  "Developmental stage (D+J)",
                  "Details on stage or age (and some experiment summary) (D+J",
                  "Development stage (Zhenxia)",
                  "Development stage (Brian)",
                  "Age stage (Brian)",
                  "NLM cell type and anatomy",
                  "Tissue (Zhenxia)",
                  "Tissue (Brian)",
                  "Tissue type (D+J)",
                  "Tissue detail (D+J)",
                  "Tissue description (D+J)",
                  "System type (if not simple tissue) (D+J)",
                  "System (D+J)",
                  "More system description (D+J)",
                  "NLM strain, treatment or note",
                  "Genotype (Zhenxia)",
                  "Background genotype (Brian)",
                  "Treatment perturbagen (Brian)",
                  "Treatment conditions (Brian)",
                  "Sample focus (Brian)",
                  "Cell type (Zhenxia)",
                  "Cell_type (Brian)",
                  "Cell line (D+J)",
                  "cell line details (D+J)",
                  "NLM cell line",
                  "Sample type (Zhenxia)",
                  "Sample type (Brian)",
                  "Tissue position (Brian)",
                  "Notes and flags (Brian)",
                  "Problem in SRA consistency or download",
                  "This SRR is a doublon of another (need to remove from SRA)",
                 ]

        # Import sample table
        self.sampleTable = pd.read_csv(sampleTable, sep='\t', quotechar='"',  na_values='NULL', comment='#', header=None, names=header, encoding='ISO-8859-1')

        # Make sure lane and replicate are strings
#         self.sampleTable.replicate = self.sampleTable.replicate.map(lambda x: str(x))
#         self.sampleTable.lane = self.sampleTable.lane.map(lambda x: str(x))

        # Create map of sample_ID to a more useful sample_name
#         self.sampleID2sampleName = self.sampleTable[['sample_id', 'sample_name']].set_index(\
#                 'sample_id').to_dict('dict')['sample_name']

        # Make useful sample lists
        ## target
#         self.dam = self.sampleTable[(self.sampleTable.target == 'dam')].sample_id.unique().tolist()
#         self.polII = self.sampleTable[(self.sampleTable.target == 'polII')].sample_id.unique().tolist()


        ## Create a dict for each experimental unit with all sample_ids
#         self.exp_unit = {}
#         for idx, grp in self.sampleTable.groupby(['sex', 'tissue', 'driver', 'target']):
#             self.exp_unit['_'.join(idx)] = sorted(list(set(grp.sample_id.values.tolist())), key=self.sortFn)

    def __repr__(self):
#         print('{0} rows, {1} columns'.format(*self.sampleTable.shape))
#         qgrid.QGridWidget(df=self.sampleTable.set_index(['sample_id', 'run'])).export()
        display(self.sampleTable)
        return ''

    def getSamples(self, categories):
        """ Given a list of categories return the intersected list.

        This is a helper method that allows you to refine a list of samples
        based on different criteria. It is a similar idea to the exp_unit, but
        allows you to directly create arbitrary lists by providing different
        categories.

        Example
        -------
        >>> sTable = remapDesign('../../../config/damid_sample_info.csv')
        >>> samples = sTable.getSamples(['male', 'testis', 'c587' , 'polII'])
        >>> assert samples == sTable.exp_unit['male_testis_c587_polII']

        """
        sets = [set(self.__dict__[x]) for x in categories]
        return sorted(list(sets[0].intersection(*sets[1:])))

    @staticmethod
    def sortFn(x):
        """ Sort function for sorting DamID sample_ids """
        if x.startswith('DamID'):
            return ('A', int(x.split('_')[1]))
        else:
            m = x.split('_')[0]
            return (m[0], int(m[1:])+100)
