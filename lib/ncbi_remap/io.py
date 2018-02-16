#!/usr/bin/env python
""" Set of helper scripts for file handling """
from typing import Union

import pandas as pd
from IPython.display import display


class HDFStore(object):
    """A facade interface to pandas.HDFStore"""
    def __init__(self, path, mode='a', srx_srr=None, *args, **kwargs):
        self._store = pd.HDFStore(path, mode, *args, **kwargs)
        if self._store.__contains__('ids'):
            return

        if srx_srr is not None:
            self.initialize_store(srx_srr)
        else:
            print(
                'To initialize the data store run '
                'HDFStore.initialize_store or set `srx_srr` when '
                'opening the data store.'
            )

    def __getitem__(self, key: str):
        return self._store[key]

    def __setitem__(self, key: str, value: Union[pd.Series, pd.DataFrame]):
        self._store[key] = value

    def __delitem__(self, key: str):
        del self._store[key]

    def contains(self, key):
        return self._store.__contains__(key)

    def close(self):
        self._store.close()

    @property
    def root(self):
        return self._store.root

    def info(self):
        return self._store.info()

    @property
    def ids(self):
        return self._store['ids']

    def initialize_store(self, srx_srr: pd.DataFrame):
        """Run this once to create basic data store.

        Parameters
        ----------
        srx_srr : pd.DataFrame
            A dataframe with an index with SRX, SRR.

        """
        self._store['ids'] = srx_srr
        self._store['prealn/queue'] = self.ids

    def initialize_data_table(self, key, names, fill_bool=False):
        """Create a new empty datatable.

        Parameters
        ----------
        key : str
            Data store key.
        names : str or List[str]
            Name or names of the columns. If a single name is given then a
            pd.Series will be created, otherwise a pd.DataFrame.
        fill_bool : bool
            If true use the boolean value to fill missing instead of
            default 'NA'.

        """
        if fill_bool:
            fill = False
        else:
            fill = 'NA'

        if isinstance(names, (list, tuple)) and len(names) == 1:
            names = names[0]

        if isinstance(names, str):
            self._store[key] = pd.Series(name=names,
                                         index=self.ids.index).fillna(fill)
        else:
            self._store[key] = pd.DataFrame(columns=names,
                                            index=self.ids.index).fillna(fill)

    def update_ids_from_db(self, current_srx_srr):
        """

        Parameters
        ----------
        current_srx_srr : pd.DataFrame
            An updated set SRX, SRR ides. Note if an id is not in this list it
            will be removed from the data store.

        """
        orig_srx_srr = self.ids

        removed_ids = orig_srx_srr[
            ~orig_srx_srr.index.isin(current_srx_srr.index)
        ]

        new_ids = current_srx_srr[
            ~current_srx_srr.index.isin(orig_srx_srr.index)
        ]

        self._store['ids'] = current_srx_srr
        for key in self._store.keys():
            if key == '/ids':
                continue

            if 'complete' in key:
                continue

            self._append_new_ids(key, new_ids)
            self._drop_removed_ids(key, removed_ids)

    def _append_new_ids(self, key: str, new_srx_srr: pd.DataFrame, fill='NA'):
        """Append new sra ids."""
        if new_srx_srr.shape[0] == 0:
            return

        if not self._store.__contains__(key):
            self._store[key] = new_srx_srr
            return

        orig = self._store[key]
        new_dat = new_srx_srr.copy()

        try:
            for col in orig.columns:
                new_dat[col] = fill
        except AttributeError:
            if orig.dtype.name == 'bool':
                fill = False

            new_dat[orig.name] = fill
            new_dat = new_dat.loc[:, orig.name]

        new = pd.concat([orig, new_dat])
        self._store[key] = new[~new.index.duplicated('first')]

    def _drop_removed_ids(self, key: str, removed_srx_srr: pd.DataFrame):
        """Remove sra ids that are no longer present.

        Somtime after running a database update SRA removes ids. I want remove
        them from my work too.

        """

        if removed_srx_srr.shape[0] == 0:
            return

        orig = self._store[key]
        curr_ids = orig[~orig.index.isin(removed_srx_srr.index)]
        self._store[key] = curr_ids

    def _tuple_to_frame(self, srx_srr: tuple):
            return pd.DataFrame(index=pd.MultiIndex.from_tuples(
                [srx_srr], names=['srx', 'srr']))

    def pop(self, key: str,
                  srx_srr: Union[tuple, pd.DataFrame]):

        if isinstance(srx_srr, tuple):
            srx_srr = self._tuple_to_frame(srx_srr)

        self._drop_removed_ids(key, srx_srr)

    def push(self, key: str,
                   srx_srr: Union[tuple, pd.DataFrame]):

        if isinstance(srx_srr, tuple):
            srx_srr = self._tuple_to_frame(srx_srr)

        self._append_new_ids(key, srx_srr)

    def get_srxs(self, key):
        return self._store[key].index.get_level_values('srx').unique().tolist()

    def get_srrs(self, key):
        return self._store[key].index.get_level_values('srr').unique().tolist()


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
    if store.__contains__(key) & (force is True):
        # If the store exists delete
        curr = data
    elif store.__contains__(key):
        curr = store[key].copy()
        curr = curr.append(data, ignore_index=True).drop_duplicates()

    store[key] = curr


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

    Builds a query with the current kwargs, if the query matches then the
    record is removed.

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
    df = store[key]
    subset = df[~df.srr.isin(srrs)].copy()
    store[key] = subset


def add_data_columns(store, key, **kwargs):
    dat = store[key].copy()
    add_table(store, key, data=dat, force=True, **kwargs)


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
