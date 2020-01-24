#!/usr/bin/env python
""" Set of helper scripts for file handling """
import re
from collections.abc import Mapping
from typing import Union
from yaml import full_load

import numpy as np
import pandas as pd
from IPython.display import display


def patterns(fname):
    with open(fname) as fh:
        return full_load(fh)


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
        self._store.put('prealn/queue', self.ids, data_columns=True,
                        format='table')

    def initialize_data_table(self, key, names):
        """Create a new empty datatable.

        Parameters
        ----------
        key : str
            Data store key.
        names : str or List[str]
            Name or names of the columns. If a single name is given then a
            pd.Series will be created, otherwise a pd.DataFrame.

        """
        if isinstance(names, (list, tuple)) and len(names) == 1:
            names = names[0]

        if isinstance(names, str):
            self._store[key] = pd.Series(name=names,
                                         index=self.ids.index)
        else:
            self._store[key] = pd.DataFrame(columns=names,
                                            index=self.ids.index)

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

    def _append_new_ids(self, key: str, new_srx_srr: pd.DataFrame,
                        fill=np.nan):
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


class LazyDict(Mapping):
    """LazyDictionary which runs a function when called.

    Based on answer here: https://stackoverflow.com/questions/16669367/setup-dictionary-lazily

    Examples
    --------
    >>> settings = LazyDict({'expensive1': (expensive_to_compute, 1), 'expensive2': (expensive_to_compute, 2)})

    """
    def __init__(self, *args, **kwargs):
        self._func_dict = dict(*args, **kwargs)
        self._dat_dict = dict()

    def __getitem__(self, key):
        if self._dat_dict.get(key):
            return self._dat_dict[key]

        func, arg = self._func_dict.__getitem__(key)
        self._dat_dict[key] = func(arg)
        return self._dat_dict[key]

    def __iter__(self):
        return iter(self._func_dict)

    def __len__(self):
        return len(self._func_dict)

    def keys(self):
        return self._func_dict.keys()


class GffRow(object):
    def __init__(self, row):
        (
            self.seqid,
            self.source,
            self.type,
            self.start,
            self.end,
            self.score,
            self.strand,
            self.phase,
            self.attributes,
        ) = row.strip().split("\t")
        self.is_gene = self.type == "gene"
        self.parsed_attributes = self.parse_attributes()

    def parse_attributes(self):
        parsed_attributes = {}
        for attr in self.attributes.strip(";").split(";"):
            mm = re.search(r"(?P<key>.*?)(\s+|=)(?P<value>.*)", attr.strip())
            if mm:
                parsed_attributes[mm.group("key").strip()] = mm.group("value").strip().strip('"')
        return parsed_attributes

    def __getitem__(self, key):
        return self.parsed_attributes[key]
