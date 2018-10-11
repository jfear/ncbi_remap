"""Helpers for NLP."""
from copy import copy
from yaml import load
from pathlib import Path

from nltk.util import ngrams

from .io import LazyDict


def load_yaml_lookup(fname):
    with open(fname) as fh:
        dat = load(fh)

    # extend dictionary mappings to remove '-' from key names.
    dat_extend = copy(dat)
    for k, v in dat.items():
        dat_extend[k.replace('-', '')] = v
        dat_extend[k.replace('-', ' ')] = v

    return dat_extend


def load_text(fname):
    with open(fname) as fh:
        return fh.read().strip().split('\n')


cwd = Path(__file__).parent.resolve()

lookups = LazyDict({
    'library_strategy': (load_yaml_lookup, Path(cwd, '../lookups/library_strategy.yaml')),
    'tissue': (load_yaml_lookup, Path(cwd, '../lookups/tissue.yaml')),
    'sex': (load_yaml_lookup, Path(cwd, '../lookups/sex.yaml')),
    'dev_stage': (load_yaml_lookup, Path(cwd, '../lookups/dev_stage.yaml')),
    'cell_type': (load_yaml_lookup, Path(cwd, '../lookups/cell_type.yaml')),
    'stopwords': (load_text, Path(cwd, '../lookups/stopwords.txt')),
})
