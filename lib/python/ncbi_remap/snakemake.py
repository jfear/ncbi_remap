#!/usr/bin/env python
"""Set of helpers for use with snakemake."""
import os
import numpy as np


def wrapper_for(path):
    URI = '../lcdb-wf/wrappers/wrappers'
    return 'file:' + os.path.join(URI, path)


def put_flag(fname, flag):
    with open(fname, 'w') as fh:
        fh.write(flag)


def get_flag(fname):
    with open(fname) as fh:
        return fh.read().strip()


def clean_types(dic):
    if isinstance(dic, dict):
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
    elif isinstance(dic, list):
        newLis = []
        for row in dic:
            newLis.append(clean_types(row))
        return newLis


