#!/usr/bin/env python
"""Set of helpers for use with snakemake."""


def wrapper_for(path):
    return 'file:' + path


def put_flag(fname, flag):
    with open(fname, 'w') as fh:
        fh.write(flag)


def get_flag(fname):
    with open(fname) as fh:
        return fh.read().strip()

