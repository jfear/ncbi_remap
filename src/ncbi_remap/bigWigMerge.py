#!/usr/bin/env python
"""Script to merge and normalize bigwigs."""
import argparse
from argparse import RawDescriptionHelpFormatter as Raw
from logging import DEBUG, INFO
import multiprocessing as mp
from collections import OrderedDict
from itertools import groupby
from functools import partial

import numpy as np
import pyBigWig

from .logging import logger


def parse_region(string):
    string = string.replace(',', '').replace(':', '|').replace('_', '|').replace('-', '|')
    chrom, start, end = string.split('|')
    return OrderedDict([('chrom', chrom), ('start', int(start)), ('end', int(end))])


def arguments():
    """Pull in command line arguments."""
    DESCRIPTION = """\
    """

    parser = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=Raw)

    parser.add_argument("--input", dest="input", nargs='*', action='store', required=True,
                        help="A list of input BigWigs.")

    parser.add_argument("--output", dest="output", action='store', required=True,
                        help="Output BigWig")

    parser.add_argument("--chromSizes", dest="chromSizes", action='store', required=True,
                        help="chromSizes")

    parser.add_argument("--region", dest="region", action='store', required=False,
                        help="Region to focus on.")

    parser.add_argument("--threads", dest="nthreads", type=int, action='store', default=1, required=False,
                        help="Number of threads to use.")

    parser.add_argument("--agg", dest="func", type=str, action='store', default='sum', required=False,
                        choices=['sum', 'mean', 'median', 'max', '75', '95'],
                        help="Type of aggregatition to use [sum default].")

    parser.add_argument("--scale", dest="scale", action='store_true', required=False,
                        help="Scale reads by the number of files input.")

    parser.add_argument("--threshold", dest="threshold", type=int, action='store', required=False,
                        help="Set positons below this value to 0.")

    parser.add_argument("--debug", dest="debug", action='store_true', required=False,
                        help="Turn on debug output.")

    args = parser.parse_args()

    if args.region:
        args.region = parse_region(args.region)

    if args.debug:
        logger.setLevel(DEBUG)
    else:
        logger.setLevel(INFO)

    return args


def parse_chromSizes(fname):
    res = []
    with open(fname) as fh:
        for line in fh.readlines():
            chrom, size = line.strip().split('\t')
            res.append((chrom, int(size)))
    return res


def parse_bigWig(fname, chrom, start, end):
    with pyBigWig.open(fname) as bw:
        values = bw.values(chrom, start, end, numpy=True)
    return values


class Sum(object):
    def __init__(self, args, start, end):
        self.start = start
        self.end = end
        self.cnt = 0
        self.values = np.zeros(end - start, dtype=np.float64)
        self.lock = mp.Lock()
        self.args = args

    def add(self, array):
        self.lock.acquire()
        self.cnt += 1
        self.values += array
        self.lock.release()

    def scale_values(self):
            self.values /= self.cnt

    def theshold_values(self):
            self.values[self.values < self.threshold] = 0

    def adjust(self):
        if self.args.scale:
            self.scale_values()

        if self.args.threshold:
            self.threshold_values()


class BedGraphEntry(object):
    def __init__(self):
        self.chroms = []
        self.starts = []
        self.ends = []
        self.values = []

    def add_entry(self, chrom, start, end, values):
        _pos = start
        for k, v in groupby(values):
            self.chroms.append(chrom)
            self.values.append(k)
            self.starts.append(_pos)
            _pos += len(list(v))
            self.ends.append(_pos)

    def clean_dtypes(self):
        self.values = [float(val) for val in self.values]

    def __iter__(self):
        return (x for x in [self.chroms, self.starts, self.ends, self.values])


def create_entry_w_callback(args, chrom, start, end):
    pool = mp.Pool(processes=args.nthreads)
    sumArray = Sum(args, start, end)
    for input in args.input:
        pool.apply_async(parse_bigWig, (input, chrom, start, end), callback=sumArray.add)
    pool.close()
    pool.join()

    sumArray.adjust()
    return sumArray.values


def unpacking_apply_along_axis(arguments):
    func1d, axis, arr = arguments
    return np.apply_along_axis(func1d, axis, arr)


def aggregate(args, dat):
    if args.func == 'sum':
        func = np.sum
    elif args.func == 'mean':
        func = np.mean
    elif args.func == 'median':
        func = np.median
    elif args.func == 'max':
        func = np.max
    elif args.func == '75':
        func = partial(np.percentile, q=75)
    elif args.func == '95':
        func = partial(np.percentile, q=95)

    return func(dat, axis=1)


def create_entry(args, chrom, start, end, chunksize=1e6):
    _start = start
    chunks = []
    while True:
        _end = int(min(end, _start + chunksize))
        logger.debug(f'Coords: {start}-{_end}')
        pool = mp.Pool(processes=args.nthreads)
        res = np.column_stack(
            pool.map_async(partial(parse_bigWig, chrom=chrom, start=_start, end=_end), args.input).get(999999)
        )
        pool.close()
        pool.join()

        chunks.append(aggregate(args, res))

        if _end == end:
            break
        else:
            _start = _end

    return np.concatenate(chunks)


def main():
    args = arguments()
    logger.debug('input: %s', args.input)
    logger.debug('region: %s', args.region)
    logger.debug('output: %s', args.output)
    logger.debug('threads: %d', args.nthreads)

    # grab the chromosomes sizes for use as the BigWig header
    chromSizes = parse_chromSizes(args.chromSizes)
    bedGraph = BedGraphEntry()
    if args.region:
        chrom = args.region['chrom']
        start = args.region['start']
        end = args.region['end']

        logger.info(f'Merging {chrom}: {start:,}-{end:,}')
        entry = create_entry(args, chrom, start, end)
        bedGraph.add_entry(chrom, start, end, entry)
    else:
        for chromSize in chromSizes:
            chrom = chromSize[0]
            start = 0
            end = chromSize[1]
            logger.info(f'Merging {chrom}: {start:,}-{end:,}')
            entry = create_entry(args, chrom, start, end)
            bedGraph.add_entry(chrom, start, end, entry)

    bedGraph.clean_dtypes()
    with pyBigWig.open(args.output, 'w') as out:
        out.addHeader(chromSizes)
        out.addEntries(*bedGraph)

    logger.info('Script Complete')
