import os
import re
from io import StringIO
import pandas as pd
from collections import OrderedDict


def get_file_regex(filepattern):
    """Given a file name pattern, build a python regex.

    This was taken from snakemake.
    """
    _wildcard_regex = re.compile(
        r"""
        \{
            (?=(   # This lookahead assertion emulates an 'atomic group'
                   # which is required for performance
                \s*(?P<name>\w+)                    # wildcard name
                (\s*,\s*
                    (?P<constraint>                 # an optional constraint
                        ([^{}]+ | \{\d+(,\d+)?\})*  # allow curly braces to nest one level
                    )                               # ...  as in '{w,a{3,5}}'
                )?\s*
            ))\1
        \}
        """, re.VERBOSE)

    f = []
    last = 0
    wildcards = set()
    for match in _wildcard_regex.finditer(filepattern):
        f.append(re.escape(filepattern[last:match.start()]))
        wildcard = match.group("name")
        if wildcard in wildcards:
            if match.group("constraint"):
                raise ValueError(
                    "Constraint regex must be defined only in the first "
                    "occurence of the wildcard in a string.")
            f.append("(?P={})".format(wildcard))
        else:
            wildcards.add(wildcard)
            f.append("(?P<{}>{})".format(wildcard, match.group("constraint") if
                                         match.group("constraint") else ".+"))
        last = match.end()
    f.append(re.escape(filepattern[last:]))
    f.append("$")  # ensure that the match spans the whole file
    return "".join(f)


def get_files(files, pattern):
    """Given a list of files and a naming pattern returns a generator with list of files.

    Returns
    -------
    generator: (sample_id, file name)
    """

    file_pattern = get_file_regex(pattern)
    for file in files:
        m = re.match(file_pattern, file)
        if m:
            d = m.groupdict()
            yield (d['sample'], file)


def parse_files(files, pattern, parser):
    """Parses a set of files form some output.

    Given a list of files and file naming pattern, builds a list of files to
    parse. Then given a parsing function, parses the data into a pandas
    dataframe.
    """
    dfs = []
    for file in get_files(files, pattern):
        dfs.append(parser(*file))
    return pd.concat(dfs)


# Parsing functions

def parse_fqscreen(sample, file):
    """Parser for fastq screen.

    Adapted from multiqc.

    Returns
    -------
    pandas.DataFrame: A single row dataframe.
    """
    with open(file, 'r') as fh:
        parsed = OrderedDict()
        for l in fh:
            fqs = re.search(r"^(\S+)\s+(\d+)\s+(\d+)\s+([\d\.]+)\s+(\d+)\s+([\d\.]+)\s+(\d+)\s+([\d\.]+)\s+(\d+)\s+([\d\.]+)\s+(\d+)\s+([\d\.]+)$", l)
            if fqs:
                org = fqs.group(1)
                parsed[(org, 'reads_processed', 'count')] = int(fqs.group(2))
                parsed[(org, 'unmapped', 'count')] = int(fqs.group(3))
                parsed[(org, 'unmapped', 'percent')] = float(fqs.group(4))
                parsed[(org, 'one_hit_one_library', 'count')] = int(fqs.group(5))
                parsed[(org, 'one_hit_one_library', 'percent')] = float(fqs.group(6))
                parsed[(org, 'multiple_hits_one_library', 'count')] = int(fqs.group(7))
                parsed[(org, 'multiple_hits_one_library', 'percent')] = float(fqs.group(8))
                parsed[(org, 'one_hit_multiple_libraries', 'count')] = int(fqs.group(9))
                parsed[(org, 'one_hit_multiple_libraries', 'percent')] = float(fqs.group(10))
                parsed[(org, 'multiple_hits_multiple_libraries', 'count')] = int(fqs.group(11))
                parsed[(org, 'multiple_hits_multiple_libraries', 'percent')] = float(fqs.group(12))

        if len(parsed) == 0:
            return None
        else:
            return pd.DataFrame(parsed, index=[sample])


def parse_inferExperiment(sample, file):
    """Parse rseqc infer expeirment."""
    with open(file, 'r') as fh:
        parsed = OrderedDict()
        for l in fh:
            fqs = re.search(r"^(.+?):\s+([\d\.]+)$", l)
            if fqs:
                parsed[fqs.group(1)] = float(fqs.group(2))

        if len(parsed) == 0:
            return None
        else:
            return pd.DataFrame(parsed, index=[sample])


def parse_geneBodyCoverage(sample, file):
    """Parse rseqc genebody coverage."""
    with open(file, 'r') as fh:
        lines = fh.readlines()
        header = lines[0].strip().split('\t')[1:]
        values = lines[1].strip().split('\t')[1:]
        parsed = OrderedDict()
        for k, v in zip(header, values):
            parsed[int(k)] = float(v)
        if len(parsed) == 0:
            return None
        else:
            return pd.DataFrame(parsed, index=[sample])


def parse_bamStat(sample, file):
    """Parse rseqc bam stat."""
    with open(file, 'r') as fh:
        parsed = OrderedDict()
        for l in fh:
            fqs = re.search(r"^(.+?):\s*(\d+)$", l)
            if fqs:
                parsed[fqs.group(1)] = int(fqs.group(2))

        if len(parsed) == 0:
            return None
        else:
            return pd.DataFrame(parsed, index=[sample])


def parse_tin(sample, file):
    """Parse rseqc tin."""
    with open(file, 'r') as fh:
        lines = fh.readlines()
        header = lines[0].strip().split('\t')[1:]
        values = lines[1].strip().split('\t')[1:]
        parsed = OrderedDict()
        for k, v in zip(header, values):
            parsed[k] = float(v)
        if len(parsed) == 0:
            return None
        else:
            return pd.DataFrame(parsed, index=[sample])


def parse_picardCollect_summary(sample, file):
    """Parser for picard collectRNAMetrics summary."""
    with open(file, 'r') as fh:
        for l in fh:
            if l.startswith('#'):
                continue
            if l.startswith('PF_BASES'):
                parsed = l
                parsed += next(fh)
                break

        if len(parsed) == 0:
            return None
        else:
            df = pd.read_csv(StringIO(parsed), sep='\t')
            df.index = [sample]
            return df


def parse_picardCollect_hist(sample, file):
    """Parser for picard collectRNAMetrics summary."""
    with open(file, 'r') as fh:
        for l in fh:
            if l.startswith('#'):
                continue
            if l.startswith('normalized'):
                parsed = l
                while True:
                    try:
                        parsed += next(fh)
                    except StopIteration:
                        break
                break

        if len(parsed) == 0:
            return None
        else:
            df = pd.read_csv(StringIO(parsed), sep='\t', index_col=0).T
            df.index = [sample]
            return df


