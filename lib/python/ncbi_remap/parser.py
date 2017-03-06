import os
import re
from io import StringIO
import pandas as pd
from collections import OrderedDict

from lcdblib.parse.fastqc import FastQC


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


def parse_files(files, pattern, parser, **kwargs):
    """Parses a set of files form some output.

    Given a list of files and file naming pattern, builds a list of files to
    parse. Then given a parsing function, parses the data into a pandas
    dataframe.
    """
    dfs = []
    for file in get_files(files, pattern):
        dfs.append(parser(*file, **kwargs))
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
    parsed = ''
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


def parse_dupradar(sample, file):
    """Parser for picard collectRNAMetrics summary."""
    df = pd.read_csv(file, sep='\t', index_col=0)
    df.columns = ['FBgn', 'geneLength', 'allCountsMulti', 'filteredCountsMulti', 'dupRateMulti',
                  'dupsPerIdMulti', 'RPKMulti', 'RPKMMulti', 'allCounts', 'filteredCounts',
                  'dupRate', 'dupsPerId', 'RPK', 'RPKM', 'mhRate']
    df['sample'] = sample
    return df.set_index(['sample', 'FBgn'])


def parse_featureCounts_counts(sample, file):
    """Parser for picard collectRNAMetrics summary."""
    df = pd.read_csv(file, sep='\t', comment='#')
    df.columns = ['FBgn', 'chr', 'start', 'end', 'strand', 'length', 'count']
    df['sample'] = sample
    return df.set_index(['sample', 'FBgn'])


def parse_featureCounts_summary(sample, file):
    """Parse rseqc bam stat."""
    with open(file, 'r') as fh:
        parsed = OrderedDict()
        for l in fh:
            fqs = re.search(r"^(.+?)\s+(\d+)$", l)
            if fqs:
                parsed[fqs.group(1)] = int(fqs.group(2))
        if len(parsed) == 0:
            return None
        else:
            return pd.DataFrame(parsed, index=[sample])


def parse_bamtools_stats(sample, file):
    """Parse bamtools stats."""
    with open(file, 'r') as fh:
        parsed = OrderedDict()
        for l in fh:
            fqs = re.search(r"^(.+?):\s+(\d+).*$", l)
            if fqs:
                parsed[fqs.group(1)] = int(fqs.group(2))
        if len(parsed) == 0:
            return None
        else:
            df = pd.DataFrame(parsed, index=[sample])
            df['Percent Mapped'] = df['Mapped reads'] / df['Total reads'] * 100
            df['Percent Forward'] = df['Forward strand'] / df['Total reads'] * 100
            df['Percent Reverse'] = df['Reverse strand'] / df['Total reads'] * 100
            df['Percent Failed QC'] = df['Failed QC'] / df['Total reads'] * 100
            df['Percent Duplicates'] = df['Duplicates'] / df['Total reads'] * 100
            df['Percent Paired-end'] = df['Paired-end reads'] / df['Total reads'] * 100
            return df


def parse_picard_markduplicate_metrics(sample, file):
    """Parser for picard markduplicates."""
    df = pd.read_csv(file, sep='\t', comment='#')
    df['sample'] = sample
    return df.set_index('sample')


def parse_samtools_idxstats(sample, file):
    """Parser for samtools idxstats."""
    df = pd.read_csv(file, sep='\t', header=None)
    df.columns = ['chrom', 'length', '# mapped reads', '# unmapped reads']
    df['sample'] = sample
    return df.set_index('sample')


def parse_samtools_stats(sample, file):
    """Parse rseqc samtools stats."""
    with open(file, 'r') as fh:
        parsed = OrderedDict()
        for l in fh:
            if l.startswith('SN'):
                fqs = re.search(r"^SN\s+(.+?):\s+([\d\.]+)\s.*$", l)
                if fqs:
                    if '.' in fqs.group(2):
                        parsed[fqs.group(1)] = float(fqs.group(2))
                    else:
                        parsed[fqs.group(1)] = int(fqs.group(2))

        if len(parsed) == 0:
            return None
        else:
            return pd.DataFrame(parsed, index=[sample])


def parse_md5(sample, file):
    """Parser for md5sum."""
    df = pd.read_csv(file, sep='\s+', header=None)
    df.columns = ['md5', 'file']
    df['sample'] = sample
    df.drop('file', axis=1, inplace=True)
    return df.set_index('sample')


def parse_libsize(sample, file):
    """Parser for md5sum."""
    with open(file, 'r') as fh:
        lsize = int(fh.read().strip())
        return pd.DataFrame({'libsize': lsize}, index=[sample])


def parse_atropos(sample, file):
    """Parse atropos."""
    with open(file, 'r') as fh:
        parsed = OrderedDict()
        block = None
        for l in fh:
            if l.startswith('==='):
                block = re.search(r"^=== (.+) ===$", l).group(1)
            else:
                if block == 'Summary':
                    l = l.replace(',', '')
                    fqs = re.search(r"^(.+?):\s+([\d\.]+)\s.*$", l)
                    if fqs:
                        parsed[fqs.group(1)] = int(fqs.group(2))
                else:
                    fqs = re.search(r"^.*Trimmed:\s+(\d+)\s.*$", l)
                    if fqs:
                        key = 'Number {} trimmed'.format(block)
                        parsed[key] = int(fqs.group(1))
        if len(parsed) == 0:
            return None
        else:
            return pd.DataFrame(parsed, index=[sample])


def parse_fastqc(sample, file, field=''):
    """Parse fastqc."""
    if field:
        return FastQC.parse_from_zip(sample, file)[field]
    else:
        return FastQC.parse_from_zip(sample, file)

def parse_fastqc_seq_quality(sample, file):
    """Parse fastqc."""
    fqc = parse_fastqc(sample, file, field='Per sequence quality scores')
    return fqc.df
