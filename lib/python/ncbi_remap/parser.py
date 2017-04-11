import os
import re
from io import StringIO
import numpy as np
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
        if os.path.exists(file[1]):
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
    """Parser for dupradar."""
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
        subBlock = None
        for l in fh:
            l = l.replace(',', '')

            if l.startswith('==='):
                block = re.search(r"^=== (.+) ===$", l).group(1)
            else:
                # The summary block is a little different between SE and PE reads
                if block == 'Summary':
                    fqs = re.search(r"^(\w+[\s\w\(\)-]+?):\s+([\d\.]+)\s.*$", l)
                    fqs2 = re.search(r"^\s+([\w\s-]+?):\s+([\d\.]+)\s.*$", l)

                    if fqs:
                        subBlock = fqs.group(1)
                        parsed[subBlock] = int(fqs.group(2))
                    elif fqs2:
                        read = '_'.join([subBlock, fqs2.group(1)])
                        parsed[read] = int(fqs2.group(2))
                else:
                    fqs = re.search(r"^.*Trimmed:\s+(\d+)\s.*$", l)
                    if fqs:
                        key = 'Number {} trimmed'.format(block)
                        parsed[key] = int(fqs.group(1))
        if len(parsed) == 0:
            return None
        else:
            df = pd.DataFrame(parsed, index=[sample])

            if [x for x in df.columns if 'Read 1' in x]:
                # PE
                df['pct_read1_adapters'] = df['Total read pairs processed_Read 1 with adapter'] / df['Total read pairs processed'] * 100

                df['pct_read2_adapters'] = df['Total read pairs processed_Read 2 with adapter'] / df['Total read pairs processed'] * 100
            else:
                # SE
                df['pct_read1_adapters'] = df['Reads with adapters'] / df['Total reads processed'] * 100

            return df


def parse_hisat2(sample, file):
    """Parse hisat2."""
    with open(file, 'r') as fh:
        parsed = OrderedDict()
        header = ['# Reads', '# Unpaired', '# Unaligned', '# Uniquely Aligned', '# Multimappers', 'Pct Alignment']
        for i, l in enumerate(fh):
            fqs = re.search(r"^\s*([\d\.]+?)[%\s].*$", l)
            if fqs:
                if i == 5:
                    parsed[header[i]] = float(fqs.group(1))
                else:
                    parsed[header[i]] = int(fqs.group(1))

        if len(parsed) == 0:
            return None
        else:
            return pd.DataFrame(parsed, index=[sample])


def split_ranges(df):
    """Split ranges into bases.

    Fastqc sometimes collapses bases into ranges, this splits them back out.
    """
    rows = []
    for i, row in df.iterrows():
        try:
            if '-' in i:
                start, end = [int(x) for x in i.split('-')]
                for j in range(start, end + 1):
                    curr_row = row.copy()
                    curr_row.name = j
                    rows.append(curr_row)
            else:
                row.name = int(i)
                rows.append(row)
        except TypeError:
            rows.append(row)

    df = pd.concat(rows, axis=1).T
    df.index.name = 'base'
    return df


def parse_fastqc(sample, file, field=''):
    """Parse fastqc zip file.

    Takes a zip file and makes a FastQC object.

    Parameters
    ----------
    sample: str
        Sample name which will be added as row index.
    file: str
        Path to the fastqc zip file.
    field: str
        Name of specific Fastqc section to return. Look at
        lcdblib.parse.fastqc.FastQC.keys() for a list of possible names.

    Returns
    -------
    lcdblib.parse.fastqc.FastQC or lcdblib.parse.fastqc.FastQCBlock if `field`
    is provided.

    """
    if field:
        return FastQC.parse_from_zip(sample, file)[field]
    else:
        return FastQC.parse_from_zip(sample, file)


def parse_fastqc_per_seq_quality(sample, file):
    """Parse fastqc."""
    fqc = parse_fastqc(sample, file, field='Per sequence quality scores')
    df = fqc.df
    wide = df.T
    wide['sample'] = sample
    return wide.set_index('sample')


def parse_fastqc_per_base_seq_quality(sample, file):
    """Parse fastqc base quality"""
    fqc = parse_fastqc(sample, file, field='Per base sequence quality')
    df = fqc.df['Mean'].copy().to_frame()
    splitRanges = split_ranges(df)
    splitRanges['sample'] = sample
    return splitRanges.set_index(append=True, keys='sample').swaplevel()


def parse_fastqc_adapter_content(sample, file):
    """Parse fastqc."""
    fqc = parse_fastqc(sample, file, field='Adapter Content')
    df = fqc.df
    splitRanges = split_ranges(df)
    splitRanges['sample'] = sample
    return splitRanges.set_index(append=True, keys='sample').swaplevel()


def parse_fastqc_per_base_seq_content(sample, file):
    """Parse fastqc."""
    fqc = parse_fastqc(sample, file, field='Per base sequence content')
    df = fqc.df
    splitRanges = split_ranges(df)
    splitRanges['sample'] = sample
    return splitRanges.set_index(append=True, keys='sample').swaplevel()


def parse_fastqc_sequence_length(sample, file):
    """Parse fastqc."""
    fqc = parse_fastqc(sample, file, field='Sequence Length Distribution')
    df = fqc.df
    df['sample'] = sample
    return df.set_index(append=True, keys='sample').swaplevel()


def parse_fastqc_overrepresented_seq(sample, file):
    """Parse fastqc."""
    fqc = parse_fastqc(sample, file, field='Overrepresented sequences')
    df = fqc.df
    df['sample'] = sample
    return df.set_index(append=True, keys='sample').swaplevel()


def parse_fastqc_basic_stats(sample, file):
    """Parse fastqc."""
    fqc = parse_fastqc(sample, file, field='Basic Statistics')
    df = fqc.df.T
    df['sample'] = sample
    return df.set_index('sample')


def parse_fastqc_kmer_content(sample, file):
    """Parse fastqc."""
    fqc = parse_fastqc(sample, file, field='Kmer Content')
    df = fqc.df
    df.reset_index(inplace=True)
    df.set_index('Max Obs/Exp Position', inplace=True)
    splitRanges = split_ranges(df)
    splitRanges.index.name = 'Max Obs/Exp Position'
    splitRanges.reset_index(inplace=True)
    splitRanges['sample'] = sample
    return splitRanges.sort_values(['Sequence', 'Max Obs/Exp Position']).set_index(['sample', 'Sequence'])


def parse_fastqc_per_base_n_content(sample, file):
    """Parse fastqc."""
    fqc = parse_fastqc(sample, file, field='Per base N content')
    df = fqc.df
    splitRanges = split_ranges(df)
    splitRanges['sample'] = sample
    return splitRanges.set_index(append=True, keys='sample').swaplevel()


def parse_fastqc_per_seq_gc_content(sample, file):
    """Parse fastqc."""
    fqc = parse_fastqc(sample, file, field='Per sequence GC content')
    df = fqc.df
    df['sample'] = sample
    return df.set_index(append=True, keys='sample').swaplevel()


def parse_fastqc_seq_dup_level(sample, file):
    """Parse fastqc."""
    fqc = parse_fastqc(sample, file, field='Sequence Duplication Levels')
    df = fqc.df
    df['sample'] = sample
    return df.set_index(append=True, keys='sample').swaplevel()