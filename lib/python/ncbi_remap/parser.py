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
    generator: (srx, srr, file name)
    """

    file_pattern = get_file_regex(pattern)
    for file in files:
        m = re.match(file_pattern, file)
        if m:
            d = m.groupdict()
            yield (d['srx'], d['srr'], file)


def parse_files(files, pattern, parser, **kwargs):
    """Parses a set of files form some output.

    Given a list of files and file naming pattern, builds a list of files to
    parse. Then given a parsing function, parses the data into a pandas
    dataframe.
    """
    dfs = []
    for file in get_files(files, pattern):
        if os.path.exists(file[1]):
            try:
                dfs.append(parser(*file, **kwargs))
            except KeyError:
                pass
            except Exception as e:
                print(file[1])
                raise e

    return pd.concat(dfs)


# Parsing functions
def parse_fastq_summary(srx, srr, file):
    """Parser for fastq summary table.

    Returns
    -------
    pandas.DataFrame: A single row dataframe.
    """
    df = pd.read_csv(file, sep='\t')
    df.index = [srx, srr]
    df.index.names = ['srx', 'srr']
    return df


def parse_fastq_screen(srx, srr, file):
    """Parser for fastq screen.

    Adapted from multiqc.

    Returns
    -------
    pandas.DataFrame: A single row dataframe.
    """
    with open(file, 'r') as fh:
        header = ['srx', 'srr', 'reference', 'type', 'value']
        parsed = []

        for l in fh:
            fqs = re.search(r"^(\S+)\s+(\d+)\s+(\d+)\s+-?([\d\.]+)\s+(\d+)\s+([\d\.]+)\s+(\d+)\s+([\d\.]+)\s+(\d+)\s+([\d\.]+)\s+(\d+)\s+([\d\.]+)$", l)
            if fqs:
                org = fqs.group(1)
                parsed.append((srx, srr, org, 'reads_processed_count', int(fqs.group(2))))
                parsed.append((srx, srr, org, 'unmapped_count', int(fqs.group(3))))
                parsed.append((srx, srr, org, 'unmapped_percent', float(fqs.group(4))))
                parsed.append((srx, srr, org, 'one_hit_one_library_count', int(fqs.group(5))))
                parsed.append((srx, srr, org, 'one_hit_one_library_percent', float(fqs.group(6))))
                parsed.append((srx, srr, org, 'multiple_hits_one_library_count', int(fqs.group(7))))
                parsed.append((srx, srr, org, 'multiple_hits_one_library_percent', float(fqs.group(8))))
                parsed.append((srx, srr, org, 'one_hit_multiple_libraries_count', int(fqs.group(9))))
                parsed.append((srx, srr, org, 'one_hit_multiple_libraries_percent', float(fqs.group(10))))
                parsed.append((srx, srr, org, 'multiple_hits_multiple_libraries_count', int(fqs.group(11))))
                parsed.append((srx, srr, org, 'multiple_hits_multiple_libraries_percent', float(fqs.group(12))))

        if len(parsed) == 0:
            return None
        else:
            df = pd.DataFrame(parsed, columns=header)
            udf = df.set_index(['srx', 'srr', 'reference', 'type']).unstack()
            udf.columns = udf.columns.droplevel()
            return udf.reset_index().set_index(['srx', 'srr'])


def parse_picardCollect_summary(srx, srr, file):
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
            df.index = [srx, srr]
            df.index.names = ['srx', 'srr']
            df.replace('?', np.nan, inplace=True)
            return df


def parse_picardCollect_hist(srx, srr, file):
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
            df.index = [srx, srr]
            df.index.names = ['srx', 'srr']
            return df


def parse_featureCounts_counts(srx, srr, file):
    """Parser for subread feature counts."""
    df = pd.read_csv(file, sep='\t', comment='#')
    df.columns = ['FBgn', 'chr', 'start', 'end', 'strand', 'length', 'count']
    df['srx'] = srx
    df['srr'] = srr
    df.set_index(['srx', 'srr'], inplace=True)
    return df[['FBgn', 'count']]


def parse_featureCounts_jcounts(srx, srr, file):
    """Parser for subread feature jcounts."""
    header = ['PrimaryGene', 'SecondaryGenes', 'Site1_chr', 'Site1_location', 'Site1_strand', 'Site2_chr', 'Site2_location', 'Site2_strand', 'count']
    df = pd.read_csv(file, sep='\t', header=None, names=header, skiprows=1)
    df['srx'] = srx
    df['srr'] = srr
    df.set_index(['srx', 'srr'], inplace=True)
    return df


def parse_featureCounts_summary(srx, srr, file):
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
            df = pd.DataFrame(parsed, index=[srx, srr])
            df.index.names = ['srx', 'srr']
            return df


def parse_bamtools_stats(srx, srr, file):
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
            df = pd.DataFrame(parsed, index=[srx, srr])
            df.index.names = ['srx', 'srr']
            df['Percent Mapped'] = df['Mapped reads'] / df['Total reads'] * 100
            df['Percent Forward'] = df['Forward strand'] / df['Total reads'] * 100
            df['Percent Reverse'] = df['Reverse strand'] / df['Total reads'] * 100
            df['Percent Failed QC'] = df['Failed QC'] / df['Total reads'] * 100
            df['Percent Duplicates'] = df['Duplicates'] / df['Total reads'] * 100
            df['Percent Paired-end'] = df['Paired-end reads'] / df['Total reads'] * 100
            return df


def parse_picard_markduplicate_metrics(srx, srr, file):
    """Parser for picard markduplicates."""
    with open(file, 'r') as fh:
        for l in fh:
            if l.startswith('## METRICS CLASS'):
                dat = next(fh)
                dat += next(fh)
                break
    df = pd.read_csv(StringIO(dat), sep='\t', comment='#')
    df['srx'] = srx
    df['srr'] = srr
    df.set_index(['srx', 'srr'], inplace=True)
    return df


def parse_samtools_idxstats(srx, srr, file):
    """Parser for samtools idxstats."""
    df = pd.read_csv(file, sep='\t', header=None)
    df.columns = ['chrom', 'length', '# mapped reads', '# unmapped reads']
    df['srx'] = srx
    df['srr'] = srr
    df.set_index(['srx', 'srr'], inplace=True)
    return df


def parse_samtools_stats(srx, srr, file):
    """Parse rseqc samtools stats."""
    with open(file, 'r') as fh:
        parsed = OrderedDict()
        for l in fh:
            if l.startswith('SN'):
                fqs = re.search(r"^SN\s+(.+?):\s+([\d\.]+)\s.*$", l)
                if fqs:
                    name = fqs.group(1).replace(' ', '_')
                    value = fqs.group(2)

                    if ('.' in value) or ('average' in name):
                        parsed[name] = float(value)
                    else:
                        parsed[name] = int(value)

        if len(parsed) == 0:
            return None
        else:
            df = pd.DataFrame(parsed, index=[srx, srr])
            df.index.names = ['srx', 'srr']
            return df


def parse_md5(srx, srr, file):
    """Parser for md5sum."""
    df = pd.read_csv(file, sep='\s+', header=None)
    df.columns = ['md5', 'file']
    df['srx'] = srx
    df['srr'] = srr
    df.set_index(['srx', 'srr'], inplace=True)
    df.drop('file', axis=1, inplace=True)
    return df


def parse_libsize(srx, srr, file):
    """Parser for md5sum."""
    with open(file, 'r') as fh:
        lsize = int(fh.read().strip())
        df = pd.DataFrame({'libsize': lsize}, index=[srx, srr])
        df.index.names = ['srx', 'srr']
        return df


def parse_atropos(srx, srr, file):
    """Parse atropos."""
    with open(file, 'r') as fh:
        parsed = OrderedDict()
        block = None
        subBlock = None
        cnts = {}
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
                    elif l.startswith('length'):
                        # This will pull out the length count tables and make dataframes
                        cnts[block] = l
                        try:
                            while True:
                                l = next(fh)
                                if l.startswith('\n'):
                                    break
                                cnts[block] += l
                        except StopIteration:
                            pass
                        cnts[block] = pd.read_table(StringIO(cnts[block]))
                        cnts[block]['adapter'] = block
                        cnts[block]['srx'] = srx
                        cnts[block]['srr'] = srr
                        cnts[block].set_index(['srx', 'srr', 'adapter', 'length'], inplace=True)

        if len(parsed) == 0:
            return None
        else:
            df = pd.DataFrame(parsed, index=[srx, srr])
            df.index.names = ['srx', 'srr']

            if [x for x in df.columns if 'Read 1' in x]:
                # PE
                df['pct_read1_adapters'] = df['Total read pairs processed_Read 1 with adapter'] / df['Total read pairs processed'] * 100

                df['pct_read2_adapters'] = df['Total read pairs processed_Read 2 with adapter'] / df['Total read pairs processed'] * 100
            else:
                # SE
                df['pct_read1_adapters'] = df['Reads with adapters'] / df['Total reads processed'] * 100

            return df, pd.concat(cnts.values())


def parse_hisat2(srx, srr, file):
    """Parse hisat2."""
    with open(file, 'r') as fh:
        parsed = OrderedDict([
            ('num_reads', np.nan),
            ('num_reads_paired', np.nan),
            ('num_reads_unpaired', np.nan),
            ('num_concordant_reads_unaligned', np.nan),
            ('num_concordant_reads_uniquely_aligned', np.nan),
            ('num_concordant_multimappers', np.nan),
            ('num_discordant_reads_aligned', np.nan),
            ('num_unaligned', np.nan),
            ('num_uniquely_aligned', np.nan),
            ('num_multimappers', np.nan),
            ('per_alignment', np.nan)
        ])

        header = {
            'reads; of these:': 'num_reads',
            'were paired; of these:': 'num_reads_paired',
            'were unpaired; of these:': 'num_reads_unpaired',
            'aligned concordantly 0 times': 'num_concordant_reads_unaligned',
            'aligned concordantly exactly 1 time': 'num_concordant_reads_uniquely_aligned',
            'aligned concordantly >1 times': 'num_concordant_multimappers',
            'aligned discordantly 1 time': 'num_discordant_reads_aligned',
            'aligned 0 times': 'num_unaligned',
            'aligned exactly 1 time': 'num_uniquely_aligned',
            'aligned >1 times': 'num_multimappers',
            'overall alignment rate': 'per_alignment',
        }

        for l in fh:
            row = l.strip()
            if row.startswith('Warning') or row.startswith('---'):
                continue

            fqs = re.search(r"^([\d\.]+)[%]*\s([\(\)\d\.%]*\s)?(.*)$", row)
            if fqs:
                name = fqs.group(3)
                value = fqs.group(1)
                if name in header:
                    head = header[name]
                    parsed[head] = float(value)

        if len(parsed) == 0:
            return None
        else:
            df = pd.DataFrame(parsed, index=[0])
            df['srx'] = srx
            df['srr'] = srr
            return df.set_index(['srx', 'srr'])


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
    df.index.names = 'base'
    return df


def parse_fastqc(srx, srr, file, field=''):
    """Parse fastqc zip file.

    Takes a zip file and makes a FastQC object.

    Parameters
    ----------
    srx: str
        SRX name which will be added as row index.
    srr: str
        SRR name which will be added as row index.
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
        return FastQC.parse_from_zip(srx, srr, file)[field]
    else:
        return FastQC.parse_from_zip(srx, srr, file)


def parse_fastqc_per_seq_quality(srx, srr, file):
    """Parse fastqc."""
    fqc = parse_fastqc(srx, srr, file, field='Per sequence quality scores')
    df = fqc.df
    wide = df.T
    wide['srx'] = srx
    wide['srr'] = srr
    return wide.set_index(['srx', 'srr'])


def parse_fastqc_per_base_seq_quality(srx, srr, file):
    """Parse fastqc base quality"""
    fqc = parse_fastqc(srx, srr, file, field='Per base sequence quality')
    df = fqc.df['Mean'].copy().to_frame()
    splitRanges = split_ranges(df)
    splitRanges['srx'] = srx
    splitRanges['srr'] = srr
    return splitRanges.set_index(append=True, keys=['srx', 'srr']).swaplevel()


def parse_fastqc_adapter_content(srx, srr, file):
    """Parse fastqc."""
    fqc = parse_fastqc(srx, srr, file, field='Adapter Content')
    df = fqc.df
    splitRanges = split_ranges(df)
    splitRanges['srx'] = srx
    splitRanges['srr'] = srr
    return splitRanges.set_index(append=True, keys=['srx', 'srr']).swaplevel()


def parse_fastqc_per_base_seq_content(srx, srr, file):
    """Parse fastqc."""
    fqc = parse_fastqc(srx, srr, file, field='Per base sequence content')
    df = fqc.df
    splitRanges['srx'] = srx
    splitRanges['srr'] = srr
    return splitRanges.set_index(append=True, keys=['srx', 'srr']).swaplevel()


def parse_fastqc_sequence_length(srx, srr, file):
    """Parse fastqc."""
    fqc = parse_fastqc(srx, srr, file, field='Sequence Length Distribution')
    df = fqc.df
    df['srx'] = srx
    df['srr'] = srr
    return df.set_index(append=True, keys=['srx', 'srr']).swaplevel()


def parse_fastqc_overrepresented_seq(srx, srr, file):
    """Parse fastqc."""
    fqc = parse_fastqc(srx, srr, file, field='Overrepresented sequences')
    df = fqc.df
    df['srx'] = srx
    df['srr'] = srr
    return df.set_index(append=True, keys=['srx', 'srr']).swaplevel()


def parse_fastqc_basic_stats(srx, srr, file):
    """Parse fastqc."""
    fqc = parse_fastqc(srx, srr, file, field='Basic Statistics')
    df = fqc.df.T
    df['srx'] = srx
    df['srr'] = srr
    return df.set_index(['srx', 'srr'])


def parse_fastqc_kmer_content(srx, srr, file):
    """Parse fastqc."""
    fqc = parse_fastqc(srx, srr, file, field='Kmer Content')
    df = fqc.df
    df.reset_index(inplace=True)
    df.set_index('Max Obs/Exp Position', inplace=True)
    splitRanges = split_ranges(df)
    splitRanges.index.names = 'Max Obs/Exp Position'
    splitRanges.reset_index(inplace=True)
    splitRanges['srx'] = srx
    splitRanges['srr'] = srr
    return splitRanges.sort_values(['Sequence', 'Max Obs/Exp Position']).set_index(['srx', 'srr', 'Sequence'])


def parse_fastqc_per_base_n_content(srx, srr, file):
    """Parse fastqc."""
    fqc = parse_fastqc(srx, srr, file, field='Per base N content')
    df = fqc.df
    splitRanges = split_ranges(df)
    splitRanges['srx'] = srx
    splitRanges['srr'] = srr
    return splitRanges.set_index(append=True, keys=['srx', 'srr']).swaplevel()


def parse_fastqc_per_seq_gc_content(srx, srr, file):
    """Parse fastqc."""
    fqc = parse_fastqc(srx, srr, file, field='Per sequence GC content')
    df = fqc.df
    df['srx'] = srx
    df['srr'] = srr
    return df.set_index(append=True, keys=['srx', 'srr']).swaplevel()


def parse_fastqc_seq_dup_level(srx, srr, file):
    """Parse fastqc."""
    fqc = parse_fastqc(srx, srr, file, field='Sequence Duplication Levels')
    df = fqc.df
    df['srx'] = srx
    df['srr'] = srr
    return df.set_index(append=True, keys=['srx', 'srr']).swaplevel()
