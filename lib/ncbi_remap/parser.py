import re
from io import StringIO
from collections import OrderedDict

import numpy as np
import pandas as pd


# Parsing functions
def parse_fastq_summary(fname):
    """Parser for fastq summary table.

    Returns
    -------
    pandas.DataFrame: A single row dataframe.
    """
    dtypes = {
        'libsize_R1': np.float64,
        'avgLen_R1': np.float64,
        'libsize_R2': np.float64,
        'avgLen_R2': np.float64,
    }

    return pd.read_csv(fname, sep='\t', dtype=dtypes)


def parse_fastq_screen(fname):
    """Parser for fastq screen.

    Adapted from multiqc.

    Returns
    -------
    pandas.DataFrame: A single row dataframe.
    """
    with open(fname, 'r') as fh:
        header = ['reference', 'type', 'value']
        parsed = []

        for l in fh:
            fqs = re.search(r"^(\S+)\s+(\d+)\s+(\d+)\s+-?([\d\.]+)\s+(\d+)\s+([\d\.]+)\s+(\d+)\s+([\d\.]+)\s+(\d+)\s+([\d\.]+)\s+(\d+)\s+([\d\.]+)$", l)
            if fqs:
                org = fqs.group(1)
                parsed.append((org, 'reads_processed_count',
                               int(fqs.group(2))))
                parsed.append((org, 'unmapped_count',
                               int(fqs.group(3))))
                parsed.append((org, 'unmapped_percent',
                               float(fqs.group(4))))
                parsed.append((org, 'one_hit_one_library_count',
                               int(fqs.group(5))))
                parsed.append((org, 'one_hit_one_library_percent',
                               float(fqs.group(6))))
                parsed.append((org, 'multiple_hits_one_library_count',
                               int(fqs.group(7))))
                parsed.append((org, 'multiple_hits_one_library_percent',
                               float(fqs.group(8))))
                parsed.append((org, 'one_hit_multiple_libraries_count',
                               int(fqs.group(9))))
                parsed.append((org, 'one_hit_multiple_libraries_percent',
                               float(fqs.group(10))))
                parsed.append((org, 'multiple_hits_multiple_libraries_count',
                               int(fqs.group(11))))
                parsed.append((org, 'multiple_hits_multiple_libraries_percent',
                               float(fqs.group(12))))

        if len(parsed) == 0:
            return None
        else:
            df = pd.DataFrame(parsed, columns=header)
            udf = df.set_index(['reference', 'type']).unstack()
            udf.columns = udf.columns.droplevel()
            return udf.reset_index()


def parse_picardCollect_summary(fname):
    """Parser for picard collectRNAMetrics summary."""
    dtypes = {
        'CODING_BASES': np.int64,
        'CORRECT_STRAND_READS': np.int64,
        'IGNORED_READS': np.int64,
        'INCORRECT_STRAND_READS': np.int64,
        'INTERGENIC_BASES': np.int64,
        'INTRONIC_BASES': np.int64,
        'LIBRARY': np.float64,
        'MEDIAN_3PRIME_BIAS': np.float64,
        'MEDIAN_5PRIME_BIAS': np.float64,
        'MEDIAN_5PRIME_TO_3PRIME_BIAS': np.float64,
        'MEDIAN_CV_COVERAGE': np.float64,
        'NUM_R1_TRANSCRIPT_STRAND_READS': np.float64,
        'NUM_R2_TRANSCRIPT_STRAND_READS': np.float64,
        'NUM_UNEXPLAINED_READS': np.float64,
        'PCT_CODING_BASES': np.float64,
        'PCT_CORRECT_STRAND_READS': np.float64,
        'PCT_INTERGENIC_BASES': np.float64,
        'PCT_INTRONIC_BASES': np.float64,
        'PCT_MRNA_BASES': np.float64,
        'PCT_R1_TRANSCRIPT_STRAND_READS': np.float64,
        'PCT_R2_TRANSCRIPT_STRAND_READS': np.float64,
        'PCT_RIBOSOMAL_BASES': np.float64,
        'PCT_USABLE_BASES': np.float64,
        'PCT_UTR_BASES': np.float64,
        'PF_ALIGNED_BASES': np.int64,
        'PF_BASES': np.int64,
        'READ_GROUP': np.float64,
        'RIBOSOMAL_BASES': np.float64,
        'SAMPLE': np.float64,
        'UTR_BASES': np.int64,
    }
    with open(fname, 'r') as fh:
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
            return pd.read_csv(StringIO(parsed), sep='\t', na_values='?',
                               dtype=dtypes)


def parse_picardCollect_hist(fname):
    """Parser for picard collectRNAMetrics summary."""
    parsed = ''
    with open(fname, 'r') as fh:
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
            df.reset_index(drop=True, inplace=True)
            df.columns = [f'pos_{x}' for x in range(101)]
            return df


def parse_featureCounts_counts(fname):
    """Parser for subread feature counts."""
    df = pd.read_csv(fname, sep='\t', comment='#')
    df.columns = ['FBgn', 'chr', 'start', 'end', 'strand', 'length', 'count']
    df.set_index('FBgn', inplace=True)
    return df['count'].to_frame()


def parse_featureCounts_jcounts(fname):
    """Parser for subread feature jcounts."""
    header = [
        'PrimaryGene', 'SecondaryGenes', 'Site1_chr', 'Site1_location',
        'Site1_strand', 'Site2_chr', 'Site2_location', 'Site2_strand',
        'count',
    ]
    df = pd.read_csv(fname, sep='\t', header=None,
                     names=header, skiprows=1, low_memory=False)
    return df.drop('SecondaryGenes', axis=1)


def parse_featureCounts_summary(fname):
    """Parse rseqc bam stat."""
    with open(fname, 'r') as fh:
        parsed = OrderedDict()
        for l in fh:
            fqs = re.search(r"^(.+?)\s+(\d+)$", l)
            if fqs:
                parsed[fqs.group(1)] = int(fqs.group(2))
        if len(parsed) == 0:
            return None
        else:
            return pd.DataFrame(parsed, index=[0])


def parse_bamtools_stats(fname):
    """Parse bamtools stats."""
    with open(fname, 'r') as fh:
        parsed = OrderedDict()
        for l in fh:
            fqs = re.search(r"^(.+?):\s+(\d+).*$", l)
            if fqs:
                parsed[fqs.group(1).replace("'", "")] = int(fqs.group(2))
        if len(parsed) == 0:
            return None
        else:
            df = pd.DataFrame(parsed, index=[0])
            df['Percent Mapped'] = df['Mapped reads'] / df['Total reads'] * 100
            df['Percent Forward'] = df['Forward strand'] / df['Total reads'] * 100
            df['Percent Reverse'] = df['Reverse strand'] / df['Total reads'] * 100
            df['Percent Failed QC'] = df['Failed QC'] / df['Total reads'] * 100
            df['Percent Duplicates'] = df['Duplicates'] / df['Total reads'] * 100
            df['Percent Paired-end'] = df['Paired-end reads'] / df['Total reads'] * 100
            return df


def parse_picard_markduplicate_metrics(fname):
    """Parser for picard markduplicates."""
    with open(fname, 'r') as fh:
        for l in fh:
            if l.startswith('## METRICS CLASS'):
                dat = next(fh)
                dat += next(fh)
                break
    return pd.read_csv(StringIO(dat), sep='\t', comment='#', na_values='?')


def parse_samtools_idxstats(fname):
    """Parser for samtools idxstats."""
    df = pd.read_csv(fname, sep='\t', header=None)
    df.columns = ['chrom', 'length', '# mapped reads', '# unmapped reads']
    return df.set_index('chrom')


def parse_samtools_stats(fname):
    """Parse rseqc samtools stats."""
    with open(fname, 'r') as fh:
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
            return pd.DataFrame(parsed, index=[0])


def parse_atropos(fname):
    """Parse atropos."""
    with open(fname, 'r') as fh:
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

        if len(parsed) == 0:
            return None
        else:
            df = pd.DataFrame(parsed, index=[0])

            if [x for x in df.columns if 'Read 1' in x]:
                # PE
                df['pct_read1_adapters'] = df['Total read pairs processed_Read 1 with adapter'] / df['Total read pairs processed'] * 100

                df['pct_read2_adapters'] = df['Total read pairs processed_Read 2 with adapter'] / df['Total read pairs processed'] * 100
            else:
                # SE
                df['pct_read1_adapters'] = df['Reads with adapters'] / df['Total reads processed'] * 100

            return df, pd.concat(cnts.values())


def parse_hisat2(fname):
    """Parse hisat2."""
    with open(fname, 'r') as fh:
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
            return pd.DataFrame(parsed, index=[0])
