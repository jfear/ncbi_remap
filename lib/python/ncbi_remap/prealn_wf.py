"""Output lists of samples meeting different criteria."""
import pandas as pd

# Constants
LIBSIZE_CUTOFF = 1e5      # 100,000 reads
READLEN_CUTOFF = 30       # 30 bp
STRAND_CUTOFF1 = 0.75     # 75% stranded reads
STRAND_CUTOFF2 = 0.95     # 95% stranded reads
UNALIGN_CUTOFF = 0.50     # 50% Reads Unaligned
CONTAMINATION_CUTOFF = 40 # 50% reads mapping to DM6


def remove_rows(df, column, values):
    if values is None:
        return df

    return df[~df[column].isin(values)]


def keep_rows(df, column, values):
    if values is None:
        return df

    return df[df[column].isin(values)]


def srr_per_srx(store):
    """Output counts of SRRs per SRX.

    Parameters
    ----------
    store : pd.HDFStore
        HDF5 data store.
    cutoff : float
        Cutoff criteria to use.
    filter_srrs : list-like
        A list-like of SRRs to remove.

    Returns
    -------
    pd.DataFrames
        List of samples that meet criteria

    """
    return store['prealn/complete'].groupby('srx').count()


def libsize(store, cutoff=1e5, filter_srrs=None, keep_srrs=None):
    """Output lists of sample that meet libsize criteria.

    Parameters
    ----------
    store : pd.HDFStore
        HDF5 data store.
    cutoff : float
        Cutoff criteria to use.
    filter_srrs : list-like
        A list-like of SRRs to remove.

    Returns
    -------
    pd.DataFrames
        List of samples that meet criteria

    """
    df = remove_rows(store['prealn/workflow/fastq'], 'srr', filter_srrs)
    df = keep_rows(df, 'srr', keep_srrs)
    df['libsize'] = df[['libsize_R1', 'libsize_R2']].max(axis=1)

    return df.loc[df['libsize'] >= cutoff, ['srx', 'srr']]



def readlen(store, cutoff=30, filter_srrs=None, keep_srrs=None):
    """Output lists of sample that meet libsize criteria.

    Parameters
    ----------
    store : pd.HDFStore
        HDF5 data store.
    cutoff : float
        Cutoff criteria to use.
    filter_srrs : list-like
        A list-like of SRRs to remove.

    Returns
    -------
    pd.DataFrames
        List of samples that meet criteria

    """
    df = remove_rows(store['prealn/workflow/fastq'], 'srr', filter_srrs)
    df = keep_rows(df, 'srr', keep_srrs)
    df['len'] = df[['libsize_R1', 'libsize_R2']].max(axis=1)

    return df.loc[df['len'] >= cutoff, ['srx', 'srr']]


def strandedness(store, cutoff=.75, filter_srrs=None, keep_srrs=None):
    """Output lists of sample strandedness.

    Parameters
    ----------
    store : pd.HDFStore
        HDF5 data store.
    cutoff : float
        Cutoff criteria to use
    filter_srrs : list-like
        A list-like of SRRs to remove.
    keep_srrs : list-like
        A list-like of SRRs to keep.

    Returns
    -------
    tuple of pd.DataFrames
        Tuple with (first_stranded, second_stranded, unstranded).

    """
    first = remove_rows(store['prealn/workflow/collectrnaseqmetrics/first'], 'srr', filter_srrs)
    second = remove_rows(store['prealn/workflow/collectrnaseqmetrics/second'], 'srr', filter_srrs)

    first = keep_rows(first, 'srr', filter_srrs)
    second = keep_rows(second, 'srr', filter_srrs)

    f = first.loc[first['PCT_CORRECT_STRAND_READS'] >= cutoff, ['srx', 'srr']]
    s = second.loc[second['PCT_CORRECT_STRAND_READS'] >= cutoff, ['srx', 'srr']]
    u = first[~(first.srr.isin(f.srr) | first.srr.isin(s.srr))]
    return f, s, u


def mappability(store, cutoff=30, filter_srrs=None, keep_srrs=None):
    """Output lists of sample that meet mappability criteria.

    Parameters
    ----------
    store : pd.HDFStore
        HDF5 data store.
    cutoff : float
        Cutoff criteria to use.
    filter_srrs : list-like
        A list-like of SRRs to remove.

    Returns
    -------
    pd.DataFrames
        List of samples that meet criteria

    """

    se = store['prealn/workflow/hisat2/SE'][['srx', 'srr', 'num_reads', 'num_uniquely_aligned']]
    se['prop_unique_aligned'] = se['num_uniquely_aligned'] / se['num_reads']

    pe = store['prealn/workflow/hisat2/PE'][['srx', 'srr', 'num_reads', 'num_concordant_reads_uniquely_aligned']]
    pe['prop_unique_aligned'] = pe['num_concordant_reads_uniquely_aligned'] / pe['num_reads']

    df = pd.concat([se[['srx', 'srr', 'prop_unique_aligned']], pe[['srx', 'srr', 'prop_unique_aligned']]])

    df = remove_rows(df, 'srr', filter_srrs)
    df = keep_rows(df, 'srr', keep_srrs)

    return df.loc[df['prop_unique_aligned'] >= cutoff, ['srx', 'srr']]


def contamination(store, cutoff=50, filter_srrs=None, keep_srrs=None):
    """Output lists of sample that meet contamination criteria.

    Parameters
    ----------
    store : pd.HDFStore
        HDF5 data store.
    cutoff : float
        Cutoff criteria to use.
    filter_srrs : list-like
        A list-like of SRRs to remove.

    Returns
    -------
    pd.DataFrames
        List of samples that meet criteria

    """

    df = store['prealn/workflow/fastq_screen']
    df = df[['srx', 'srr', 'reference', 'one_hit_one_library_percent']].set_index(['srx', 'srr', 'reference']).unstack()
    df.columns = df.columns.droplevel(0)
    df.reset_index(inplace=True)

    df = remove_rows(df, 'srr', filter_srrs)
    df = keep_rows(df, 'srr', keep_srrs)

    return df.loc[df['dm6'] >= cutoff, ['srx', 'srr']]
