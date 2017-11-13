"""Output lists of samples meeting different criteria."""
import numpy as np
import pandas as pd
import scipy
from .logging import logger

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
    df = store['prealn/workflow/fastq'].copy()
    df.reset_index(inplace=True)
    df = remove_rows(df, 'srr', filter_srrs)
    df = keep_rows(df, 'srr', keep_srrs)
    df['libsize'] = df[['libsize_R1', 'libsize_R2']].max(axis=1)

    return df.loc[df['libsize'] >= cutoff, ['srx', 'srr']]


def libsize_cnts(store, filter_srrs=None, keep_srrs=None):
    """Output lists of sample that meet libsize criteria.

    Parameters
    ----------
    store : pd.HDFStore
        HDF5 data store.
    filter_srrs : list-like
        A list-like of SRRs to remove.

    Returns
    -------
    tuple of floats
       (minimum, median, max) of libsize.

    """
    df = store['prealn/workflow/fastq'].copy()
    df.reset_index(inplace=True)
    df = remove_rows(df, 'srr', filter_srrs)
    df = keep_rows(df, 'srr', keep_srrs)
    df['libsize'] = df[['libsize_R1', 'libsize_R2']].max(axis=1)
    return df.libsize.min(), df.libsize.median(), df.libsize.max()


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
    df = store['prealn/workflow/fastq'].copy()
    df.reset_index(inplace=True)
    df = remove_rows(df, 'srr', filter_srrs)
    df = keep_rows(df, 'srr', keep_srrs)
    df['len'] = df[['avgLen_R1', 'avgLen_R2']].max(axis=1)

    return df.loc[df['len'] >= cutoff, ['srx', 'srr']]


def readlen_cnts(store, filter_srrs=None, keep_srrs=None):
    """Output lists of sample that meet libsize criteria.

    Parameters
    ----------
    store : pd.HDFStore
        HDF5 data store.
    filter_srrs : list-like
        A list-like of SRRs to remove.

    Returns
    -------
    tuple of floats
       (minimum, median, mode, max) of libsize.

    """
    df = store['prealn/workflow/fastq'].copy()
    df.reset_index(inplace=True)
    df = remove_rows(df, 'srr', filter_srrs)
    df = keep_rows(df, 'srr', keep_srrs)
    df['len'] = df[['avgLen_R1', 'avgLen_R2']].max(axis=1)

    return df.len.min(), df.len.median(), df.len.mode().iloc[0], df.len.max()


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
    first = store['prealn/workflow/collectrnaseqmetrics/first'].copy()
    first.reset_index(inplace=True)
    second = store['prealn/workflow/collectrnaseqmetrics/second'].copy()
    second.reset_index(inplace=True)

    first = remove_rows(first, 'srr', filter_srrs)
    second = remove_rows(second, 'srr', filter_srrs)

    first = keep_rows(first, 'srr', keep_srrs)
    second = keep_rows(second, 'srr', keep_srrs)

    f = first.loc[first['PCT_CORRECT_STRAND_READS'] >= cutoff, ['srx', 'srr']]
    s = second.loc[second['PCT_CORRECT_STRAND_READS'] >= cutoff, ['srx', 'srr']]
    u = first[~(first.srr.isin(f.srr) | first.srr.isin(s.srr))]
    return f, s, u


def mappability(store, cutoff=.50, filter_srrs=None, keep_srrs=None):
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

    se = store['prealn/workflow/hisat2'][['num_reads', 'num_reads_unpaired', 'num_unaligned']].copy()
    se.dropna(inplace=True)
    se['prop_unaligned'] = se['num_unaligned'] / se['num_reads']

    pe = store['prealn/workflow/hisat2'][['num_reads', 'num_reads_paired', 'num_concordant_reads_unaligned']].copy()
    pe.dropna(inplace=True)
    pe['prop_unaligned'] = pe['num_concordant_reads_unaligned'] / pe['num_reads']

    df = pd.concat([se, pe])
    df.reset_index(inplace=True)

    df = remove_rows(df, 'srr', filter_srrs)
    df = keep_rows(df, 'srr', keep_srrs)

    return df.loc[df['prop_unaligned'] <= cutoff, ['srx', 'srr']]


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

    df = store['prealn/workflow/fastq_screen'].copy()
    df.reset_index(inplace=True)
    df = df[['srx', 'srr', 'reference', 'one_hit_one_library_percent']].set_index(['srx', 'srr', 'reference']).unstack()
    df.columns = df.columns.droplevel(0)
    df.reset_index(inplace=True)

    df = remove_rows(df, 'srr', filter_srrs)
    df = keep_rows(df, 'srr', keep_srrs)

    return df.loc[df['dm6'] >= cutoff, ['srx', 'srr']]

def calculate_multi_corr_method1(df, srx, method):
    """Calculate correlation to the median.

    Parameters
    ----------
    df : pd.DataFrame
        Data frame with genes as rows and samples as columns.
    srx : str
        SRA experiment accession (SRX).
    method : str
        Which method to use for calculating correlation. ['spearman' or 'pearson']

    Returns
    -------
    list of tuples
        Where the tuples are (srx, srr1:median, correlation).

    """
    corrs = []
    # Correlation function
    if method == 'spearman':
        corr_func = scipy.stats.spearmanr
    else:
        corr_func = scipy.stats.pearsonr

    # Calcualte median values
    med = df.median(axis=1)
    for col in df.columns:
        corrs.append([srx, ':'.join([col, 'median']), corr_func(df[col], med).correlation])

    return corrs


def calculate_multi_corr_method2(df, srx, method, cutoff=0.9):
    """Use pairwise correlations with cutoff.

    Uses pairwise correlations along with a cutoff to identify samples that
    have outliers.

    1) Calculates all pairwise correlations
    2) Removes pairwise correlations that are less than cutoff
    3) Re-calculates pairwise correlations among remaining samples.
    4) Returns list of tuples with (srx, srr1:srr2, correlation) for the
    remaining srrs.

    Note: A sample is kept if it has a correlation with any other sample above
    cutoff. It is possible that if there are two groups of samples that are not
    related, all of the samples would be kept.

    Parameters
    ----------
    df : pd.DataFrame
        Data frame with genes as rows and samples as columns.
    srx : str
        SRA experiment accession (SRX).
    method : str
        Which method to use for calculating correlation. ['spearman' or 'pearson']
    cutoff : float
        Cutoff to use to filter low correlated samples.

    Returns
    -------
    list of tuples
        Where the tuples are (srx, srr1:srr2, correlation).

    """
    # Figure out what samples are correlated above cutoff.
    corr = df.corr(method=method)
    np.fill_diagonal(corr.values, 0)
    stacked = corr.stack()
    filtered = stacked[stacked >= cutoff].index.get_level_values(1).unique()

    if filtered.shape[0] > 1:
        # Remove "bad" samples and recalculate correlation.
        corr = df[filtered].corr(method=method)

        # If we still have bad correlation then just call it bad.
        if corr.min().min() < cutoff - (cutoff * .1):
            return [(df.index.name, None, None)]

        # Remove duplicates, i.e. get ride of lower triangle of correlation matrix
        np.fill_diagonal(corr.values, 0)
        corr = pd.DataFrame(np.triu(corr.values), columns=corr.columns, index=corr.index)

        # Stack and format for output
        stacked = corr.stack().to_frame()
        stacked.columns = ['corr']
        stacked = stacked[stacked['corr'] > 0].copy()
        stacked['srrs'] = stacked.apply(lambda x: ':'.join(x.name), axis=1)
        stacked['srx'] = srx
        stacked.reset_index(drop=True, inplace=True)
        stacked.set_index('srx', inplace=True)
        return stacked[['srrs', 'corr']].to_records().tolist()
    else:
        return [(srx, None, None)]


def calculate_corr_among_srr(store, srx, method='spearman', show_warn=True, multi='m1', **kwargs):
    """Calculate correlation among SRR for an SRX.

    Parameters
    ----------
    store : pd.HDFStore
        HDFStore with feature_counts in it.
    srx : str
        SRA experiment accession (SRX).
    method : str
        Method to use for calculating correlation, either 'pearson' or 'spearman'.
    show_warn : bool
        If true then show warnings. If False then hide warnings.
    multi : str
        If 'm1' then use the correlation with the median. If 'm2' then use
        pairwise correlations.
    kwargs
        Additional kwargs which will be passed to calculate_multi_corr_method1
        and 2.

    Returns
    -------
    list of tuples
        Returns a list of tuples where each tuple contains (srx, srr1:srr2,
        correlation). Where srr1:srr2 is the comparison for which correlation
        was calculated. For SRX with >2 SRRs, srr2 is the median of all SRRs.

    """
    # Import data
    df = get_feature_counts(store, srx, show_warn=show_warn)

    # If no data then return None
    if df is None:
        return [(srx, None, None)]

    # Remove SRX from column index
    if isinstance(df.columns[0], tuple):
        df.columns = df.columns.droplevel(0)

    # Calculate correlation
    if df.shape[1] == 1:
        return [(srx, ':'.join([df.columns[0], 'None']), None)]

    elif df.shape[1] == 2:
        return [(srx, ':'.join(df.columns), df.corr(method=method).values[0, 1])]

    elif df.shape[1] > 2:
        corrs = []
        if multi == 'm2':
            cutoff = kwargs.get('cutoff', 0.9)
            corrs.extend(calculate_multi_corr_method2(df, srx, method, cutoff=cutoff))
        else:
            corrs.extend(calculate_multi_corr_method1(df, srx, method))
        return corrs

    # Something is wrong output warning and return None
    if show_warn:
        logger.warn('Something wrong with {}'.format(srx))

    return [(srx, None, None)]


def get_feature_counts(store, srx, show_warn=True):
    """Import feature counts tables for all SRR for a given SRX(s).

    Parameters
    ----------
    store : pd.HDFStore
        HDFStore with feature_counts in it.
    srx : str or list
        SRA experiment accession(s) (SRX).
    show_warn : bool
        If true then show warnings. If False then hide warnings.

    Returns
    -------
    pd.DataFrame | None
        If there are data then it outputs a pd.DataFrame, otherwise it returns None.

    """

    # Get data.
    df = store.select('prealn/workflow/feature_counts/counts', 'srx == srx').unstack().T

    if (df.shape == (0, 0)):
        if show_warn:
            logger.warn('Missing Feature Counts: {}'.format(srx))
        return None

    return df
