"""Output lists of samples meeting different criteria."""
from itertools import combinations, product
import numpy as np
import pandas as pd
import scipy
import scipy.stats
from .logging import logger
from .normalization import cpm

# Constants
LIBSIZE_CUTOFF = 1e5      # 100,000 reads
READLEN_CUTOFF = 30       # 30 bp
STRAND_CUTOFF1 = 0.75     # 75% stranded reads
STRAND_CUTOFF2 = 0.95     # 95% stranded reads
UNALIGN_CUTOFF = 0.50     # 50% Reads Unaligned
CONTAMINATION_CUTOFF = 40 # 50% reads mapping to DM6
SPEARMAN_CUTOFF = .95  # Minimum spearman r


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
    u = first.loc[~(first.srr.isin(f.srr) | first.srr.isin(s.srr)), ['srx', 'srr']]
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


def srx_reproducibility_score(store, srx, method='spearman', multi='pairwise', TH=1, show_warn=True, **kwargs):
    """Calculate reproducibility among SRR for an SRX.

    This can either be the correlation or the SERE.

    Parameters
    ----------
    store : pd.HDFStore
        HDFStore with feature_counts in it.
    srx : str
        SRA experiment accession (SRX).
    method : str
        Method to use for calculating reproducibility. Option: 'pearson',
        'spearman', 'sere'.
    multi : str
        If 'median' then use the correlation with the median. Else do pairwise
        correlations.
    TH : int
        Minimum number of reads for a gene to be considered expressed.
    show_warn : bool
        If true then show warnings. If False then hide warnings.
    kwargs
        Additional kwargs which will be passed to either sere_iter or corr_iter
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

    # If only one column then return None
    if df.shape[1] == 1:
        return [(srx, ':'.join([df.columns[0], 'None']), None)]

    # more than two columns then output score
    if method == 'sere':
        return [(srx, pair, score) for pair, score in sere_iter(df, TH=TH, **kwargs)]

    if method == 'spearman' or method == 'pearson':
        return [(srx, pair, score) for pair, score in corr_iter(df, TH=TH, method=method, multi=multi, **kwargs)]

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


def corr_iter(df, method='spearman', multi='median', TH=1, normalize=False, normalize_kw=None, **kwargs):
    """Generator function to provide list of correlations.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with counts.
    method : str
        Method to use for calculating reproducibility. Option: 'pearson',
        'spearman'.
    multi : str
        If 'median' then use the correlation with the median. Else do pairwise
        correlations.
    TH : int
        Minimum number of reads for a gene to be considered expressed.
    normalize : None or function
        Normalization function to use.
    normalize_kw : None or dict
        kwargs to pass to the normalization function.

    Yields
    ------
    tuples
        Yields a tuples where each tuple contains (srr1:srr2, correlation).
        Where srr1:srr2 is the comparison for which correlation was calculated.

    """
    # Normalize
    if normalize_kw is None:
        normalize_kw = {}

    if normalize:
        _df = normalize(df, **normalize_kw)
    else:
        _df = df.copy()

    columns = _df.columns

    # Filter out low expressing genes
    mask = _df.sum(axis = 1) > TH
    _df = _df[mask]

    # Either compare to median or pairwise.
    if multi == 'median':
        _df['med'] = _df.median(axis=1)
        combos = product(columns, ['med'])
    else:
        combos = combinations(columns, 2)

    # Which correlation to calculate
    if method =='spearman':
        func = scipy.stats.spearmanr
    elif method == 'pearson':
        func = scipy.stats.pearsonr

    for x1, x2 in combos:
        score = func(_df[x1], _df[x2])[0]
        yield '{}:{}'.format(x1, x2), score


def sere_score(df, TH=1):
    """Calculate SERE (Single-paramter quality control and sample comparison for RNA-Seq).

    A method to compare RNA-seq samples and estimates reproducibility.
    Replicates should have a SERE value of near 1. This code is adapted from an
    R version [here](https://github.com/eigenv/SERE/blob/master/sere.R). Which is based on the paper:

    Schulze, Stefan K., Rahul Kanwar, Meike Gölzenleuchter, Terry M. Therneau,
    and Andreas S. Beutler. 2012. “SERE: Single-Parameter Quality Control and
    Sample Comparison for RNA-Seq.” BMC Genomics 13 (October):524.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with raw RNA-seq counts.
    TH : int
        Minimum number of reads for a gene to be considered expressed.

    Returns
    -------
    float
        SERE score.
    """
    num_samples = df.shape[1]

    total_counts = df.sum().sum()
    sample_counts = df.sum()
    row_counts = df.sum(axis=1)

    # Filter out non-expressed genes
    mask = row_counts > TH
    row_counts = row_counts[mask]
    df = df[mask]

    num_genes = df.shape[0]

    # Make expected count table
    def _estimate(idx):
        return row_counts * sample_counts.iloc[idx] / total_counts

    expected_counts = pd.concat([_estimate(x) for x in range(num_samples)], axis=1)
    expected_counts.columns = df.columns

    # sere score
    dispersion = ((df - expected_counts) ** 2 / expected_counts).sum().sum()
    sere = np.sqrt(dispersion / (num_genes * (num_samples - 1)))

    return sere


def sere_matrix(df, TH=1):
    """Calculate pairwise SERE.

    Uses SERE to estimate reproducibility for all pairwise combinations.
    This code is adapted from an R version
    [here](https://github.com/eigenv/SERE/blob/master/sere.R).

    Parameters
    -----------
    df : pd.DataFrame
        DataFrame with raw RNA-seq counts.
    TH : int
        Minimum number of reads for a gene to be considered expressed.

    Returns
    -------
    pd.DataFrame
        Matrix of SERE scores for all pairwise comparisons.

    """
    number_samples = df.shape[1]
    columns = df.columns

    # Distance matrix
    distance = np.full((num_samples, num_samples), np.nan)
    for i in range(num_samples):
        for j in range(i, num_samples):
            distance.iloc[i, j] = sere_score(df.loc[:, [i, j]], TH)
            distance.iloc[j, i] = distance.iloc[i, j]

    return pd.DataFrame(distance, index=columns, columns=columns)


def sere_iter(df, TH=1, **kwargs):
    """Generator function to provide list of SERE scores.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with counts.
    TH : int
        Minimum number of reads for a gene to be considered expressed.

    Yields
    ------
    tuples
        Yields a tuples where each tuple contains (srr1:srr2, SERE).
        Where srr1:srr2 is the comparison for which SERE was calculated.

    """
    columns = df.columns
    combos = combinations(columns, 2)
    for x1, x2 in combos:
        score = sere_score(df[[x1, x2]], TH=TH)
        yield '{}:{}'.format(x1, x2), score

