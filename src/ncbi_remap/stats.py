import numpy as np
import pandas as pd

from scipy.stats import chi2_contingency, norm
from scipy.stats.contingency import margins
from statsmodels.stats.multitest import multipletests

def adjusted_residuals(observed, expected):
    n = observed.sum().sum()
    rsum, csum = margins(observed)
    v = csum * rsum * (n - rsum) * (n - csum) / n ** 3
    return (observed - expected) / np.sqrt(v)


def make_big_table(obs, expected, resid, adj_resid, cell_chisqs, cell_qvals, cell_flags):
    expected = pd.DataFrame(expected, index=obs.index, columns=obs.columns)
    cell_chisqs = pd.DataFrame(cell_chisqs, index=obs.index, columns=obs.columns)
    cell_qvals = pd.DataFrame(
        cell_qvals.reshape(resid.shape), index=adj_resid.index, columns=adj_resid.columns
    )
    cell_flags = pd.DataFrame(
        cell_flags.reshape(resid.shape), index=adj_resid.index, columns=adj_resid.columns
    )

    obs["type"] = "observed"
    obs = obs.set_index("type", append=True)

    expected["type"] = "expected"
    expected = expected.set_index("type", append=True)

    resid["type"] = "residual"
    resid = resid.set_index("type", append=True)

    adj_resid["type"] = "adj std residual"
    adj_resid = adj_resid.set_index("type", append=True)

    cell_chisqs["type"] = "X^2"
    cell_chisqs = cell_chisqs.set_index("type", append=True)

    cell_qvals["type"] = "fdr q-value"
    cell_qvals = cell_qvals.set_index("type", append=True)

    cell_flags["type"] = "flag_sig"
    cell_flags = cell_flags.set_index("type", append=True)

    _df = pd.concat(
        [obs, expected, resid, adj_resid, cell_chisqs, cell_qvals, cell_flags]
    ).reset_index(level="type")
    _df["type"] = pd.Categorical(
        _df["type"],
        ordered=True,
        categories=[
            "observed",
            "expected",
            "residual",
            "adj std residual",
            "X^2",
            "fdr q-value",
            "flag_sig",
        ],
    )
    return _df.set_index("type", append=True).sort_index().round(4)


def run_chisq(df, **kwargs):
    """ A helper function to run post hoc tests on chi^2.

    Parameters
    ----------
    df: a cross tabulation table
    kwargs: kwargs passed to statsmodels.stats.multitest.multipletests

    Returns
    -------
    str with chi-squre test results
    pandas.DataFrame with chi-square post-hoc tests.

    """
    obs = df.copy()
    stat, pval, degrees, expected = chi2_contingency(obs)
    chi2_res = f"chi^2: {stat:,.4f}, p-value: {pval:,.4f}, df: {degrees:,}"
    resid = obs - expected
    adj_resid = adjusted_residuals(df, expected)
    cell_chisqs = resid ** 2 / expected
    cell_flags, cell_qvals, _, _ = multipletests(
        norm.pdf(adj_resid).flatten(), method="fdr_bh", **kwargs
    )
    return chi2_res, make_big_table(obs, expected, resid, adj_resid, cell_chisqs, cell_qvals, cell_flags)
