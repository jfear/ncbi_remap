"""Distribution sample submssion by year."""
from typing import Union, Tuple

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.formula.api import ols

from ncbi_remap.plotting import NcbiPlotter, update_kwargs


PLOT_DEFAULTS = dict()
AXES_DEFAULTS = dict(ylabel="% Genes On (> 5 reads)", xlabel="Samples (SRX)")


def get_data(file_name: str) -> pd.Series:
    return (
        (pd.read_table(file_name, index_col=0) > 5)
        .mean()
        .mul(100)
        .sort_values(ascending=True)
        .rename("prop_genes_on")
        .to_frame()
        .assign(Samples=lambda x: range(x.shape[0]))
    )


class Plot(NcbiPlotter):
    def __init__(
        self,
        file_name: str,
        plot_kwargs: Union[None, dict] = None,
        ax_kwargs: Union[None, dict] = None,
        ax: plt.Axes = None,
    ):
        """
        Example
        ----------
        >>> from ncbi_remap.plotting.genes_on import Plot
        >>> plot = Plot("../../../output/paper-wf/srx_prealn.tsv")
        """
        self.update_figsize()
        self.ax = ax or self.get_ax()

        self.plot_kwargs = update_kwargs(PLOT_DEFAULTS, plot_kwargs)
        self.axes_kwargs = update_kwargs(AXES_DEFAULTS, ax_kwargs)

        self.data = get_data(file_name)
        self.plot()
        self.tweak()

    def plot(self):
        plt.plot(self.data.Samples, self.data.prop_genes_on)

    def tweak(self):
        self.ax.set(**self.axes_kwargs)
        sns.despine(ax=self.ax, left=True)
