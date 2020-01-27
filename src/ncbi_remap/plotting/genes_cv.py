"""Distribution sample submssion by year."""
from typing import Union, Tuple

import pandas as pd
from scipy.stats import variation
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.formula.api import ols

from ncbi_remap.plotting import NcbiPlotter, update_kwargs


PLOT_DEFAULTS = dict()
AXES_DEFAULTS = dict(ylabel="Coefficient of Variation (Genes)", xlabel="FBgn")


def get_data(file_name: str) -> pd.Series:
    return (
        pd.read_table(file_name, index_col=0)
        .T.assign(cv=lambda x: x.apply(variation, axis=1))
        .cv.sort_values(ascending=True)
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
        >>> from ncbi_remap.plotting.genes_cv import Plot
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
        plt.plot(self.data, **self.plot_kwargs)

    def tweak(self):
        self.ax.set(**self.axes_kwargs)
        sns.despine(ax=self.ax, left=True)
