""""""
from typing import Union, Tuple, List
from functools import partial

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import NcbiPlotter, update_kwargs


STRATEGY_PLOT_DEFAULTS = dict(palette=["black", "red"], s=12, linewidth=.1, rasterized=True, legend=None)
SELECTION_PLOT_DEFAULTS = dict(palette=["black", "red"], s=12, linewidth=.1, rasterized=True, legend=None)
STRATEGY_AXES_DEFAULTS = dict(xlabel="Samples", xticklabels=[], ylabel="-(Isolation Score)", title="Strategy")
SELECTION_AXES_DEFAULTS = dict(xlabel="Samples", xticklabels=[], ylabel="", title="Selection")


def get_data(file_name: str) -> pd.DataFrame:
    return pd.read_parquet(file_name)


class Plot(NcbiPlotter):
    def __init__(
        self,
        file_name: str,
        strategy_plot_kwargs: Union[None, dict] = None,
        selection_plot_kwargs: Union[None, dict] = None,
        strategy_axes_kwargs: Union[None, dict] = None,
        selection_axes_kwargs: Union[None, dict] = None,
    ):
        """
        Example
        ----------
        >>> from ncbi_remap.plotting.library_size import Plot
        >>> plot = Plot("../../../output/paper-wf/srx_prealn.tsv")
        """
        self.update_figsize(wmul=2)
        self.fig_ncols = 2
        self.ax1, self.ax2 = self.get_ax(sharey=True)

        self.strategy_plot_kwargs = update_kwargs(STRATEGY_PLOT_DEFAULTS, strategy_plot_kwargs)
        self.selection_plot_kwargs = update_kwargs(SELECTION_PLOT_DEFAULTS, selection_plot_kwargs)

        self.strategy_axes_kwargs = update_kwargs(STRATEGY_AXES_DEFAULTS, strategy_axes_kwargs)
        self.selection_axes_kwargs = update_kwargs(SELECTION_AXES_DEFAULTS, selection_axes_kwargs)

        self.data = get_data(file_name)
        self.plot()
        self.tweak()

    def plot(self):
        x = range(self.data.shape[0])

        sns.scatterplot(
            x,
            self.data.library_strategy_outlier_score * -1,
            hue=self.data.library_strategy_flag_outlier,
            ax=self.ax1,
            **self.strategy_plot_kwargs
        )

        sns.scatterplot(
            x,
            self.data.library_selection_outlier_score * -1,
            hue=self.data.library_selection_flag_outlier,
            ax=self.ax2,
            **self.selection_plot_kwargs
        )

    def tweak(self):
        self.ax1.set(**self.strategy_axes_kwargs)
        sns.despine(ax=self.ax1, left=True)

        self.ax2.set(**self.selection_axes_kwargs)
        sns.despine(ax=self.ax2, left=True)

