""""""
from typing import Union, Tuple, List
from functools import partial

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import NcbiPlotter, update_kwargs


STRATEGY_PLOT_DEFAULTS = dict(color="lightgray", edgecolor="k")
SELECTION_PLOT_DEFAULTS = dict(color="lightgray", edgecolor="k")
STRATEGY_AXES_DEFAULTS = dict(xlabel="Importance", ylabel="Feature", title="Strategy")
SELECTION_AXES_DEFAULTS = dict(xlabel="Importance", ylabel="", title="Selection")


def get_data(file_name: str) -> Tuple[pd.Series, pd.Series]:
    return (
        pd.read_table(
            file_name,
            header=None,
            names=[
                "Feature",
                "Strategy",
                "Selection",
            ]
        )
        .assign(Feature=lambda x: x.Feature.str.lower())
        .sort_values("Strategy", ascending=False)
        .head(20)
    )


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
        self.update_figsize()
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
        sns.barplot("Strategy", "Feature", data=self.data, ax=self.ax1, **self.strategy_plot_kwargs)
        sns.barplot("Selection", "Feature", data=self.data, ax=self.ax2, **self.selection_plot_kwargs)

    def tweak(self):
        self.ax1.set(**self.strategy_axes_kwargs)
        self.ax2.set(**self.selection_axes_kwargs)

