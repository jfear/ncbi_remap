"""Read Length boxplot by year.

Creates a boxplot of read length by year
"""
from typing import Union, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import NcbiPlotter, update_kwargs


PLOT_DEFAULTS = dict()
AXES_DEFAULTS = dict(ylabel="Read Length (bp)", xlabel="Year")


def get_data(file_name: str) -> Tuple[pd.Series, pd.Series]:
    df = (
        pd.read_table(file_name, usecols=["srx", "date_created", "avgLen_R1", "avgLen_R2"])
        .assign(date_created=lambda x: pd.to_datetime(x.date_created))
        .assign(year=lambda x: x.date_created.dt.year)
        .set_index(["srx", "date_created", "year"])
    )

    max_read_length = df.dropna().max(axis=1).rename("avg_read_length")

    return max_read_length.reset_index()


class Plot(NcbiPlotter):
    def __init__(
        self,
        file_name: str,
        plot_kwargs: Union[None, dict] = None,
        ax_kwargs: Union[None, dict] = None,
        ax: plt.Axes = None,
    ):
        """Plot Read Length.

        Example
        ----------
        >>> from ncbi_remap.plotting.read_length import Plot
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
        sns.boxplot(
            x="year", y="avg_read_length", data=self.data, ax=self.ax, **self.plot_kwargs
        )

    def tweak(self):
        self.ax.set(**self.axes_kwargs)
        sns.despine(ax=self.ax, left=True)
