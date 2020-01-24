"""Distribution sample submssion by year."""
from typing import Union, Tuple

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.formula.api import ols

from ncbi_remap.plotting import NcbiPlotter, update_kwargs


PLOT_DEFAULTS = dict()
AXES_DEFAULTS = dict(ylabel="Samples (Thousands)", xlabel="Year")


def get_data(file_name: str) -> Tuple[pd.Series, pd.Series]:
    return (
        pd.read_table(file_name, usecols=["srx", "date_created"])
        .assign(date_created=lambda x: pd.to_datetime(x.date_created))
        .assign(year=lambda x: x.date_created.dt.year)
        .groupby("year")
        .size()
        .rename("num_samples")
        .reset_index()
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
        >>> from ncbi_remap.plotting.sample_submission import Plot
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
        sns.regplot(x="year", y="num_samples", data=self.data, ax=self.ax, **self.plot_kwargs)

    def tweak(self):
        self.ax.set(**self.axes_kwargs)
        self.ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, pos: f"{y/1000:.0f}"))
        sns.despine(ax=self.ax, left=True)

    def stats(self):
        model = ols("num_samples ~ year", data=df).fit()
        model.summary()
