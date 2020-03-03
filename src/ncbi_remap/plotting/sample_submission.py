"""Distribution sample submssion by year."""
from typing import Union, Tuple

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.formula.api import ols

from ncbi_remap.plotting import NcbiPlotter, update_kwargs


BAR_PLOT_DEFAULTS = dict()
REG_PLOT_DEFAULTS = dict()
BAR_AXES_DEFAULTS = dict(ylabel="Cumulative Samples\n(Thousands)", xlabel="Year")
REG_AXES_DEFAULTS = dict(ylabel="Samples (Thousands)", xlabel="Year")


def get_data(file_name: str) -> Tuple[pd.Series, pd.Series]:
    return (
        pd.read_table(file_name, usecols=["srx", "date_created"])
        .assign(date_created=lambda x: pd.to_datetime(x.date_created))
        .assign(year=lambda x: x.date_created.dt.year)
        .groupby("year")
        .size()
        .rename("num_samples")
        .reset_index()
        .assign(cum_sum=lambda x: x.num_samples.cumsum())
    )


class Plot(NcbiPlotter):
    def __init__(
        self,
        file_name: str,
        bar_plot_kwargs: Union[None, dict] = None,
        reg_plot_kwargs: Union[None, dict] = None,
        bar_ax_kwargs: Union[None, dict] = None,
        reg_ax_kwargs: Union[None, dict] = None,
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
        self.ax2 = self.ax.twinx()
        self.bar_plot_kwargs = update_kwargs(BAR_PLOT_DEFAULTS, bar_plot_kwargs)
        self.reg_plot_kwargs = update_kwargs(REG_PLOT_DEFAULTS, reg_plot_kwargs)
        self.bar_axes_kwargs = update_kwargs(BAR_AXES_DEFAULTS, bar_ax_kwargs)
        self.reg_axes_kwargs = update_kwargs(REG_AXES_DEFAULTS, reg_ax_kwargs)

        self.data = get_data(file_name)
        self.plot()
        self.add_fill()
        self.tweak()

    def plot(self):
        sns.lineplot(x="year", y="cum_sum", data=self.data, ax=self.ax, **self.bar_plot_kwargs)
        sns.regplot(x="year", y="num_samples", data=self.data, ax=self.ax2, **self.reg_plot_kwargs)

    def add_fill(self):
        x = self.data.year
        y = self.data.cum_sum
        self.ax.fill_between(x, y, color=self.bar_plot_kwargs.get("color", "C0"))

    def tweak(self):
        self.ax.set(**self.bar_axes_kwargs)
        self.ax2.set(**self.reg_axes_kwargs)
        self.ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, pos: f"{y/1000:.0f}"))
        self.ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, pos: f"{y/1000:.0f}"))
        sns.despine(ax=self.ax, left=True)
        sns.despine(ax=self.ax2, left=True)

        self.color_axis(self.ax, self.bar_plot_kwargs["color"])
        self.color_axis(self.ax2, self.reg_plot_kwargs["color"])

    def stats(self):
        model = ols("num_samples ~ year", data=df).fit()
        model.summary()

    @staticmethod
    def color_axis(ax, color):
        ax.tick_params(axis="y", color=color, labelcolor=color)
        ax.yaxis.label.set_color(color)
