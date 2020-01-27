"""Distribution sample submssion by year."""
from typing import Union, Tuple

import pandas as pd
from pandas.plotting import register_matplotlib_converters
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.formula.api import ols

from ncbi_remap.plotting import NcbiPlotter, update_kwargs

register_matplotlib_converters()

PLOT_DEFAULTS = dict()
AXES_DEFAULTS = dict(ylabel="% Aligned", xlabel="Year")


def get_data(file_name: str) -> Tuple[pd.Series, pd.Series]:
    return (
        pd.read_table(file_name, usecols=["srx", "date_created", "pct_uniquely_aligned"])
        .assign(date_created=lambda x: pd.to_datetime(x.date_created))
        .assign(year=lambda x: x.date_created.dt.year)
        .set_index("srx")
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
        self.ax_joint = None
        self.ax_marg_x = None
        self.ax_marg_y = None

        self.plot_kwargs = update_kwargs(PLOT_DEFAULTS, plot_kwargs)
        self.axes_kwargs = update_kwargs(AXES_DEFAULTS, ax_kwargs)

        self.data = get_data(file_name)
        self.plot()
        self.tweak()

    def plot(self):
        sns.pointplot(
            x="year", y="pct_uniquely_aligned", data=self.data, ax=self.ax, **self.plot_kwargs
        )

    def tweak(self):
        self.ax.set(**self.axes_kwargs)
        plt.setp(self.ax.get_xticklabels(), rotation=45, ha="right")
        sns.despine(ax=self.ax, left=True)
