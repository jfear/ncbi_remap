"""Library layout bar graph.

Creates a bargraph of library layout counts and adds simple annotation
"""
from typing import Union, Tuple

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import NcbiPlotter, update_kwargs

ORDER = ["SE", "PE", "PE->SE"]
OFFSET = 300

PLOT_DEFAULTS = dict(order=ORDER)
AXES_DEFAULTS = dict(yticks=[], ylabel="# SRX", xlabel="", title="Library Layout")
ANNOT_DEFAULTS = dict(ha="center")


def get_data(file_name: str) -> Tuple[pd.Series, pd.Series]:
    """Import library layout data."""
    data = (
        pd.read_parquet(file_name, columns=["library_layout"])
        .squeeze()
        .replace({"keep_R1": "PE->SE", "keep_R2": "PE->SE"})
    )

    counts = data.value_counts()

    return data, counts

class LayoutPlot(NcbiPlotter):
    def __init__(
        self,
        file_name: str,
        plot_kwargs: Union[None, dict] = None,
        ax_kwargs: Union[None, dict] = None,
        annot_kwargs: Union[None, dict] = None,
        ax: plt.Axes = None,
    ):
        """Plot library layout.

        Example
        ----------
        >>> from ncbi_remap.plotting import LayoutPlot
        >>> layout_plot = LayoutPlot("../../../output/paper-wf/srx_data.parquet")
        """
        self.ax = ax or self.get_ax()
        self.plot_kwargs = update_kwargs(PLOT_DEFAULTS, plot_kwargs)
        self.axes_kwargs = update_kwargs(AXES_DEFAULTS, ax_kwargs)
        self.annot_kwargs = update_kwargs(ANNOT_DEFAULTS, annot_kwargs)

        self.data, self.counts = get_data(file_name)
        self.plot()
        self.add_annotation()
        self.tweak()


    def plot(self):
        sns.countplot(self.data, ax=self.ax, **self.plot_kwargs)

    def add_annotation(self):
        counts = self.counts
        kwargs = self.annot_kwargs
        self.ax.text(0, counts["SE"] + OFFSET, f"{counts['SE']:,}", **kwargs)
        self.ax.text(1, counts["PE"] + OFFSET, f"{counts['PE']:,}", **kwargs)
        self.ax.text(2, counts["PE->SE"] + OFFSET, f"{counts['PE->SE']:,}", **kwargs)

    def tweak(self):
        self.ax.set(ylim=(0, self.counts.max() + OFFSET * 10), **self.axes_kwargs)
        self.ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, pos: f"{x:,d}"))
        sns.despine(ax=self.ax, left=True, bottom=True)
