""""""
from typing import Union, Tuple, List
from functools import partial

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import NcbiPlotter, update_kwargs


PLOT_DEFAULTS = dict(hist=False)
AXES_DEFAULTS = dict(xlabel="% Same Strand", ylabel="")


def get_data(file_name: str, threshold) -> Tuple[pd.Series, pd.Series]:
    df = (
        pd.read_table(file_name, usecols=["srx", "pct_first_strand_reads"], index_col="srx")
        .squeeze()
        .dropna()
    )

    counts = pd.Series([
        (threshold <= df).sum(),
        (((100 - threshold) < df) & (df < threshold)).sum(),
        (df <= (100 - threshold)).sum(),
    ], index=["First", "Unstranded", "Second"])

    return df, counts


class Plot(NcbiPlotter):
    def __init__(
        self,
        file_name: str,
        plot_kwargs: Union[None, dict] = None,
        ax_kwargs: Union[None, dict] = None,
        threshold: Union[int, float] = 75,
    ):
        """
        Example
        ----------
        >>> from ncbi_remap.plotting.strandedness import Plot
        >>> plot = Plot("../../../output/paper-wf/srx_prealn.tsv")
        """
        self.update_figsize()
        self.ax = self.get_ax()

        self.plot_kwargs = update_kwargs(PLOT_DEFAULTS, plot_kwargs)
        self.axes_kwargs = update_kwargs(AXES_DEFAULTS, ax_kwargs)

        self.threshold = threshold

        self.data, self.counts = get_data(file_name, self.threshold)
        self.plot()
        self.tweak()

    def plot(self):
        sns.distplot(self.data, ax=self.ax, **self.plot_kwargs)

    def tweak(self):
        self.ax.set(**self.axes_kwargs)

        _, ylim = self.ax.get_ylim()
        ylim = ylim + ylim * .1
        self.ax.set_ylim(0, ylim)

        yloc = .96
        
        # Upper
        xloc = .8
        self.ax.text(xloc, yloc, f"{self.counts['First']:,}", ha="center", va="top", transform=self.ax.transAxes)
        self.ax.axvline(self.threshold, ls='--')

        # Middle
        xloc = .5
        self.ax.text(xloc, yloc, f"{self.counts['Unstranded']:,}", ha="center", va="top", transform=self.ax.transAxes)

        # Lower
        xloc = .2
        self.ax.text(xloc, yloc, f"{self.counts['Second']:,}", ha="center", va="top", transform=self.ax.transAxes)
        self.ax.axvline(100 - self.threshold, ls='--')
