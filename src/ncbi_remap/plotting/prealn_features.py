""""""
from typing import Union, Tuple, List
from functools import partial

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import NcbiPlotter, update_kwargs


GRID_DEFAULTS = dict(vars=["% Coding (Log)", "% UTR (Log)", "% Intronic (Log)", "% Intergenic (Log)"])
DIAG_PLOT_DEFAULTS = dict()
OFFDIAG_PLOT_DEFAULTS = dict(shade=True, shade_lowest=False)
DIAG_AXES_DEFAULTS = dict()
OFFDIAG_AXES_DEFAULTS = dict()


def get_data(file_name: str) -> Tuple[pd.Series, pd.Series]:
    return (
        pd.read_table(
            file_name,
            usecols=[
                "srx",
                "pct_coding_bases",
                "pct_utr_bases",
                "pct_intronic_bases",
                "pct_intergenic_bases",
            ],
            index_col="srx"
        )
        .rename(columns=lambda x: rename_col(x))
        .fillna(0)
        .apply(lambda x: np.log10(x+1), axis=1)
    )

def rename_col(x: str):
    _, name, _ = x.split("_")
    if name == "utr":
        name = "UTR"
    else:
        name = name.title()
    return f"% {name} (Log)"

class Plot(NcbiPlotter):
    def __init__(
        self,
        file_name: str,
        grid_kwargs: Union[None, dict] = None,
        diag_plot_kwargs: Union[None, dict] = None,
        offdiag_plot_kwargs: Union[None, dict] = None,
    ):
        """
        Example
        ----------
        >>> from ncbi_remap.plotting.library_size import Plot
        >>> plot = Plot("../../../output/paper-wf/srx_prealn.tsv")
        """
        self.update_figsize()
        self.graph = None
        self.axes = None

        self.grid_kwargs = update_kwargs(GRID_DEFAULTS, grid_kwargs)
        self.diag_plot_kwargs = update_kwargs(DIAG_PLOT_DEFAULTS, diag_plot_kwargs)
        self.offdiag_plot_kwargs = update_kwargs(OFFDIAG_PLOT_DEFAULTS, offdiag_plot_kwargs)

        self.data = get_data(file_name)
        self.plot()

    def plot(self):
        g = sns.PairGrid(self.data, **self.grid_kwargs) # type: sns.PairGrid
        g.map_offdiag(sns.kdeplot, **self.offdiag_plot_kwargs)
        g.map_diag(plt.hist, **self.diag_plot_kwargs)
        self.graph = g
        self.axes = self.graph.axes # type: List[plt.Axes]
