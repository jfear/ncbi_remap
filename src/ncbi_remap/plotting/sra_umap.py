from typing import Union, Optional

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import NcbiPlotter, update_kwargs

PLOT_DEFAULTS = dict(alpha=0.5, rasterized=True, palette=["C0", "C1", "lightgray", "r"], s=8)
AXES_DEFAULTS = dict(ylabel="UMAP 2", xlabel="UMAP 1")
CATEGORIES = ["RNA-Seq", "ChIP-Seq", "Other", "Outlier"]


def update_labels(x):
    if x in CATEGORIES:
        return x
    else:
        return "Other"


class Plot(NcbiPlotter):
    def __init__(
        self,
        embeddings: pd.DataFrame,
        labels: pd.Series,
        outlier_flag: Optional[pd.Series],
        plot_kwargs: Union[None, dict] = None,
        ax_kwargs: Union[None, dict] = None,
        ax: plt.Axes = None,
    ):
        """
        Example
        ----------
        >>> %load_ext autoreload
        >>> %autoreload 2

        >>> from ncbi_remap.plotting.sra_umap import Plot
        >>> plot = Plot(
        ...     embeddings_file="../../../output/library_strategy-wf/umap_prealn_features_embeddings.tsv",
        ...     labels_file="../../../output/library_strategy-wf/sra_strategy_selection.parquet",
        ... )

        """

        self.update_figsize()
        self.ax = ax or self.get_ax()

        self.embeddings = embeddings
        self.labels = labels
        self.outlier_flag = outlier_flag
        self.plot_kwargs = update_kwargs(PLOT_DEFAULTS, plot_kwargs)
        self.axes_kwargs = update_kwargs(AXES_DEFAULTS, ax_kwargs)

        self.data = None

        self.wrangle()
        self.plot()
        self.tweak()

    def wrangle(self):
        labels = self.labels.map(update_labels)
        self.data = self.embeddings.join(labels, how="inner")

        if self.outlier_flag is None:
            return

        outlier_srxs = self.outlier_flag[self.outlier_flag].index
        self.data.loc[self.data.index.isin(outlier_srxs), self.labels.name] = "Outlier"

    def plot(self):
        """Plot scatter plot of RNA-Seq, ChIP-Seq, Other.

        I plot each one separately, so I can set the z-order
        """
        self.zorder = {k: v for k, v in zip(CATEGORIES, [10, 5, 0, 11])}
        self.colors = {k: v for k, v in zip(CATEGORIES, self.plot_kwargs.pop("palette"))}
        for category in CATEGORIES:
            self.plot_scatter(category)
        self.ax.legend(loc="upper center", bbox_to_anchor=[1, 1])

    def plot_scatter(self, category):
        dd = self.data.query(f"{self.labels.name} == '{category}'")
        sns.scatterplot(
            dd.UMAP1,
            dd.UMAP2,
            label=category,
            zorder=self.zorder[category],
            color=self.colors[category],
            ax=self.ax,
            **self.plot_kwargs,
        )

    def tweak(self):
        self.ax.set(**self.axes_kwargs)
        sns.despine(ax=self.ax, left=True, bottom=True)
