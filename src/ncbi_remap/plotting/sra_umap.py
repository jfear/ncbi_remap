from typing import Union
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import NcbiPlotter, update_kwargs

PLOT_DEFAULTS = dict(
    alpha=0.5,
    rasterized=True,
    palette=["C0", "C1", "lightgray"],
    s=8,
)
AXES_DEFAULTS = dict(ylabel="UMAP 2", xlabel="UMAP 1")
CATEGORIES = ["RNA-Seq", "ChIP-Seq", "Other"]


def get_data(embeddings_file, labels_file) -> pd.DataFrame:
    embeddings = pd.read_table(embeddings_file, index_col=0)
    labels = pd.read_parquet(labels_file)
    return embeddings.join(labels, how="inner")


class Plot(NcbiPlotter):
    def __init__(
        self,
        embeddings_file: str,
        labels_file: str,
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

        self.plot_kwargs = update_kwargs(PLOT_DEFAULTS, plot_kwargs)
        self.axes_kwargs = update_kwargs(AXES_DEFAULTS, ax_kwargs)

        self.data = get_data(embeddings_file, labels_file)

        self.add_color_categories()
        self.plot()
        self.tweak()

    def add_color_categories(self):
        self.data["color"] = "Other"
        for category in CATEGORIES[:-1]:
            self.data.loc[self.data.library_strategy == category, "color"] = category

    def plot(self):
        """Plot scatter plot of RNA-Seq, ChIP-Seq, Other.

        I plot each one separately, so I can set the z-order
        """
        self.zorder = {k: v for k, v in zip(CATEGORIES, [10, 5, 0])}
        self.colors = {k: v for k, v in zip(CATEGORIES, self.plot_kwargs.pop("palette"))}
        for category in CATEGORIES:
            self.plot_scatter(category)
        self.ax.legend()

    def plot_scatter(self, category):
        dd = self.data.query(f"color == '{category}'")
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
