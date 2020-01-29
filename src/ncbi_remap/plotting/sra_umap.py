from typing import Union, Tuple, List

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import umap

from ncbi_remap.plotting import NcbiPlotter, update_kwargs

PLOT_DEFAULTS = dict(palette="spectral")
AXES_DEFAULTS = dict(ylabel="UMAP 2", xlabel="UMAP 1")


def get_features(file_name: str) -> Tuple[pd.Series, pd.Series]:
    return pd.read_parquet(file_name).fillna(0)


def get_labels(file_name: str, srxs: list) -> List[str]:
    return pd.read_parquet(file_name).library_strategy.reindex(srxs)


class Plot(NcbiPlotter):
    def __init__(
        self,
        features_file: str,
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
        ... features_file="../../../output/library_strategy-wf/prealn_feature_set.parquet",
        ... labels_file="../../../output/library_strategy-wf/sra_strategy_selection.parquet")
        """
        self.update_figsize()
        self.ax = ax or self.get_ax()

        self.plot_kwargs = update_kwargs(PLOT_DEFAULTS, plot_kwargs)
        self.axes_kwargs = update_kwargs(AXES_DEFAULTS, ax_kwargs)

        self.features = get_features(features_file)
        self.labels = get_labels(labels_file, self.features.index)
        self.run_umap()

        # self.plot()
        # self.tweak()

    def run_umap(self):
        self.mapper = umap.UMAP().fit(self.features)

    def plot(self):
        plt.scatter(x=self.x, y=self.y, c=self.labels, ax=self.ax, **self.plot_kwargs)

    def tweak(self):
        self.ax.set(**self.axes_kwargs)
        sns.despine(ax=self.ax, left=True)
