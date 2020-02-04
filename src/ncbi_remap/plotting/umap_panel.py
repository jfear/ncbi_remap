from typing import Union, Tuple
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import NcbiPlotter, update_kwargs

PANEL_DEFAULTS = dict(col_wrap=4, palette=["lightgray", "C0"])
PLOT_DEFAULTS = dict(rasterized=True, s=8, linewidth=.2)


class Plot(NcbiPlotter):
    def __init__(
        self,
        embeddings_file: str,
        labels_file: str,
        column: str = "library_strategy",
        panel_kwargs: Union[None, dict] = None,
        plot_kwargs: Union[None, dict] = None,
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

        self.embeddings_file = embeddings_file
        self.labels_file = labels_file
        self.column = column
        self.panel_kwargs = update_kwargs(PANEL_DEFAULTS, panel_kwargs)
        self.plot_kwargs = update_kwargs(PLOT_DEFAULTS, plot_kwargs)

        self.embedding_columns = None
        self.data = None
        self.panel = None

        self.get_data()

        self.plot()
        self.tweak()

    def get_data(self):
        labels = pd.read_parquet(self.labels_file)
        embeddings = pd.read_table(self.embeddings_file, index_col=0)

        self.embedding_columns = embeddings.columns
        self.data = (
            pd.get_dummies(labels[self.column])
            .reset_index()
            .melt(id_vars="srx", var_name=self.column, value_name="dummy")
            .set_index("srx")
        ).join(embeddings, how="inner")

    def plot(self):
        x, y = self.embedding_columns
        self.panel = sns.FacetGrid(
            data=self.data, col=self.column, hue="dummy", **self.panel_kwargs
        )
        self.panel.map(sns.scatterplot, x, y, **self.plot_kwargs)

    def tweak(self):
        self.panel.set_titles("{col_name}")
