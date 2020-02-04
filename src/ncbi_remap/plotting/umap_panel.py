from typing import Union, Tuple, Optional
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import NcbiPlotter, update_kwargs

PANEL_DEFAULTS = dict(col_wrap=4, palette=["lightgray", "C0", "C1"])
PLOT_DEFAULTS = dict(rasterized=True, s=8, linewidth=0.2)


class Plot(NcbiPlotter):
    def __init__(
        self,
        embeddings: pd.DataFrame,
        labels: pd.Series,
        outlier_flag: Optional[pd.Series] = None,
        panel_kwargs: Optional[dict] = None,
        plot_kwargs: Optional[dict] = None,
    ):
        """
        Example
        ----------
        >>> %load_ext autoreload
        >>> %autoreload 2

        >>> import pandas as pd
        >>> from ncbi_remap.plotting.umap_panel import Plot
        >>> embeddings = pd.read_parquet("../../../output/library_strategy-wf/umap_prealn_features_embeddings.tsv")
        >>> labels = pd.read_parquet("../../../output/library_strategy-wf/prealn_outliers.parquet")
        >>> plot = Plot(
        ...     embeddings=embeddings,
        ...     labels=labels.library_strategy,
        ...     outlier_flag=labels.library_strategy_flag_outlier
        ... )

        """

        self.update_figsize()

        self.embeddings = embeddings
        self.embedding_columns = embeddings.columns
        self.labels = labels
        self.label_source = labels.name
        self.outlier_flag = outlier_flag
        self.panel_kwargs = update_kwargs(PANEL_DEFAULTS, panel_kwargs)
        self.plot_kwargs = update_kwargs(PLOT_DEFAULTS, plot_kwargs)

        self.data = None
        self.panel = None

        self.wrangle()
        self.plot()
        self.tweak()

    def wrangle(self):
        self.data = (
            pd.get_dummies(self.labels)
            .reset_index()
            .melt(id_vars="srx", var_name=self.label_source, value_name="dummy")
            .set_index("srx")
        ).join(self.embeddings, how="inner")

        if self.outlier_flag is None:
            return

        # Make outliers 2
        outlier_srxs = self.outlier_flag[self.outlier_flag].index
        mask = (self.data.index.isin(outlier_srxs)) & (self.data.dummy == 1)
        self.data.loc[mask, "dummy"] = 2

    def plot(self):
        x, y = self.embedding_columns
        self.panel = sns.FacetGrid(
            data=self.data, col=self.label_source, hue="dummy", **self.panel_kwargs
        )
        self.panel.map(sns.scatterplot, x, y, **self.plot_kwargs)

    def tweak(self):
        self.panel.set_titles("{col_name}")
