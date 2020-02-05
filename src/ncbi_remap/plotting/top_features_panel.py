from typing import Union, Tuple, Optional

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import NcbiPlotter, update_kwargs

STRATEGY_ORDER = ["Other", "RNA-Seq", "ChIP-Seq"]
PANEL_DEFAULTS = dict(
    sharex=False,
    sharey=False,
    col_wrap=3,
    hue_order=STRATEGY_ORDER,
    palette=["lightgray", "C0", "C1"],
)
PLOT_DEFAULTS = dict()


def simplify_strategy(x):
    if (x == "RNA-Seq") | (x == "ChIP-Seq"):
        return x
    return "Other"


class Plot(NcbiPlotter):
    def __init__(
        self,
        features_file: str,
        importance_file: str,
        labels_file: str,
        panel_kwargs: Optional[dict] = None,
        plot_kwargs: Optional[dict] = None,
    ):
        """
        Example
        ----------
        >>> %load_ext autoreload
        >>> %autoreload 2

        >>> import pandas as pd
        >>> from ncbi_remap.plotting.top_features_panel import Plot
        >>> plot = Plot(
                features_file="../../../output/library_strategy-wf/scaled_prealn_feature_set.parquet",
                importance_file="../../../output/library_strategy-wf/random_forest_feature_importance.tsv",
                labels_file="../../../output/library_strategy-wf/summarized_metadata.parquet",
        ... )

        """

        self.update_figsize()

        self.features_file = features_file
        self.importance_file = importance_file
        self.labels_file = labels_file
        self.panel_kwargs = update_kwargs(PANEL_DEFAULTS, panel_kwargs)
        self.plot_kwargs = update_kwargs(PLOT_DEFAULTS, plot_kwargs)

        self.data = None
        self.panel = None

        self.load_data()
        self.plot()
        self.tweak()

    def load_data(self):
        labels = (
            pd.read_parquet(self.labels_file)
            .Fear_et_al_library_strategy.rename("library_strategy")
            .map(simplify_strategy)
        )

        important_features = (
            pd.read_table(self.importance_file, header=None, index_col=0).index[:6].tolist()
        )

        self.data = (
            pd.read_parquet(self.features_file)
            .reindex(columns=important_features)
            .join(labels)
            .reset_index()
            .melt(id_vars=["srx", "library_strategy"], var_name="feature", value_name="value")
        )

    def plot(self):
        self.panel = sns.FacetGrid(
            data=self.data, col="feature", hue="library_strategy", **self.panel_kwargs
        )
        self.panel.map(sns.kdeplot, "value", **self.plot_kwargs)

    def tweak(self):
        self.panel.set_titles("{col_name}")
        self.panel.axes[2].legend(loc="best")


    def savefig(self, file_name, **kwargs):
        self.panel.savefig(file_name, **kwargs)
