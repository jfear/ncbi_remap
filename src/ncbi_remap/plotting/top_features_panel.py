from typing import Union, Tuple, Optional

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import NcbiPlotter, update_kwargs

STRATEGY_ORDER = ["Other", "RNA-Seq", "ChIP-Seq"]
PANEL_DEFAULTS = dict(
    sharex=False,
    sharey=False,
    hue_order=STRATEGY_ORDER,
    palette=["lightgray", "C0", "C1"],
    col_wrap=3,
    aspect=4/3
)
PLOT_DEFAULTS = dict()
FEATURE_NAME_MAPPER = {
    "PCT_UTR_BASES": "% UTR",
    "PCT_INTRONIC_BASES": "% Intronic",
    "PCT_MRNA_BASES": "% mRNA",
    "PCT_INTERGENIC_BASES": "% Intergenic",
    "PCT_CODING_BASES": "% Coding",
    "pos_0": "First Base Coverage",
}


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
        n_features: int = 6,
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
        self.n_features = n_features

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
            pd.read_table(self.importance_file, header=None, index_col=0).index[:self.n_features].tolist()
        )

        self.data = (
            pd.read_parquet(self.features_file)
            .reindex(columns=important_features)
            .join(labels)
            .reset_index()
            .melt(id_vars=["srx", "library_strategy"], var_name="feature", value_name="value")
        )

        self.data.feature = self.data.feature.map(FEATURE_NAME_MAPPER)

    def plot(self):
        self.panel = sns.FacetGrid(
            data=self.data, col="feature", hue="library_strategy", height=self.fig_height, **self.panel_kwargs
        )
        self.panel.map(sns.kdeplot, "value", **self.plot_kwargs)

    def tweak(self):
        self.panel.set_titles("{col_name}", va="top")
        self.panel.set_xlabels("")
        self.panel.set_xticklabels([])
        self.panel.set_yticklabels([])
        self.panel.axes[2].legend(loc="upper left")
        plt.subplots_adjust(hspace=0.18, wspace=0.05)


    def savefig(self, file_name, **kwargs):
        self.panel.savefig(file_name, **kwargs)
