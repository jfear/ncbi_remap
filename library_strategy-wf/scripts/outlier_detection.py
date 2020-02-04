import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.ensemble import IsolationForest
from sklearn.preprocessing import StandardScaler

np.random.seed(42)
plt.rcParams["figure.figsize"] = (5, 5)
PLOT_DEFAULTS = dict(x="UMAP1", y="UMAP2", s=12, rasterized=True, linewidth=0.2)


def main():
    features = pd.read_parquet(snakemake.input.features).dropna()
    scaler = StandardScaler().fit(features)
    srxs = features.index

    labels = pd.read_parquet(snakemake.input.labels).reindex(srxs)
    strategy = pd.get_dummies(labels.library_strategy)
    selection = pd.get_dummies(labels.library_selection)

    umap = pd.read_table(snakemake.input.umap, index_col=0).join(labels, how="inner")

    # Find outliers in RNA-Seq subset
    rnaseq_srxs = labels.query("library_strategy == 'RNA-Seq'").index
    rnaseq_features = features.reindex(rnaseq_srxs)
    rnaseq_umap = umap.reindex(rnaseq_srxs)

    X = scaler.transform(rnaseq_features)
    n_samples, n_features = X.shape
    max_features = int(np.ceil(np.sqrt(n_features)))

    clf = IsolationForest(max_features=max_features, bootstrap=True).fit(X)
    sample_scores = clf.score_samples(X)
    flag_outlier = np.where(clf.predict(X) == -1, True, False)
    n_outliers = flag_outlier.sum()
    rnaseq_umap["flag_outlier"] = flag_outlier

    sns.scatterplot(hue="flag_outlier", data=rnaseq_umap, **PLOT_DEFAULTS)

    # Where does this map on entire dataset
    umap["flag_outlier"] = False
    umap.loc[rnaseq_umap[rnaseq_umap.flag_outlier].index, "flag_outlier"] = True

    sns.scatterplot(hue="flag_outlier", data=umap, **PLOT_DEFAULTS)



if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                features="../../output/library_strategy-wf/prealn_feature_set.parquet",
                umap="../../output/library_strategy-wf/umap_prealn_features_embeddings.tsv",
                labels="../../output/library_strategy-wf/sra_strategy_selection.parquet",
            )
        )

    main()
