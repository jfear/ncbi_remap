"""Identify outliers in each library strategy and selection group.

Uses Isolation forest to find outliers in each group. Provides an outlier
score (close to -1 is outlier) and a flag (True is outlier).

"""
import os
from typing import Generator

import numpy as np
import pandas as pd
from sklearn.ensemble import IsolationForest

# Setting the numpy seed was not working, so I am directly using random_state in class.
RANDOM_STATE = 42


def main():
    scaled_features = pd.read_parquet(snakemake.input.features)
    labels = pd.read_parquet(snakemake.input.labels).reindex(scaled_features.index)

    library_strategy = pd.concat(groupby_classify(labels.library_strategy, scaled_features))
    library_selection = pd.concat(groupby_classify(labels.library_selection, scaled_features))
    df = library_strategy.join(library_selection)

    df.to_parquet(snakemake.output[0])


def groupby_classify(labels, features) -> Generator[pd.DataFrame, None, None]:
    """Find outliers for each group separately.
    
    Parameters
    ----------
    labels : pd.Series
        A Series of either library strategy or selection.
    features : pd.DataFrame
        A DataFrame with scaled features.
    """
    type_ = labels.name
    for label, data in labels.groupby(labels):
        srxs = data.index
        feature_subset = features.reindex(srxs)

        X = feature_subset.values
        n_samples, n_features = X.shape
        max_features = int(np.ceil(np.sqrt(n_features)))

        clf = IsolationForest(
            max_features=max_features, bootstrap=True, random_state=RANDOM_STATE
        ).fit(X)
        sample_scores = clf.score_samples(X)
        flag_outlier = np.where(clf.predict(X) == -1, True, False)

        df = pd.DataFrame(
            index=srxs, columns=[type_, f"{type_}_outlier_score", f"{type_}_flag_outlier"]
        )
        df[type_] = label
        df[f"{type_}_outlier_score"] = sample_scores
        df[f"{type_}_flag_outlier"] = flag_outlier

        yield df


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                features="../../output/library_strategy-wf/scaled_prealn_feature_set.parquet",
                labels="../../output/library_strategy-wf/sra_strategy_selection.parquet",
            )
        )

    main()
