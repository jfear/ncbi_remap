import pandas as pd

def load_features(library_strategy):
    selected_srxs = pd.read_parquet(snakemake.input.labels).query("library_strategy == library_strategy").index
    return pd.read_parquet(snakemake.input.features).reindex(selected_srxs).dropna(how="all")


def groupby_classify(labels, features):
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