"""Identify outliers in each library strategy and selection group.

Uses Isolation forest to find outliers in each group. Provides an outlier
score (close to -1 is outlier) and a flag (True is outlier).

"""
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import IsolationForest
import shap
import joblib

RANDOM_STATE = np.random.RandomState(42)


def main():
    # Build RNA-Seq anomaly detection model
    rnaseq_features = pd.read_parquet(snakemake.input.rnaseq)
    X_train, X_test = train_test_split(rnaseq_features, random_state=RANDOM_STATE)
    iso = IsolationForest(n_estimators=100, random_state=RANDOM_STATE).fit(X_train)

    joblib.dump(iso, snakemake.output.iso)
    joblib.dump(X_train, snakemake.output.train)
    joblib.dump(X_test, snakemake.output.test)

    # Detect outliers in Train and Test and make sure similar proportions.
    pred_train = np.mean(np.where(iso.predict(X_train) == -1, True, False))
    pred_test = np.mean(np.where(iso.predict(X_test) == -1, True, False))
    assert np.isclose(pred_test, pred_train, atol=0.01)

    # Test model by mixing in other types of library strategy
    est_features = pd.read_parquet(snakemake.input.est)
    wgs_features = pd.read_parquet(snakemake.input.wgs)
    chip_features = pd.read_parquet(snakemake.input.chip)
    pd.concat(
        [
            mixin_contamination(iso, X_test, est_features, "EST"),
            mixin_contamination(iso, X_test, wgs_features, "WGS"),
            mixin_contamination(iso, X_test, chip_features, "ChIP-Seq"),
        ]
    ).to_parquet(snakemake.output.mixin)

    # Detect rnaseq outliers and save out
    outliers = df_outliers(iso, rnaseq_features)
    outliers.to_parquet(snakemake.output.outliers)

    # Save a list of good rnaseq samples
    inliers = outliers[~outliers.squeeze()].index.tolist()
    joblib.dump(inliers, snakemake.output.rnaseq_inliers)


def mixin_contamination(model, rnaseq: pd.DataFrame, other: pd.DataFrame, source: str):
    """Mix in different proportions of 'Other' samples.

    Mixes different proportions (0 - 1) of `other` and measures proportion of
    outliers.
    """
    n = min(rnaseq.shape[0], other.shape[0])
    res = []
    for rna_frac in [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]:
        other_frac = 1 - rna_frac
        mixture = pd.concat(
            [rnaseq.sample(n=int(n * rna_frac)), other.sample(n=int(n * other_frac))]
        )
        prediction = model.predict(mixture)
        prop_outliers = np.mean(np.where(prediction == -1, 1, 0))
        res.append((other_frac, prop_outliers, source))

    return pd.DataFrame(
        res, columns=["Prop Contamination", "Prop Outliers", "Contamination Source"]
    )


def df_outliers(model, features: pd.DataFrame) -> pd.DataFrame:
    outliers = pd.DataFrame(
        np.where(model.predict(features) == -1, True, False),
        index=features.index,
        columns=["rnaseq_outliers"],
    )
    return outliers


if __name__ == "__main__":
    if "snakemake" not in locals() or not hasattr(snakemake, "scriptdir"):
        from ncbi_remap.mock import MockSnake

        snakemake = MockSnake(
            input=dict(
                rnaseq="../../output/library_strategy-wf/rnaseq_features.parquet",
                est="../../output/library_strategy-wf/est_features.parquet",
                wgs="../../output/library_strategy-wf/wgs_features.parquet",
                chip="../../output/library_strategy-wf/chip_features.parquet",
            ),
        )

    main()
