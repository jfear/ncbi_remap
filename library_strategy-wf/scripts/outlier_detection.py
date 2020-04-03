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
    labels_from_sra = pd.read_parquet(snakemake.input.labels)
    strategy2srx = {k: v.index for k, v in labels_from_sra.groupby("library_strategy")}

    features = pd.read_parquet(snakemake.input.features)

    rnaseq_features = features.reindex(strategy2srx["RNA-Seq"]).dropna(how="all")
    rnaseq_train, rnaseq_test = train_test_split(rnaseq_features, random_state=RANDOM_STATE)
    model = IsolationForest(n_estimators=100, random_state=RANDOM_STATE).fit(rnaseq_train)

    # Detect outliers and save out
    estimate_feature_importance(model, rnaseq_features)
    outliers = save_outliers(model, rnaseq_features)
    save_list_of_good_rnaseq(outliers)

    # Validate model using ChIP-Seq and WGS mixin
    chip_features = features.reindex(strategy2srx["ChIP-Seq"]).dropna(how="all")
    wgs_features = features.reindex(strategy2srx["WGS"]).dropna(how="all")
    pd.concat(
        [
            mixin_contamination(rnaseq_test, chip_features, "ChIP-Seq", model),
            mixin_contamination(rnaseq_test, wgs_features, "WGS", model),
        ]
    ).to_parquet(snakemake.output.mixin)


def estimate_feature_importance(model, rnaseq_features):
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(rnaseq_features)
    joblib.dump(explainer, snakemake.output.shap_model)
    joblib.dump(rnaseq_features, snakemake.output.shap_features)
    joblib.dump(shap_values, snakemake.output.shap_values)


def save_outliers(model, rnaseq_features):
    outliers = pd.DataFrame(
        np.where(model.predict(rnaseq_features) == -1, True, False),
        index=rnaseq_features.index,
        columns=["rnaseq_outliers"],
    )
    outliers.to_parquet(snakemake.output.outliers)
    return outliers


def save_list_of_good_rnaseq(outliers):
    rnaseq = outliers[~outliers.squeeze()].index.tolist()
    joblib.dump(rnaseq, snakemake.output.rnaseq_samples)


def mixin_contamination(rnaseq: pd.DataFrame, other: pd.DataFrame, source, model):
    n = min(rnaseq.shape[0], other.shape[0])
    res = []
    for rna_frac in np.arange(0.1, 1.1, step=0.1)[::-1]:
        _rna_frac = np.round(rna_frac, 2)
        _other_frac = 1 - _rna_frac
        mixture = pd.concat(
            [rnaseq.sample(n=int(n * _rna_frac)), other.sample(n=int(n * _other_frac))]
        )
        prediction = model.predict(mixture)
        prop_outliers = np.mean(np.where(prediction == -1, 1, 0))
        res.append((_other_frac, prop_outliers, source))

    return pd.DataFrame(
        res, columns=["Prop Contamination", "Prop Outliers", "Contamination Source"]
    )


if __name__ == "__main__":
    if "snakemake" not in locals() or not hasattr(snakemake, "scriptdir"):
        from ncbi_remap.mock import MockSnake

        snakemake = MockSnake(
            input=dict(
                features="../../output/library_strategy-wf/scaled_prealn_feature_set.parquet",
                labels="../../output/library_strategy-wf/sra_strategy_selection.parquet",
            ),
            output=dict(shap_model="", shap_values="", outliers="", mixin=""),
        )

    main()
