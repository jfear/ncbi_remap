"""Use RandomForest to classify library strategy and selection.

Here I use a MultiOutput Random Forest Classifier to simultaneously classify
library strategy and library selection based on features from the
pre-alignment workflow. I ignore OTHER samples during training, but include
them in the final classification.

"""
import os
from collections import namedtuple
from typing import List

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split, learning_curve
from sklearn.preprocessing import LabelEncoder
from sklearn.multioutput import MultiOutputClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report

np.random.seed(42)

Labels = namedtuple("Labels", "y_enc labeled_srxs strategy_encoder selection_encoder")


def main():
    # Import Features and Labels
    features = pd.read_parquet(snakemake.input["features"]).fillna(0)
    Y = read_and_encode_labels(features.index)
    X = features[features.index.isin(Y.labeled_srxs)].values

    # Split into training and test data
    # I was tyring to keep class proportions similar (stratify), but some
    # classes have so few members it is not possible.
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y.y_enc)

    # Train library strategy/selection random forest classifier
    forest = RandomForestClassifier(n_estimators=1000, oob_score=True)
    multi_target_forest = MultiOutputClassifier(forest, n_jobs=snakemake.threads)
    multi_target_forest.fit(X_train, Y_train)

    # Test classifier performance
    overall_score = multi_target_forest.score(X_test, Y_test)
    Y_pred = multi_target_forest.predict(X_test)

    with open(snakemake.output["metrics"], "w") as fh:
        fh.write(f"Overall Score: {overall_score}")
        fh.write(f"Library Strategy OOB: {multi_target_forest.estimators_[0].oob_score_}")
        fh.write(
            classification_report(
                Y_test[:, 0], Y_pred[:, 0], target_names=Y.strategy_encoder.classes_
            )
        )
        fh.write(f"{'-' * 80}")
        fh.write(f"Library Selection OOB: {multi_target_forest.estimators_[1].oob_score_}")
        fh.write(
            classification_report(
                Y_test[:, 1], Y_pred[:, 1], target_names=Y.selection_encoder.classes_
            )
        )

    # Save out feature importance
    feature_importance = (
        pd.concat(
            [
                pd.Series(
                    multi_target_forest.estimators_[0].feature_importances_,
                    name="Importance (Library Strategy)",
                    index=features.columns,
                ),
                pd.Series(
                    multi_target_forest.estimators_[1].feature_importances_,
                    name="Importance (Library Selection)",
                    index=features.columns,
                ),
            ],
            sort=False,
            axis=1,
        )
        .rename_axis("Feature")
        .sort_values("Importance (Library Strategy)", ascending=False)
    )

    feature_importance.to_csv(snakemake.output["feature_importance"], sep="\t", header=False)

    # Make predictions on entire dataset including the OTHERs
    Y_pred_all = multi_target_forest.predict(features.values)
    Y_pred_all_decode = pd.DataFrame(
        np.vstack(
            [
                Y.strategy_encoder.inverse_transform(Y_pred_all[:, 0]),
                Y.selection_encoder.inverse_transform(Y_pred_all[:, 1]),
            ]
        ).T,
        columns=["pred_library_strategy", "pred_library_selection"],
        index=features.index,
    )

    Y_pred_all_decode.to_csv(snakemake.output["predicted_labels"])


def read_and_encode_labels(srxs) -> Labels:
    Y_all = (
        pd.read_parquet(snakemake.input["labels"])
        .reindex(srxs)
        .replace({"other": "OTHER"})
        # Set classes with fewer than 40 samples to OTHER
        .apply(drop_rare_classes, axis=0)
    )

    # Remove OTHER because they are not informative
    Y = Y_all[~(Y_all == "OTHER").any(axis=1)]

    # Encode Library Strategy
    strategy_encoder = LabelEncoder()
    Y_strategy_enc = strategy_encoder.fit_transform(Y.library_strategy)

    # Encode Library strategy
    selection_encoder = LabelEncoder()
    Y_selection_enc = selection_encoder.fit_transform(Y.library_selection)

    # Put everything together
    Y_enc = np.vstack([Y_strategy_enc, Y_selection_enc]).T

    return Labels(Y_enc, Y.index, strategy_encoder, selection_encoder)


def drop_rare_classes(labels: pd.Series, threshold=40) -> pd.Series:
    label_counts = labels.value_counts()
    labels_to_drop = label_counts.index[label_counts < threshold].tolist()
    return labels.replace({k: "OTHER" for k in labels_to_drop})


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                features="../../output/library_strategy-wf/prealn_feature_set.parquet",
                labels="../../output/library_strategy-wf/sra_strategy_selection.parquet",
            ),
            threads=10
        )

    main()
