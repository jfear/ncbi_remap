"""Identify outliers in each library strategy and selection group.

Uses Isolation forest to find outliers in each group. Provides an outlier
score (close to -1 is outlier) and a flag (True is outlier).

"""
import sys

import numpy as np
import pandas as pd
import joblib

RANDOM_STATE = np.random.RandomState(42)

sys.path.insert(0, "../src")
from ncbi_remap.iforest import SraIsolationForest


def main():
    # Build RNA-Seq anomaly detection model
    rnaseq_features = pd.read_parquet(snakemake.input.rnaseq)
    iso = SraIsolationForest(
        rnaseq_features, random_state=RANDOM_STATE, iso_kwargs=dict(n_estimators=100)
    )
    joblib.dump(iso, snakemake.output.iso)

    # Check that proportion of outliers in Train and Test is similar.
    assert np.isclose(iso.prop_outliers_test, iso.prop_outliers_train, atol=0.01)

    # Save a list of good rnaseq samples
    inliers = iso.inliers(rnaseq_features).index.tolist()
    joblib.dump(inliers, snakemake.output.rnaseq_inliers)

    # Test model by mixing in other types of library strategy
    pd.concat(
        [
            mixin_contamination(iso, "est"),
            mixin_contamination(iso, "wgs"),
            mixin_contamination(iso, "chip"),
        ]
    ).to_parquet(snakemake.output.mixin)


def mixin_contamination(model, source: str):
    """Mix in different proportions of 'Other' samples.

    Mixes different proportions (0 - 1) of `other` and measures proportion of
    outliers.
    """
    name_mapper = {"est": "EST", "wgs": "WGS", "chip": "ChIP-Seq"}
    rnaseq = model.X_test
    other = pd.read_parquet(snakemake.input[source])
    n = min(rnaseq.shape[0], other.shape[0])

    res = []
    for rna_frac in [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]:
        other_frac = 1 - rna_frac
        mixture = pd.concat(
            [rnaseq.sample(n=int(n * rna_frac)), other.sample(n=int(n * other_frac))]
        )
        prop_outliers = model.prop_outliers(mixture)
        res.append((other_frac, prop_outliers, name_mapper[source]))

    return pd.DataFrame(
        res, columns=["Prop Contamination", "Prop Outliers", "Contamination Source"]
    )


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
