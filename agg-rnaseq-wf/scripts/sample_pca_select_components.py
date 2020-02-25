import os

import numpy as np
import pandas as pd
from joblib import load
import matplotlib.pyplot as plt

from ncbi_remap import plotting


def main():
    plot_explained_varaince()

    # Keep selected dimensions
    df = pd.read_parquet(
        snakemake.input.reduced_values,
        columns=[f"PC{x + 1}" for x in range(snakemake.params.n_comp)],
    ).rename_axis("srx")

    df.to_parquet(snakemake.output.selected_values)


def plot_explained_varaince():
    reducer = load(snakemake.input.model)  # type: sklearn.decomposition.PCA

    # How many components to get 99% of variance
    cumsum_ = np.cumsum(np.sort(reducer.explained_variance_ratio_)[::-1])
    _99_threshold = np.sum(cumsum_ < 0.99)
    _90_threshold = np.sum(cumsum_ < 0.90)
    _85_threshold = np.sum(cumsum_ < 0.85)
    _80_threshold = np.sum(cumsum_ < 0.80)
    _selected_threshold = snakemake.params.n_comp
    _selected_threshold_variance = cumsum_[_selected_threshold]

    # Plot of explained variance
    plt.style.use("sra")
    fig, ax = plt.subplots()
    ax.plot(cumsum_)
    ax.axvline(_selected_threshold, ls="--", color="k")
    ax.text(
        _99_threshold + 1_000,
        0.8,
        "\n".join(
            [
                f"99% Threshod ({_99_threshold:,})",
                f"90% Threshod ({_90_threshold:,})",
                f"85% Threshod ({_85_threshold:,})",
                f"80% Threshod ({_80_threshold:,})",
                f"Selected Threshod ({_selected_threshold:,}): {_selected_threshold_variance:.3}",
            ]
        ),
    )
    plt.savefig(snakemake.output.variance)


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG"):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                model="../../output/agg-rnaseq-wf/sample_pca_reducer.pkl",
                reduced_values="../../output/agg-rnaseq-wf/sample_principal_components.parquet",
            ),
            params=dict(n_comp=4_000),
        )

    main()
