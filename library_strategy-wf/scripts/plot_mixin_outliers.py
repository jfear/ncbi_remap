"""Plot Mixin Experiment Results"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, "../src")
from ncbi_remap.plotting import style_use


def main():
    style_use(snakemake.params.get("style", "sra"))

    mixin = pd.read_parquet(snakemake.input[0])
    sns.lmplot(
        "Prop Contamination",
        "Prop Outliers",
        hue="Contamination Source",
        data=mixin,
        palette=["C1", "C2", "C4"],
        hue_order=["EST", "WGS", "ChIP-Seq"],
        legend=False,
    )
    plt.legend(loc="upper left")
    plt.savefig(snakemake.output[0])


if __name__ == "__main__":
    if "snakemake" not in locals() or not hasattr(snakemake, "scriptdir"):
        from ncbi_remap.mock import MockSnake

        snakemake = MockSnake(input="../../output/library_strategy-wf/contamination_mixin.parquet")

    main()
