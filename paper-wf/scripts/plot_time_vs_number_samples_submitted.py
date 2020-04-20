"""Plot time vs number samples submitted."""
import sys

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pymongo import MongoClient

sys.path.insert(0, "../src")
from ncbi_remap.plotting import style_use


def main():
    # Data prep
    df = get_submission_dates()

    submission_by_year = (
        df.resample("Y")
        .size()
        .rename("num_samples")
        .to_frame()
        .assign(cum_num_samples=lambda x: x.num_samples.cumsum())
        .assign(year=lambda x: x.index.year.astype(str))
        # had to use this for lmplot
        .assign(idx=lambda x: range(x.shape[0]))
    )

    # Make plot
    style_use(snakemake.params.get("style", "sra"))
    plt.rcParams["figure.figsize"] = (2.636, 1.496)

    _, ax = plt.subplots()  # type: plt.Figure, plt.Axes
    ax2 = ax.twinx()

    sns.lineplot("year", "cum_num_samples", data=submission_by_year, ax=ax)
    sns.regplot("idx", "num_samples", data=submission_by_year, lowess=True, ax=ax2, color="C1")

    # Add fill
    ax.fill_between(submission_by_year.year, submission_by_year.cum_num_samples, color="C0")

    # Clean up X-axis
    ax.xaxis.set_major_locator(plt.MaxNLocator(4))
    ax.set_xlabel("Submission Date")

    # Clean up Y-Axis 1
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, pos: f"{y/1000:.0f}"))
    ax.set_ylabel("Cumulative # Samples\n(Thousands)", color="C0")

    # Clean up Y-Axis 2
    ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, pos: f"{y/1000:.0f}"))
    ax2.set_ylabel("# Samples\n(Thousands)", color="C1")

    # Remove spines
    sns.despine(ax=ax, left=True, bottom=True)
    sns.despine(ax=ax2, left=True, bottom=True)

    plt.savefig(snakemake.output[0])


def get_submission_dates():
    try:
        client = MongoClient()
        db = client["sramongo"]
        ncbi = db["ncbi"]

        df = pd.DataFrame(
            ncbi.aggregate(
                [{"$project": {"_id": 0, "srx": 1, "submission_date": "$sra_create_date"}},]
            )
        )

        df["submission_date"] = pd.to_datetime(df.submission_date)

    finally:
        client.close()

    return df.set_index("submission_date")


if __name__ == "__main__":

    main()
