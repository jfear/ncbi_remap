import os

import pandas as pd
from pymongo import MongoClient


def main():
    pd.DataFrame(
        NCBI.aggregate(
            [{"$unwind": {"path": "$runs"}}, {"$project": {"_id": 0, "srx": 1, "srr": "$runs.srr"}}]
        )
    ).sort_values(["srx", "srr"]).to_csv(snakemake.output[0], index=False)


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(output=dict("output/srx2srr.csv"))

    try:
        CLIENT = MongoClient()
        DB = CLIENT["sramongo"]
        NCBI = DB["ncbi"]
        main()
    finally:
        CLIENT.close()
