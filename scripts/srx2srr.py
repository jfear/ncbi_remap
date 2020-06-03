import os

import pandas as pd
from pymongo import MongoClient


def main(db):
    pd.DataFrame(
        db.aggregate(
            [
                {"$unwind": {"path": "$runs"}},
                {"$match": {"runs.nspots": {"$ne": "", "$gt": 100000}}},
                {"$project": {"_id": 0, "srx": 1, "srr": "$runs.srr", "libsize": "$runs.nspots"}},
            ]
        )
    ).sort_values(["srx", "srr"]).to_csv(snakemake.output[0], index=False)


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(output=dict("output/srx2srr.csv"))

    try:
        client = MongoClient()
        db = client["sramongo"]["ncbi"]
        main(db)
    finally:
        client.close()
