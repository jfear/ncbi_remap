"""Create a list of SRXs that are annotated as RNA-Seq in the SRA"""
from pymongo import MongoClient


def main():
    rnaseq_srxs = sorted(
        list(
            {
                record["srx"]
                for record in NCBI.aggregate(
                    [
                        {"$match": {"library_strategy": "RNA-Seq"}},
                        {"$project": {"_id": 0, "srx": 1}},
                    ]
                )
            }
        )
    )

    with open(snakemake.output[0], "w") as fh:
        fh.write("\n".join(rnaseq_srxs))


if __name__ == "__main__":
    if "snakemake" not in locals() or not hasattr(snakemake, "scriptdir"):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(output=dict("output/rnaseq_srxs.pkl"))

    try:
        CLIENT = MongoClient()
        DB = CLIENT["sramongo"]
        NCBI = DB["ncbi"]
        main()
    finally:
        CLIENT.close()
