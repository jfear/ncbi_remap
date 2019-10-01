#!/usr/bin/env python
"""Interact with the pre-alignment workflow data store.

The HDF5 data store is used as a queue system as well as a data store for
various flags and datasets generated by the pre-alignment workflow. This script
provides easy access to the pre-alignment workflow queue as well as easy updating
of the queue and data. It also parses large files and saves them as a
data frame. Finally, it can scan fastq files and determine if they are
abi-solid format. Specifically this script can:

* Initialization
    * Initialize the data store from a mongo database `./prealn-store.py init`
    * Append new IDs from mongo database to data store `./prealn-store.py init --append`
* Interact with queue
    * Print a summary of the queue `./prealn-store.py queue --print`
    * Scan the `samples` directory and remove completed items from the queue
      `./prealn-store.py queue --update -j 8`
    * Check for abi-solid fastqs and remove them from the queue
      `./prealn-store.py queue --abi -j 8`
* Scan the `samples` directory and parse data and add small datasets and flags
  to the data sore and save large datasets as a data frame `./prealn-store.py data
  -j 8`

Flags checked and saved:

* LAYOUT
* STRAND
* download_bad
* alignment_bad
* quality_scores_bad
* abi_solid

Small datasets added to the data store:

* fastq summary counts
* `fastq_screen` results
* `hisat2` alignment log
* `samtools stats`
* `bamtools stats`
* `picard markduplicates`
* `picard collectRNASeqMetrics` unstranded
* `picard collectRNASeqMetrics` first stranded
* `picard collectRNASeqMetrics` second stranded
* Genic `feature counts` summary

Large datasets parsed to data frame and saved:

* `samtools idxstats`
* Genic `feature counts` counts
* Genic `feature counts` junction counts
"""

from pathlib import Path
import argparse

import numpy as np
import pandas as pd
from dask import delayed
from dask.distributed import Client

from ncbi_remap.logging import logger
from ncbi_remap.snakemake import get_patterns
from ncbi_remap.queue import dask_run_srr_checker, dask_run_srr_parser, check_indicator_file

WORKDIR = Path(__file__).absolute().parent
DATA_STORE = Path(WORKDIR, "../output/sra.h5").as_posix()


def get_options():
    parser = argparse.ArgumentParser(description=__doc__, usage="prealn_store.py <command> [args]")

    subparsers = parser.add_subparsers(dest="command", help="Commands")

    parser_init = subparsers.add_parser("init", help="Initialize data store")
    parser_init.add_argument(
        "--host",
        dest="host",
        action="store",
        default="localhost",
        help="Host running mongo database [Optional].",
    )
    parser_init.add_argument(
        "--port",
        dest="port",
        action="store",
        type=int,
        default=27017,
        help="Mongo database port [Optional].",
    )
    parser_init.add_argument(
        "-j",
        dest="n_workers",
        action="store",
        default=2,
        type=int,
        help="Number of workers [default 2]",
    )
    parser_init.add_argument(
        "--append",
        dest="append",
        action="store_true",
        help="Download and append new ids to the data store.",
    )

    parser_queue = subparsers.add_parser("queue", help="Interact with prealn/queue")
    parser_queue.add_argument(
        "-j",
        dest="n_workers",
        action="store",
        default=2,
        type=int,
        help="Number of workers [default 2]",
    )
    parser_queue.add_argument(
        "--update", dest="queue_update", action="store_true", help="Run update on prealn/queue."
    )
    parser_queue.add_argument(
        "--print", dest="queue_print", action="store_true", help="Print the prealn/queue"
    )
    parser_queue.add_argument(
        "--abi", dest="queue_abi", action="store_true", help="Check for ABI Solid data by fastq."
    )

    parser_data = subparsers.add_parser("data", help="Interact with prealn/workflow")
    parser_data.add_argument(
        "-j",
        dest="n_workers",
        action="store",
        default=2,
        type=int,
        help="Number of workers [default 2]",
    )
    return parser.parse_args()


def get_keepers():
    return [
        patterns["fastq_screen"],
        patterns["layout"],
        patterns["strand"],
        patterns["hisat2"]["summary"],
        patterns["feature_counts"]["summary"],
        patterns["samtools_stats"],
        patterns["samtools_idxstats"],
        patterns["bamtools_stats"],
        patterns["picard"]["markduplicates"]["metrics"],
    ]


def get_updated_db_ids():
    """Pull out SRX, SRR from mongo."""
    from datetime import datetime
    from pymongo import MongoClient

    mongoClient = MongoClient(host=args.host, port=args.port)
    db = mongoClient["sramongo"]
    ncbi = db["ncbi"]

    # Get ids from database with > 1k reads and ≥ 25bp read len
    df = pd.DataFrame(
        ncbi.aggregate(
            [
                {"$match": {"sra_create_date": {"$lte": datetime(2019, 3, 31)}}},
                {"$unwind": {"path": "$runs"}},
                {"$match": {"runs.srr": {"$exists": True}}},
                {
                    "$project": {
                        "_id": 0,
                        "srx": 1,
                        "srr": "$runs.srr",
                        "nspots": {
                            "$cond": [
                                {"$eq": ["$runs.nspots", ""]},
                                0,
                                {"$ifNull": ["$runs.nspots", 0]},
                            ]
                        },
                        "read_count_r1": {"$ifNull": ["$runs.read_count_r1", 0]},
                        "read_count_r2": {"$ifNull": ["$runs.read_count_r2", 0]},
                        "read_len_r1": {"$ifNull": ["$runs.read_len_r1", 0]},
                        "read_len_r2": {"$ifNull": ["$runs.read_len_r2", 0]},
                    }
                },
                {
                    # Keep samples with > 1,000 reads
                    "$match": {
                        "$or": [
                            {
                                # Keep samples with all 0 b/c we just don't know
                                "$and": [
                                    {"nspots": {"$eq": 0}},
                                    {"read_count_r1": {"$eq": 0}},
                                    {"read_count_r2": {"$eq": 0}},
                                ]
                            },
                            {
                                "$or": [
                                    {"nspots": {"$gt": 1000}},
                                    {"read_count_r1": {"$gt": 1000}},
                                    {"read_count_r2": {"$gt": 1000}},
                                ]
                            },
                        ]
                    }
                },
                {
                    # Keep samples with > 25 bp reads lens
                    "$match": {
                        "$or": [
                            {
                                # Keep samples with all 0 b/c we just don't know
                                "$and": [{"read_len_r1": {"$eq": 0}}, {"read_len_r2": {"$eq": 0}}]
                            },
                            {"$or": [{"read_len_r1": {"$gte": 25}}, {"read_len_r2": {"$gte": 25}}]},
                        ]
                    }
                },
            ]
        )
    )

    return df[["srx", "srr"]]


def dask_run_large_data(name, pattern, munger, parser):
    ids = store["prealn/complete"]
    futures = []
    for idx, dat in ids.groupby("srx"):
        srx = idx
        srrs = dat.srr.tolist()
        futures.append(delayed(munger)(srx, srrs, name, pattern, parser))

    return daskClient.gather(daskClient.compute(futures))


def check_flag_file(srx, srr, pattern):
    """Parse flag file and return value."""
    fname = Path(pattern.format(srx=srx, srr=srr))
    try:
        with fname.open() as fh:
            flag = fh.read().strip()
        return srx, srr, flag
    except FileNotFoundError:
        return srx, srr, np.nan


def check_abi_solid(srx, srr, *args, **kwargs):
    """Check Fastq is Abi Solid."""
    logger.info("Checking queue for Abi Solid")
    import gzip

    layout = Path(patterns["layout"].format(srx=srx, srr=srr))
    abi_solid = Path(patterns["abi_solid"].format(srx=srx, srr=srr))
    fname = Path(patterns["fastq"]["r1"].format(srx=srx, srr=srr))
    r2 = Path(patterns["fastq"]["r2"].format(srx=srx, srr=srr))

    if abi_solid.exists():
        return

    try:
        with layout.open() as fh:
            flag = fh.read().strip()

        if flag == "keep_R2":
            fname = r2
        with gzip.open(fname, "rt") as fh:
            _, l = fh.readline(), fh.readline()
            if l.startswith("T"):
                abi_solid.touch()
    except (KeyError, FileNotFoundError):
        pass

    logger.info("Abi Solid Complete")


def check_outputs(srx, srr, *args):
    """Check if output file is present."""
    for pattern in get_keepers():
        fname = Path(pattern.format(srx=srx, srr=srr))
        if not fname.exists():
            return np.nan, np.nan
    return srx, srr


def process_outputs():
    """Check if all output files are present and update queues."""
    logger.info("Checking outputs.")
    res = dask_run_srr_checker(store["prealn/queue"], check_outputs, None, daskClient)
    ids = pd.DataFrame(res, columns=["srx", "srr"]).dropna()

    mask = store["prealn/queue"].srr.isin(ids.srr)
    cnt = mask.sum()
    if cnt > 0:
        store.remove("prealn/queue", mask)
        logger.info("{:,} SRRs removed from queue.".format(cnt))
        store.append("prealn/complete", ids, data_columns=True, format="t")
        store.append("aln/queue", ids, data_columns=True, format="t")
    logger.info("Complete")


def process_flags():
    """Check all flag and indicator files and store results."""
    queue = store["prealn/queue"]

    # flag files
    for key in ["layout", "strand"]:
        logger.info(f"Checking {key} files.")
        if store.__contains__(key):
            ids = queue[~queue.srr.isin(store[key].index.get_level_values("srr"))]
        else:
            ids = queue
        dat = dask_run_srr_checker(ids, check_flag_file, patterns[key], daskClient)
        df = (
            pd.DataFrame(dat, columns=["srx", "srr", key])
            .set_index(["srx", "srr"])
            .iloc[:, 0]
            .dropna()
        )

        store.append(key, df, data_columns=True, format="t")

    # indicator files
    for key in ["download_bad", "alignment_bad", "quality_scores_bad", "abi_solid"]:

        queue = store["prealn/queue"]

        logger.info(f"Checking {key} files.")
        if store.__contains__(key):
            ids = queue[~queue.srr.isin(store[key].srr)]
        else:
            ids = queue
        dat = dask_run_srr_checker(ids, check_indicator_file, patterns[key], daskClient)
        df = pd.DataFrame(dat, columns=["srx", "srr"]).dropna().reset_index(drop=True)

        # remove form queue because an indicator files indicates that
        # something is something wrong. Using mask because querying a list of
        # ids was not working.
        mask = queue.srr.isin(df.srr)
        cnt = mask.sum()
        if cnt > 0:
            store_key = "prealn/" + key
            store.append(store_key, df, data_columns=True, format="t")

            logger.info("{:,} SRRs removed from queue.".format(cnt))
            store.remove("prealn/queue", mask)


def munge_as_dataframe(srx, srrs, name, pattern, parser):
    oname = Path(name.format(srx=srx))

    if oname.exists():
        return

    dfs = []
    for srr in srrs:
        fname = pattern.format(srx=srx, srr=srr)
        df = parser(fname)
        dfs.append(df.assign(srx=srx, srr=srr))

    if len(dfs) > 1:
        df = pd.concat(dfs)
    else:
        df = dfs[0]

    oname.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(oname.as_posix(), engine="pyarrow")


def munge_big_data():
    """Munge larger datasets for easier handling with dask.

    Larger tables are not performant when added to the datastore. Dask does
    well accessing multiple files on disk in parallel. Dask documentation
    suggests to use the apache parquet. Here we munge the data and save it out
    to the parquet format.
    """
    from ncbi_remap.parser import (
        parse_samtools_idxstats,
        parse_featureCounts_counts,
        parse_featureCounts_jcounts,
    )

    logger.info(f"Building samtools idxstats files")
    dask_run_large_data(
        Path(WORKDIR, "../output/prealn-wf/samtools_idx_stats/{srx}.parquet").as_posix(),
        patterns["samtools_idxstats"],
        munge_as_dataframe,
        parse_samtools_idxstats,
    )

    logger.info(f"Building counts files")
    dask_run_large_data(
        Path(WORKDIR, "../output/prealn-wf/gene_counts/{srx}.parquet").as_posix(),
        patterns["feature_counts"]["counts"],
        munge_as_dataframe,
        parse_featureCounts_counts,
    )

    logger.info(f"Building junction counts files")
    dask_run_large_data(
        Path(WORKDIR, "../output/prealn-wf/junction_counts/{srx}.parquet").as_posix(),
        patterns["feature_counts"]["jcounts"],
        munge_as_dataframe,
        parse_featureCounts_jcounts,
    )


def agg_small(name, pattern, parser):
    complete = store["prealn/complete"]
    key = "prealn/workflow/" + name
    logger.info(f"Checking {key} files.")

    if store.__contains__(key):
        srrs_not_processed = ~complete.srr.isin(store[key].index.get_level_values("srr"))
        toRun = complete[srrs_not_processed]
    else:
        toRun = complete

    dat = dask_run_srr_parser(toRun, parser, pattern, daskClient)
    dfs = [x for x in dat if x is not None]
    if len(dfs) > 0:
        df = pd.concat(dfs)
        try:
            store.append(key, df, data_columns=True, format="t")
        except ValueError:
            # Sometimes data types don't match, use pandas to coerce.
            tmp = store[key]
            store.append(key, pd.concat([tmp, df]), data_columns=True, format="t", append=False)


def add_small_data_to_store():
    """Parse smaller datasets and add them to the store.

    Smaller datasets are munged and added to their corresponding table in the
    data store.
    """
    from ncbi_remap.parser import (
        parse_fastq_summary,
        parse_fastq_screen,
        parse_hisat2,
        parse_samtools_stats,
        parse_bamtools_stats,
        parse_picard_markduplicate_metrics,
        parse_picardCollect_summary,
        parse_picardCollect_hist,
        parse_featureCounts_summary,
    )

    agg_small("fastq", patterns["fastq"]["summary"], parse_fastq_summary)

    agg_small("fastq_screen", patterns["fastq_screen"], parse_fastq_screen)

    agg_small("hisat2", patterns["hisat2"]["bam"] + ".log", parse_hisat2)

    agg_small("samtools_stats", patterns["samtools_stats"], parse_samtools_stats)

    agg_small("bamtools_stats", patterns["bamtools_stats"], parse_bamtools_stats)

    agg_small(
        "markduplicates",
        patterns["picard"]["markduplicates"]["metrics"],
        parse_picard_markduplicate_metrics,
    )

    agg_small(
        "collectrnaseqmetrics/first",
        patterns["picard"]["collectrnaseqmetrics"]["metrics"]["first"],
        parse_picardCollect_summary,
    )

    agg_small(
        "collectrnaseqmetrics/second",
        patterns["picard"]["collectrnaseqmetrics"]["metrics"]["second"],
        parse_picardCollect_summary,
    )

    agg_small(
        "collectrnaseqmetrics/unstranded",
        patterns["picard"]["collectrnaseqmetrics"]["metrics"]["unstranded"],
        parse_picardCollect_summary,
    )

    agg_small(
        "collectrnaseqmetrics/genebody",
        patterns["picard"]["collectrnaseqmetrics"]["metrics"]["unstranded"],
        parse_picardCollect_hist,
    )

    agg_small(
        "feature_counts/summary", patterns["feature_counts"]["summary"], parse_featureCounts_summary
    )


def initialize_store():
    """Initialize the data store with ids and prealn/queue."""
    logger.info("Initialize data store.")
    logger.info("Querying database for ids.")
    store["ids"] = get_updated_db_ids()
    logger.info("Adding ids to store.")
    store.put("prealn/queue", store["ids"], data_columns=True, format="table")
    logger.info("Initialization complete")


def append_store():
    """Download new ids and append to store."""
    logger.info("Append new ids to data store.")
    curr_ids = store["ids"].copy()

    logger.info("Querying database for ids.")
    ids = get_updated_db_ids()
    store["ids"] = ids

    # Sometimes an id is removed from the SRA, need to make sure to remove it
    # from the queue.
    srr_no_longer_in_store = ~curr_ids.srr.isin(ids.srr)
    removed_ids = curr_ids[srr_no_longer_in_store]

    if len(removed_ids) > 0:
        logger.info("There are {:,} ids no longer in the SRA.".format(len(removed_ids)))
        srrs = removed_ids.srr
        store.remove("prealn/queue", "srr == srrs")

    # Only add ids not already in the store.
    srr_not_in_store = ~ids.srr.isin(curr_ids.srr)
    new_ids = ids[srr_not_in_store]

    if len(new_ids) == 0:
        logger.info("There are no new ids.")
        return

    logger.info("There are {:,} new ids.".format(len(new_ids)))
    logger.info("Adding ids to store.")
    store.append("prealn/queue", new_ids, data_columns=True, format="table")
    logger.info("Initialization complete")


def print_queue():
    def parse_pairs(name, key):
        try:
            return name, store[key].shape[0]
        except KeyError:
            return name, 0

    pairs = [
        ("ids in the system", "ids"),
        ("queued", "prealn/queue"),
        ("completed", "prealn/complete"),
        ("download bad", "prealn/download_bad"),
        ("quality scores bad", "prealn/quality_scores_bad"),
        ("alignment bad", "prealn/alignment_bad"),
        ("abi solid", "prealn/abi_solid"),
    ]
    parsed = [parse_pairs(k, v) for k, v in pairs]

    if store.__contains__("layout"):
        layout = store["layout"].value_counts()
        parsed.extend(
            [
                ("Single End", layout.SE),
                ("Pair End", layout.PE),
                ("Really Single End R1", layout.keep_R1),
                ("Really Single End R2", layout.keep_R2),
            ]
        )

    if store.__contains__("strand"):
        strand = store["strand"].value_counts()

        parsed.extend(
            [
                ("first strand", strand.same_strand),
                ("second strand", strand.opposite_strand),
                ("unstranded", strand.unstranded),
            ]
        )

    report = "\nCurrent Queue Summary\n"
    for k, v in parsed:
        report += "{:,}\t\t\t{}\n".format(v, k)

    print(report)


def update_queue():
    logger.info("Updating queue.")
    process_flags()
    process_outputs()
    logger.info("Update complete")


def updated_data_tables():
    logger.info("Aggregating data tables.")
    add_small_data_to_store()
    munge_big_data()
    logger.info("Aggregation complete.")


if __name__ == "__main__":
    args = get_options()
    daskClient = Client(n_workers=args.n_workers, threads_per_worker=1)
    store = pd.HDFStore(DATA_STORE)
    patterns = get_patterns(Path(WORKDIR, "patterns.yaml").as_posix())

    if args.command == "init":
        if args.append:
            append_store()
        else:
            initialize_store()
    elif args.command == "queue":
        if args.queue_update:
            update_queue()
        elif args.queue_print:
            print_queue()
        elif args.queue_abi:
            dask_run_srr_checker(store["prealn/queue"], check_abi_solid, None, daskClient)
    elif args.command == "data":
        updated_data_tables()

    store.close()
    daskClient.close()
