from pathlib import Path
from typing import List, Dict
import sqlite3

from more_itertools import flatten
import joblib

from .db_connections import mongo, sqlite

PROJECT_DIR = Path(__file__).absolute().parents[3]

RNASEQ_SRXS_PATH = PROJECT_DIR / "output/library_strategy-wf/rnaseq_inliers.pkl"
ANATOMY_OBO = PROJECT_DIR / "data/fly_anatomy.obo"
DEVELOPMENT_OBO = PROJECT_DIR / "data/fly_development.obo"
CELL_TSV = PROJECT_DIR / "data/fly_cell_lines.tsv"


def get_rnaseq_srxs() -> List[str]:
    return joblib.load(RNASEQ_SRXS_PATH)


def _parse_obo(file_name):
    with open(file_name) as fh:
        for row in fh:
            if row.startswith("name:"):
                yield row.replace("name:", "").strip()


def get_fly_anatomy() -> List[str]:
    return [term for term in _parse_obo(ANATOMY_OBO) if not "cell-line" in term]


def get_fly_development() -> List[str]:
    return list(_parse_obo(DEVELOPMENT_OBO))


def get_fly_cell_line() -> List[str]:
    with open(CELL_TSV) as fh:
        next(fh)
        return [row.strip().split("\t")[1] for row in fh]


def sql_query_complete_biosamples() -> List[str]:
    """These samples are already in the database"""
    sql_query = "SELECT biosample FROM biometa WHERE complete = 1"
    with sqlite() as db:
        c = db.cursor()  # type: sqlite3.Cursor
        c.execute(sql_query)
        return list(flatten(c.fetchall()))


def sql_query_biosample(biosample) -> tuple:
    sql_query = "SELECT * FROM biometa WHERE biosample = ?"
    with sqlite() as db:
        cur = db.cursor()  # type: sqlite3.Cursor
        cur.execute(sql_query, (biosample,))
        return cur.fetchone()


def sql_update_biosample(results: List[tuple]):
    sql_upsert = """INSERT INTO biometa(biosample, sex, dev_stage, tissue, cell_type, perturbed, complete)
        VALUES (?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(biosample) DO UPDATE SET
            sex = excluded.sex,
            dev_stage = excluded.dev_stage,
            tissue = excluded.tissue,
            cell_type = excluded.cell_type,
            perturbed = excluded.perturbed,
            complete = excluded.complete
    """
    with sqlite() as db:  # type: sqlite3.Connection
        cur = db.cursor()  # type: sqlite3.Cursor
        cur.executemany(sql_upsert, results)
        db.commit()


def get_bioprojects(limit=100_000, skip=0) -> List[dict]:
    rnaseq_srxs = get_rnaseq_srxs()

    with mongo() as db:
        cursor = db.aggregate(
            [
                {"$match": {"srx": {"$in": rnaseq_srxs},}},
                {"$group": {"_id": "$BioProject.accn", "title": {"$first": "$BioProject.title"},}},
                {"$sort": {"_id": 1}},
                {"$skip": skip},
                {"$limit": limit},
            ]
        )

        return list(cursor)


def get_bioproject(bioproject: str) -> dict:
    rnaseq_srxs = get_rnaseq_srxs()

    with mongo() as db:
        cursor = db.aggregate(
            [
                {"$match": {"BioProject.accn": bioproject, "srx": {"$in": rnaseq_srxs},}},
                {
                    "$group": {
                        "_id": "$BioProject.accn",
                        "title": {"$first": "$BioProject.title"},
                        "description": {"$first": "$BioProject.description"},
                        "samples": {"$addToSet": "$sample"},
                    }
                },
            ]
        )
        bioproject = next(cursor)

    # Add on updated annotations from SQL database
    samples = []
    for sample in bioproject["samples"]:
        values = sql_query_biosample(sample["biosample"])
        if values:
            sample["sex"] = values[1]
            sample["dev_stage"] = values[2]
            sample["tissue"] = values[3]
            sample["cell_type"] = values[4]
            sample["perturbed"] = values[5]
            sample["complete"] = values[6]
        samples.append(sample)
    bioproject["samples"] = samples

    return bioproject

def query_term(col: str, term: str) -> dict:
    """Get a list of all samples with a given term"""
    rnaseq_srxs = get_rnaseq_srxs()

    # This allows me to query in different columns. I could not pass a column name without making
    # the app vulnerable to SQL injection.
    query_table = {
        "cell_type": "SELECT biosample FROM biometa WHERE cell_type LIKE ( '%' || ? || '%')",
        "tissue": "SELECT biosample FROM biometa WHERE tissue LIKE ( '%' || ? || '%')",
        "dev_stage": "SELECT biosample FROM biometa WHERE dev_stage LIKE ( '%' || ? || '%')",
        "sex": "SELECT biosample FROM biometa WHERE sex LIKE ( '%' || ? || '%')",
    }

    # Get all biosamples that match term
    with sqlite() as db:
        cur = db.cursor()  # type: sqlite3.Cursor
        cur.execute(query_table[col], (term, ))
        biosamples = list(flatten(cur.fetchall()))

    # Get their SRA attributes and build a fake bioproject
    with mongo() as db:
        cursor = db.aggregate(
            [
                {"$match": {"BioSample.accn": {"$in": biosamples}, "srx": {"$in": rnaseq_srxs},}},
                {
                    "$group": {
                        "_id": "$BioSample.accn",
                        "info": {"$first": "$sample"},
                    }
                },
            ]
        )

        # This helps me hijack the bioproject template system.
        mock_project = {
            "_id": f"{col} == {term}",
            "title": "",
            "description": "",
            "samples": [sample["info"] for sample in cursor]
        }

    # Add on updated annotations from SQL database
    samples = []
    for sample in mock_project["samples"]:
        values = sql_query_biosample(sample["biosample"])
        if values:
            sample["sex"] = values[1]
            sample["dev_stage"] = values[2]
            sample["tissue"] = values[3]
            sample["cell_type"] = values[4]
            sample["perturbed"] = values[5]
            sample["complete"] = values[6]
        samples.append(sample)
    mock_project["samples"] = samples

    return mock_project
