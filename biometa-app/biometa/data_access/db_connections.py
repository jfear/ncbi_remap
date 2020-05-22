from pathlib import Path
from contextlib import contextmanager
import sqlite3

from pymongo import MongoClient

PROJECT_DIR = Path(__file__).absolute().parents[3]
SQL_PATH = PROJECT_DIR / "data/biometa.db"


@contextmanager
def mongo(*args, **kwargs):
    with MongoClient() as client:
        yield client["sramongo"]["ncbi"]


@contextmanager
def sqlite(*args, **kwargs):
    conn = sqlite3.connect(SQL_PATH)
    # If biometa table not there then initialize
    cur = conn.cursor()
    cur.execute("SELECT count(*) FROM sqlite_master WHERE type='table' AND name='biometa'")
    if cur.fetchone()[0] == 0:
        cur.execute(
            "CREATE TABLE biometa (biosample text PRIMARY KEY, sex text, devel_stage text, tissue text, cell_line text, notes text, genetic integer, diet integer, chemical integer, radiation integer, temperature integer, other integer, control integer)"
        )
        conn.commit()
    try:
        yield conn
    finally:
        conn.close()
