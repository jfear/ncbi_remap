from pymongo import MongoClient


def mongo_connect(host='localhost'):
    """Connect to mongo SRA database.

    Parameters
    ----------
    host : str
        Hostname to connect to [default = localhost]

    Returns
    -------
    tuple : client, db, db['ncbi']

    """
    client = MongoClient(host=host, port=27017)
    db = client['sra']
    return client, db, db['ncbi']
