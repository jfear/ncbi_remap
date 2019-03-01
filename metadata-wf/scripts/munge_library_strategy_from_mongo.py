import pandas as pd
from pymongo import MongoClient


# Connect to MongoDB
try:
    with open('../output/.mongodb_host', 'r') as fh:
        host = fh.read().strip()
except FileNotFoundError:
    host = 'localhost'

mongoClient = MongoClient(host=host, port=27017)
db = mongoClient['sra']
ncbi = db['ncbi']

# Create table mapping SRX to library startegy
libstrat = pd.DataFrame(list(ncbi.aggregate([
    {
        '$project': {
            '_id': 0,
            'srx': '$srx',
            'library_strategy': '$sra.experiment.library_strategy'
        }
    }
]))).set_index('srx')

libstrat.to_parquet(snakemake.output[0])

# Generate table of summary counts
vcnts = libstrat.library_strategy.value_counts().to_frame()
vcnts.to_csv(snakemake.output[1], sep='\t')