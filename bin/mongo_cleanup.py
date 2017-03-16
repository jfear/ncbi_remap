import os
from shutil import rmtree
from pymongo import MongoClient
client = MongoClient(port=27022)
db = client['sra2']
remap = db['remap']

pe = remap.aggregate([
    {'$unwind': '$runs'},
    {
        '$match': {
            'runs.srr': {'$exists': 1},
            'runs.pre_aln_flags': 'PE',
        }
    },
    {
        '$project': {
            '_id': 0,
            'srx': '$_id',
            'srr': '$runs.srr'
        }
    }
])

for doc in pe:
    path = '../output/prealignment/raw/{srx}/{srr}'.format(**doc)
    rmtree(path, ignore_errors=True)

print('done')

result = remap.update_many({'$and': [{'runs.pre_aln_flags': 'PE'},{'runs.pre_aln_flags': 'complete'}]}, {'$pull': {'runs.$.pre_aln_flags': 'complete'}})

list(remap.aggregate([
    {'$unwind': '$runs'},
    {
        '$match': {
            'runs.srr': {'$exists': 1},
            '$and': [
                {'runs.pre_aln_flags': 'PE'},
                {'runs.pre_aln_flags': 'complete'},
            ]
        }
    },
    {'$count': 'count'}
]))

