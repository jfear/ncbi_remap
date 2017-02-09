import mongoengine as me
from mongoengine import connect
import sys
sys.path.insert(0, '../lib/python')
from sramongo.mongo_schema import BioSample
from ncbi_remap import mongo_schema as ms
client = connect('sra', port=27022)

ms.Run.objects().count()

oliver = [x.sample_id for x in BioSample.objects(contacts__last_name='Oliver')]
rnaseq_exp = [x.pk for x in ms.Experiment.objects(db_flags='RNASeq')]

runs = ms.Run.objects(samples__in=oliver, experiment_id__in=rnaseq_exp, read_len_r1__gte=75, tax_analysis__tax_counts__exists=1, size_MB__lt=600)
runs.count()

tax = []
for r in runs:
    totals = r.tax_analysis.total_spots
    cnts = r.tax_analysis.tax_counts['melanogaster subgroup']['total_count']
    tax.append([r.pk, cnts/totals * 100])

from pprint import pprint
pprint(sorted(tax, key=lambda x: x[1]))

