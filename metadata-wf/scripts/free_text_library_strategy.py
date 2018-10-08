"""Parse SRA study and experiment free text for library strategy information.

The SRA limits the terms used for library strategy. Sometimes authors add more
details about what strategy they used. Keep in mind that I am parsing study
level information, so there may be a mix of strategies used.

This script harnesses a lookup table that I created by hand. Essentially, after
tokenization, I did bi-gram and tri-gram searches to figure out what was used and
create the lookup table. There maybe technologies missed, but I should capture
most of the information.

"""
import numpy as np
import pandas as pd

from pymongo import MongoClient

from nltk.tokenize import regexp_tokenize
from nltk.corpus import stopwords
from nltk.stem import WordNetLemmatizer

from ncbi_remap.nlp import lookups


def get_documents(doc):
    srx = doc['_id']
    del doc['_id']
    txt = ' '.join(doc.values()).lower()\
        .replace('_', ' ')\
        .replace('hi-c', 'hic')\
        .replace('3-c', '3c')\
        .replace('4-c', '4c')\
        .replace("3'", '3prime')\
        .replace('-', ' ')\
        .replace('sequencing', 'seq')\
        .replace('sequenced', 'seq')\
        .replace('sequence', 'seq')

    # Translate based on known phrases in lookup table
    for k, v in lookups['library_strategy'].items():
        txt = txt.replace(k, v)

    return srx, txt


def main(ncbi):
    # Build a list of documents
    documents = [get_documents(x) for x in ncbi.aggregate([
        {
            '$project': {
                'title': '$sra.study.title',
                'abstract': '$sra.study.abstract',
                'type': '$sra.study.study_type',
                'exp_title': '$sra.experiment.title',
                'exp_design': '$sra.experiment.design',
                'exp_library_name': '$sra.experiment.library_name',
            }
        },
    ])]

    # tokenize documents
    eng_stops = stopwords.words('english') + lookups['stopwords']
    wordnets_lemmatizer = WordNetLemmatizer()

    tokenized_documents = []
    for doc in np.asarray(documents)[:, 1]:
        tokens = regexp_tokenize(doc, r"[\w-]+")

        lemma = []
        for token in tokens:
            if token.isnumeric():
                continue
            elif len(token) <= 1:
                continue
            elif token in eng_stops:
                continue

            lemma.append(wordnets_lemmatizer.lemmatize(token))

        # store
        tokenized_documents.append(lemma)

    # Pull out library strategy based on lookup table
    known_strategies = set(lookups['library_strategy'].values())
    doc_strategies = []
    for doc in tokenized_documents:
        token_strategies = set()
        for token in doc:
            if token in known_strategies:
                token_strategies.add(token)
        doc_strategies.append('|'.join(list(token_strategies)))

    # Output table
    df = pd.DataFrame(doc_strategies, index=np.asarray(documents)[:, 0])
    df.index.name = 'SRX'
    df.columns = ['freetxt_library_strategy']
    df.to_parquet(snakemake.output[0])


if __name__ == '__main__':
    # Connect to the database
    try:
        with open('../output/.mongodb_host', 'r') as fh:
            host = fh.read().strip()
    except FileNotFoundError:
        host = 'localhost'

    mongoClient = MongoClient(host=host, port=27017)
    db = mongoClient['sra']
    ncbi = db['ncbi']
    main(ncbi)
