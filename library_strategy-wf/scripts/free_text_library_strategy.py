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

import nltk
from nltk.tokenize import regexp_tokenize
from nltk.corpus import stopwords
from nltk.stem import WordNetLemmatizer

from ncbi_remap.nlp import lookups


def main():
    # Build a list of documents
    documents = [
        get_documents(x)
        for x in ncbi.aggregate(
            [
                {
                    "$project": {
                        "_id": False,
                        "srx": "$srx",
                        "exp_title": "$title",
                        "exp_design": "$design",
                        "exp_library_name": "$library_name",
                        "exp_construction": "$library_construction_protocol",
                    }
                }
            ]
        )
    ]

    # tokenize documents
    download_corpus()
    eng_stops = stopwords.words("english") + lookups["stopwords"]
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
    known_strategies = set(lookups["library_strategy"].values())
    doc_strategies = []
    for doc in tokenized_documents:
        token_strategies = []
        for token in doc:
            if token in known_strategies:
                token_strategies.append(token)

        # count the number of times a strategy appears and keep only the top ones.
        arr, cnt = np.unique(np.array(token_strategies), return_counts=True)
        cnts = list(zip(arr, cnt))
        # Pull out the best(s) library strategy based on number occurrences
        if len(cnts) > 0:
            keeps = []
            high = 0
            for v, c in sorted(cnts, key=lambda x: x[1], reverse=True):
                if c >= high:
                    keeps.append(v)
                    high = c

            string = "|".join(keeps)
        else:
            string = np.nan

        doc_strategies.append(string)

    # Output table
    df = pd.DataFrame(doc_strategies, index=np.asarray(documents)[:, 0])
    df.index.name = "SRX"
    df.columns = ["freetxt_library_strategy"]
    df.to_parquet(snakemake.output[0])


def download_corpus():
    try:
        nltk.data.find("corpora/stopwords")
    except LookupError:
        nltk.download("stopwords")

    try:
        nltk.data.find("corpora/wordnet")
    except LookupError:
        nltk.download("wordnet")


def get_documents(doc):
    srx = doc["srx"]
    del doc["srx"]

    for k, v in doc.items():
        if v is None:
            doc[k] = ""

    txt = (
        " ".join(doc.values())
        .lower()
        .replace("_", " ")
        .replace("hi-c", "hic")
        .replace("3-c", "3c")
        .replace("4-c", "4c")
        .replace("3'", "3prime")
        .replace("5'", "5prime")
        .replace("-", " ")
        .replace("sequencing", "seq")
        .replace("sequenced", "seq")
        .replace("sequence", "seq")
    )

    # Translate based on known phrases in lookup table
    for k, v in lookups["library_strategy"].items():
        txt = txt.replace(k, v)

    return srx, txt


if __name__ == "__main__":
    # Connect to the database
    mongoClient = MongoClient(host="localhost", port=27017)
    db = mongoClient["sramongo"]
    ncbi = db["ncbi"]

    try:
        main()
    finally:
        mongoClient.close()
