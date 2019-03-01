"""Create a list of keywords for each sample."""
from collections import defaultdict
from itertools import chain
from functools import partial

import numpy as np
import pandas as pd
from nltk.tokenize import regexp_tokenize
from nltk.corpus import stopwords
from nltk.stem import WordNetLemmatizer

from gensim.corpora.dictionary import Dictionary
from gensim.models.tfidfmodel import TfidfModel

from ncbi_remap.logging import logger
from ncbi_remap.mongo import mongo_connect
from ncbi_remap.sample_lists import get_complete

OUTPUT = '../output/geo-wf/bow.parquet'

def parse_document(doc):
    srx = doc['srx']

    sample_title = doc.get('sample_title', '')
    attrs = doc.get('attrs', [])

    string = []
    string.append(sample_title)
    for attr in attrs:
        string.append(attr['value'])

    return srx, ' '.join([str(x) for x in string])


def create_documents(samples):
    client, db, ncbi = mongo_connect()

    results = ncbi.aggregate([
        {
            '$match': {
                '_id': {'$in': samples}
            }
        },
        {
            '$project': {
                '_id': 0,
                'srx': '$srx',
                'sample_title': '$sra.sample.title',
                'attrs': '$sra.sample.attributes',
            }
        }
    ])

    srxs = []
    docs = []
    for result in results:
        srx, doc = parse_document(result)
        srxs.append(srx)
        docs.append(doc)

    return srxs, docs


def tokenize_document(doc):
    # tokens document into individual words
    tokens = regexp_tokenize(doc.lower(), r"[\w\(\);\-\+\[\]\/]+")

    # remove punctuation
    alpha_num = [token for token in tokens if token.isalnum()]

    # remove lone numbers
    no_num = [token for token in alpha_num if not token.isnumeric()]

    # remove single characters
    not_single = [token for token in no_num if len(token) > 1]

    # remove stop words
    eng_stops = stopwords.words('english') + ['drosophila', 'melanogaster']
    no_stops = [token for token in not_single if token not in eng_stops]

    # lemmatize to remove plurls and word variations
    wordnet_lemmatizer = WordNetLemmatizer()
    lemma = [wordnet_lemmatizer.lemmatize(token) for token in no_stops]

    return lemma


def drop_unique(document, total_cnts):
    """Remove tokens that are only in one document."""
    res = []
    for token_id, token_cnt in document:
        if token_cnt == total_cnts[token_id]:
            continue

        res.append((token_id, token_cnt))
    return res


def human_weights(document, dictionary):
    sorted_weights = sorted(document, key=lambda w: w[1], reverse=True)

    human = []
    for wt in sorted_weights:
        if wt[1] > 0.25:
            word = dictionary.get(wt[0])
            human.append((word, wt[1]))

    return human


def main():
    samples = get_complete()

    # Create a set of cleaned documents
    srxs, docs = create_documents(samples)
    tokenized_documents = list(map(tokenize_document, docs))

    # Create a bag of words from each document
    dictionary = Dictionary(tokenized_documents)
    corpus = [dictionary.doc2bow(document) for document in tokenized_documents]

    # Get corpus level counts
    total_cnts = defaultdict(int)
    for word_id, word_cnt in chain.from_iterable(corpus):
        total_cnts[word_id] += word_cnt

    # Remove tokens that are only in a single document
    corpus_no_unique = list(map(partial(drop_unique, total_cnts=total_cnts), corpus))

    # Create term frequency document model
    tfidf = TfidfModel(corpus_no_unique)

    # Calculate weights
    tfidf_weights = [tfidf[document] for document in corpus_no_unique]

    # pull out the best terms and make them readable
    weights = list(map(partial(human_weights, dictionary=dictionary), tfidf_weights))

    # Concatenate words in weight order
    bows = []
    for srx, wts in zip(srxs, weights):
        if len(wts) == 0:
            string = ''
        else:
            string = '|'.join(np.asarray(wts)[:, 0])
        bows.append((srx, string))

    # Build output dataframe and save
    keywords = pd.DataFrame(bows, columns=['srx', 'keywords'])
    keywords.set_index('srx', inplace=True)
    keywords.to_parquet(OUTPUT)


if __name__ == '__main__':
    main()
    logger.info('Script complete')
