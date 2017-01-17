"""Mongo document schema."""
from mongoengine import ListField, DictField, StringField, IntField, FloatField, \
    EmbeddedDocument, EmbeddedDocumentField

from sramongo import mongo_schema


# Run
class FQPreprocess(EmbeddedDocument):
    reads = IntField()


class Run(mongo_schema.Run):
    pipeline_flags = ListField(StringField(), default=list)

    fq_preprocess = EmbeddedDocumentField(FQPreprocess)


# Experiment
class Experiment(mongo_schema.Experiment):
    pass


# Sample
class Sample(mongo_schema.Sample):
    pass
