"""Mongo document schema."""
from mongoengine import ListField, DictField, StringField, IntField, FloatField, \
    EmbeddedDocument, EmbeddedDocumentField, Document

from sramongo.mongo_schema import Pubmed
class Run(EmbeddedDocument):
    srr = StringField()
    pre_aln_flags = ListField(StringField(), default=list)
    libsize = DictField()
    avgReadLen = DictField()
    md5 = DictField()

class Sample(EmbeddedDocument):
    srs = StringField()
    biosample = StringField()
    gsm = StringField()

class Contacts(EmbeddedDocument):
    first_name = StringField()
    last_name = StringField()
    email = StringField()

class Remap(Document):
    """Drosophila remap data base schema.

    Attributes
    ----------
    srx: mongoengine.StringField
        The primary identifier (SRX/ERX/DRX) for an experiment.

    sra: mongoengine.StringField
        The submission identifier (SRA/ERA/DRA).

    srp: mongoengine.StringField
        The project/study identifier (SRP/ERP/DRP).

    bioproject: mongoengine.StringField
        The BioProject identifier (SRX/ERX/DRX).

    sample: mongoengine.EmbeddedDocument
        A document containing the following sample information\:

            srs: str
                Primary identifier for a sample (SRS/ERS/DRS).
            biosample: str
                Identifier from biosample.
            gsm: str
                Identifier for a GEO samples.

    runs: mongoengine.ListField
        A list containing the following run information\:

            srr: str
                Primary identifier for a run (SRR/ERR/DRR).

            pre_aln_flags: mongoengine.ListField
                A list strings containing various flags generated during the
                pre-alignment workflow.

            libsize: mongoengine.DictField
                A document containing 'R1' or ('R1 and 'R2) with the number of
                reads in the fastq.

            avgReadLen: mongoengine.DictField
                A document containing 'R1' or ('R1 and 'R2) with the average
                read length.

            md5: mongoengine.DictField
                A document containing 'R1' or ('R1 and 'R2) with the md5sum of
                the fastq.

    papers: mongoengine.ListField
        A list papers containing the following information\:

            pubmed_id: mongoengine.StringField
                The primary identifier for Pubmed. These are the accession number
                which begin with PMID.

            title: mongoengine.StringField
                Title of the paper.

            abstract: mongoengine.StringField
                Paper abstract.

            authors: mongoengine.ListField
                List of authors.

            citation: mongoengine.StringField
                Citation information for the paper.

            date_created: mongoengine.DateTimeField
                Date the pubmed entry was created.

            date_completed: mongoengine.DateTimeField
                Date the pubmed entry was completed.

            date_revised: mongoengine.DateTimeField
                Date the pubmed entry was last updated.

    contacts: mongoengine.ListField
        A list containing the following contact information\:

            fist_name: str
            last_name: str
            email: str

    pre_aln_flags: mongoengine.ListField
        A list strings containing various flags generated during the
        pre-alignment workflow.

    """
    srx = StringField(primary_key=True)
    sra = StringField()
    srp = StringField()
    bioproject = StringField()

    sample = EmbeddedDocumentField(Sample)

    runs = ListField(EmbeddedDocumentField(Run), default=list)
    papers = ListField(EmbeddedDocumentField(Pubmed), default=list)
    contacts = ListField(EmbeddedDocumentField(Contacts), default=list)
    pre_aln_flags = ListField(StringField(), default=list)
