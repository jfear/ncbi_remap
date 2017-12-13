"""Mongo document schema."""
from mongoengine import ListField, DictField, MapField, StringField, IntField, FloatField, \
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


class FastqScreen(EmbeddedDocument):
    """Representation of FASTQ Screen results."""
    reads_processed_count = IntField()
    unmapped_count = IntField()
    unmapped_percent = FloatField()
    one_hit_one_library_count = IntField()
    one_hit_one_library_percent = FloatField()
    multiple_hits_one_library_count = IntField()
    multiple_hits_one_library_percent = FloatField()
    one_hit_multiple_libraries_count = IntField()
    one_hit_multiple_libraries_percent = FloatField()
    multiple_hits_multiple_libraries_count = IntField()
    multiple_hits_multiple_libraries_percent = FloatField()


class BamStat(EmbeddedDocument):
    non_unique = IntField()
    nonsplice_reads = IntField()
    optical_pcr_duplicates = IntField()
    proper_pair_map_to_different_chrom = IntField()
    qc_failed = IntField()
    read_1 = IntField()
    read_2 = IntField()
    reads_map_minus = IntField()
    reads_map_plus = IntField()
    reads_mapped_proper_pairs = IntField()
    splice_reads = IntField()
    total_records = IntField()
    unique = IntField()
    unmapped_reads = IntField()


class InferExperiment(EmbeddedDocument):
    same_strand = FloatField()
    opposite_strand = FloatField()
    undetermined = FloatField()


class Hisat2(EmbeddedDocument):
    """Container for Hisat2 alignemnt summary output.

    https://ccb.jhu.edu/software/hisat2/manual.shtml#alignment-summary
    """
    num_reads = IntField()
    num_reads_paired = IntField()
    num_reads_unpaired = IntField()
    num_reads_unaligned = IntField()
    num_reads_uniquely_aligned = IntField()
    num_multimappers = IntField()
    num_reads_aligned_discordantly = IntField()
    num_unaligned_mates = IntField()
    num_uniquely_aligned_mates = IntField()
    num_multimappers_mates = IntField()
    per_alignment = FloatField()


class SamtoolsStats(EmbeddedDocument):
    """The SN fields from samtools stats."""
    raw_total_sequences = IntField()
    filtered_sequences = IntField()
    sequences = IntField()
    is_sorted = IntField()
    first_fragments = IntField()
    last_fragments = IntField()
    reads_mapped = IntField()
    reads_mapped_and_paired = IntField()
    reads_unmapped = IntField()
    reads_properly_paired = IntField()
    reads_paired = IntField()
    reads_duplicated = IntField()
    reads_MQ0 = IntField()
    reads_QC_failed = IntField()
    non_primary_alignments = IntField()
    total_length = IntField()
    bases_mapped = IntField()
    bases_mapped_cigar = IntField()
    bases_trimmed = IntField()
    bases_duplicated = IntField()
    mismatches = IntField()
    error_rate = FloatField()
    average_length = FloatField()
    maximum_length = IntField()
    average_quality = FloatField()
    insert_size_average = FloatField()
    insert_size_standard_deviation = FloatField()
    inward_oriented_pairs = IntField()
    outward_oriented_pairs = IntField()
    pairs_with_other_orientation = IntField()
    pairs_on_different_chromosomes = IntField()


class SamtoolsIdxStats(EmbeddedDocument):
    chrom = StringField()
    length = IntField()
    num_mapped_reads = IntField()
    num_unmapped_reads = IntField()


class CollectRNAseqMetrics(EmbeddedDocument):
    PF_BASES = IntField()
    PF_ALIGNED_BASES = IntField()
    RIBOSOMAL_BASES = IntField()
    CODING_BASES = IntField()
    UTR_BASES = IntField()
    INTRONIC_BASES = IntField()
    INTERGENIC_BASES = IntField()
    IGNORED_READS = IntField()
    CORRECT_STRAND_READS = IntField()
    INCORRECT_STRAND_READS = IntField()
    PCT_RIBOSOMAL_BASES = FloatField()
    PCT_CODING_BASES = FloatField()
    PCT_UTR_BASES = FloatField()
    PCT_INTRONIC_BASES = FloatField()
    PCT_INTERGENIC_BASES = FloatField()
    PCT_MRNA_BASES = FloatField()
    PCT_USABLE_BASES = FloatField()
    PCT_CORRECT_STRAND_READS = FloatField()
    MEDIAN_CV_COVERAGE = FloatField()
    MEDIAN_5PRIME_BIAS = FloatField()
    MEDIAN_3PRIME_BIAS = FloatField()
    MEDIAN_5PRIME_TO_3PRIME_BIAS = FloatField()
    SAMPLE = StringField()
    LIBRARY = StringField()
    READ_GROUP = StringField()


class MarkDuplicates(EmbeddedDocument):
    LIBRARY = StringField()
    UNPAIRED_READS_EXAMINED = IntField()
    READ_PAIRS_EXAMINED = IntField()
    SECONDARY_OR_SUPPLEMENTARY_RDS = IntField()
    UNMAPPED_READS = IntField()
    UNPAIRED_READ_DUPLICATES = IntField()
    READ_PAIR_DUPLICATES = IntField()
    READ_PAIR_OPTICAL_DUPLICATES = IntField()
    PERCENT_DUPLICATION = FloatField()
    ESTIMATED_LIBRARY_SIZE = IntField()


class FeatureCounts(EmbeddedDocument):
    Assigned = IntField()
    Assigned_Junction = IntField()
    Unassigned_Ambiguity = IntField()
    Unassigned_MultiMapping = IntField()
    Unassigned_NoFeatures = IntField()
    Unassigned_Unmapped = IntField()
    Unassigned_MappingQuality = IntField()
    Unassigned_FragmentLength = IntField()
    Unassigned_Chimera = IntField()
    Unassigned_Secondary = IntField()
    Unassigned_Nonjunction = IntField()
    Unassigned_Duplicate = IntField()


class PreAlignmentWorkflow(EmbeddedDocument):
    hisat2 = EmbeddedDocumentField(Hisat2)
    samtools_stats = EmbeddedDocumentField(SamtoolsStats)
    samtools_idxstats = ListField(SamtoolsIdxStats)
    fastq_screen = MapField(EmbeddedDocumentField(FastqScreen))
    infer_experiment = EmbeddedDocumentField(InferExperiment)
    bam_stat = EmbeddedDocumentField(BamStat)
    picard_collectrnaseqmetrics = MapField(EmbeddedDocumentField(CollectRNAseqMetrics))
    picard_markduplicates = EmbeddedDocumentField(MarkDuplicates)
    featurecounts = EmbeddedDocumentField(FeatureCounts)


class AlignmentWorkflow(EmbeddedDocument):
    hisat2 = EmbeddedDocumentField(Hisat2)
    samtools_stats = EmbeddedDocumentField(SamtoolsStats)


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
    pre_aln_workflow = EmbeddedDocumentField(PreAlignmentWorkflow)
    aln_workflow = EmbeddedDocumentField(AlignmentWorkflow)
