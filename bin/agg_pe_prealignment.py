#!/usr/bin/env python
import os
import sys
import pandas as pd
from pymongo import MongoClient

from lcdblib.snakemake import helpers
from lcdblib.utils import utils

sys.path.insert(0, '../lib/python')
from ncbi_remap import parser
from ncbi_remap.parser import parse_files

patterns = {
    'atropos': {
        'r1': '../output/prealignment/raw/{experiment}/{sample}/{sample}_1.trim.clean.fastq.gz.log',
    },
    'fastq_screen': '../output/prealignment/raw/{experiment}/{sample}/{sample}_1.clean_screen.txt',
    'fastqc': {
        'zip': '../output/prealignment/raw/{experiment}/{sample}/{sample}_1_fastqc.zip',
    },
    'hisat2': '../output/prealignment/raw/{experiment}/{sample}/{sample}.fq.bam.log',
    'rseqc': {
        'bam_stat': '../output/prealignment/raw/{experiment}/{sample}/{sample}.bam_stat.txt',
        'infer_experiment': '../output/prealignment/raw/{experiment}/{sample}/{sample}.infer_experiment.txt',
    },
    'feature_counts': {
        'first': {
            'counts': '../output/prealignment/raw/{experiment}/{sample}/{sample}_1.fq.bam.counts',
            'summary': '../output/prealignment/raw/{experiment}/{sample}/{sample}_1.fq.bam.counts.summary',
        },
        'second': {
            'counts': '../output/prealignment/raw/{experiment}/{sample}/{sample}_2.fq.bam.counts',
            'summary': '../output/prealignment/raw/{experiment}/{sample}/{sample}_2.fq.bam.counts.summary',
        },
    },
    'picard': {
        'collectrnaseqmetrics': {
            'metrics': {
                'first': '../output/prealignment/raw/{experiment}/{sample}/{sample}_FIRST_READ_TRANSCRIPTION_STRAND.fq.bam.picard.collectrnaseqmetrics',
                'second': '../output/prealignment/raw/{experiment}/{sample}/{sample}_SECOND_READ_TRANSCRIPTION_STRAND.fq.bam.picard.collectrnaseqmetrics',
            },
        },
        'markduplicates': {
            'metrics': '../output/prealignment/raw/{experiment}/{sample}/{sample}.fq.bam.picard.markduplicatesmetrics',
        },
    },
    'samtools_stats': '../output/prealignment/raw/{experiment}/{sample}/{sample}.fq.bam.samtools.stats',
    'samtools_idxstats': '../output/prealignment/raw/{experiment}/{sample}/{sample}.fq.bam.samtools.idxstats',
    'bamtools_stats': '../output/prealignment/raw/{experiment}/{sample}/{sample}.fq.bam.bamtools.stats',
}

agg = {
    'atropos': '../output/prealignment/agg/pe_atropos.tsv',
    'fastq_screen': '../output/prealignment/agg/pe_fastq_screen.tsv',
    'hisat2': '../output/prealignment/agg/pe_hisat2.tsv',
    'rseqc': {
        'infer_experiment': '../output/prealignment/agg/pe_inferExperiment.tsv',
        'bam_stat': '../output/prealignment/agg/pe_bamStat.tsv',
    },
    'feature_counts': {
        'first': {
            'counts': '../output/prealignment/agg/pe_feature_counts_first.tsv',
            'summary': '../output/prealignment/agg/pe_feature_counts_summary_first.tsv',
            },
        'second': {
            'counts': '../output/prealignment/agg/pe_feature_counts_second.tsv',
            'summary': '../output/prealignment/agg/pe_feature_counts_summary_second.tsv',
            },
    },
    'picard': {
        'markduplicates': '../output/prealignment/agg/pe_picard_markduplicates.tsv',
        'collectrnaseqmetrics': {
            'first': {
                'summary': '../output/prealignment/agg/pe_picard_collectrnaseqmetrics_first_summary.tsv',
                'hist': '../output/prealignment/agg/pe_picard_collectrnaseqmetrics_first_hist.tsv',
            },
            'second': {
                'summary': '../output/prealignment/agg/pe_picard_collectrnaseqmetrics_first_summary.tsv',
                'hist': '../output/prealignment/agg/pe_picard_collectrnaseqmetrics_first_hist.tsv',
            },
        }
    },
    'samtools_stats': '../output/prealignment/agg/pe_samtools_stats.tsv',
    'samtools_idxstats': '../output/prealignment/agg/pe_samtools_idxstats.tsv',
    'bamtools_stats': '../output/prealignment/agg/pe_bamtools_stats.tsv',
    'fastqc': {
        'per_seq_quality': '../output/prealignment/agg/pe_fastqc_per_seq_quality.tsv',
        'per_base_seq_quality': '../output/prealignment/agg/pe_fastqc_per_base_seq_quality.tsv',
        'adapter_content': '../output/prealignment/agg/pe_fastqc_adapter_content.tsv',
        'per_base_seq_content': '../output/prealignment/agg/pe_fastqc_per_base_seq_content.tsv',
        'sequence_length': '../output/prealignment/agg/pe_fastqc_sequence_length.tsv',
        'overrepresented_seq': '../output/prealignment/agg/pe_fastqc_overrepresented_seq.tsv',
        'basic_stats': '../output/prealignment/agg/pe_fastqc_basic_stats.tsv',
        'kmer_content': '../output/prealignment/agg/pe_fastqc_kmer_content.tsv',
        'per_base_n_content': '../output/prealignment/agg/pe_fastqc_per_base_n_content.tsv',
        'per_seq_gc_content': '../output/prealignment/agg/pe_fastqc_per_seq_gc_content.tsv',
        'seq_dup_level': '../output/prealignment/agg/pe_fastqc_seq_dup_level.tsv',
    },
}

parser = {
    'fastq_md5': parser.parse_md5,
    'atropos': parser.parse_atropos,
    'fastq_screen': parser.parse_fqscreen,
    'hisat2': parser.parse_hisat2_pe,
    'rseqc': {
        'infer_experiment': parser.parse_inferExperiment,
        'bam_stat': parser.parse_bamStat,
    },
    'feature_counts': {
        'counts': parser.parse_featureCounts_counts,
        'summary': parser.parse_featureCounts_summary,
    },
    'picard': {
        'markduplicates': parser.parse_picard_markduplicate_metrics,
        'collectrnaseqmetrics': {
            'summary': parser.parse_picardCollect_summary,
            'hist': parser.parse_picardCollect_hist,
        }
    },
    'samtools_stats': parser.parse_samtools_stats,
    'samtools_idxstats': parser.parse_samtools_idxstats,
    'bamtools_stats': parser.parse_bamtools_stats,
    'fastqc': {
        'per_seq_quality': parser.parse_fastqc_per_seq_quality,
        'per_base_seq_quality': parser.parse_fastqc_per_base_seq_quality,
        'adapter_content': parser.parse_fastqc_adapter_content,
        'per_base_seq_content': parser.parse_fastqc_per_base_seq_content,
        'sequence_length': parser.parse_fastqc_sequence_length,
        'overrepresented_seq': parser.parse_fastqc_overrepresented_seq,
        'basic_stats': parser.parse_fastqc_basic_stats,
        'kmer_content': parser.parse_fastqc_kmer_content,
        'per_base_n_content': parser.parse_fastqc_per_base_n_content,
        'per_seq_gc_content': parser.parse_fastqc_per_seq_gc_content,
        'seq_dup_level': parser.parse_fastqc_seq_dup_level,
    },
}

def aggregate(files, pattern, fout, parser):
    if not os.path.exists(fout):
        df = parse_files(files, pattern, parser)
        df.to_csv(fout, sep='\t')

def main():
    # Connect to mongodb
    with open('.mongodb_host', 'r') as fh:
        mongo_client = MongoClient(host=fh.read().strip(), port=27022)
        db = mongo_client['sra2']
        remap = db['remap']

    # Build sample table
    samples = remap.aggregate([
        {'$unwind': '$runs'},
        {
            '$match': {
                '$and': [
                    {'runs.srr': {'$exists': 1}},
                    {'runs.pre_aln_flags': {'$eq': 'complete'}},
                    {'runs.pre_aln_flags': 'PE'},
                ]
            }
        },
        {'$project': {'_id': 0, 'experiment': '$_id', 'sample': '$runs.srr'}},
        {'$sort': {'sample': 1}},
    ])

    sample_table = pd.DataFrame(list(samples))
    targets = helpers.fill_patterns(patterns, sample_table)

    os.makedirs('../output/prealignment/agg', exist_ok=True)        # Make agg direcotry if does not exists

    aggregate(targets['atropos']['r1'], patterns['atropos']['r1'], agg['atropos'], parser['atropos'])
    aggregate(targets['fastq_screen'], patterns['fastq_screen'], agg['fastq_screen'], parser['fastq_screen'])
    aggregate(targets['hisat2'], patterns['hisat2'], agg['hisat2'], parser['hisat2'])
    aggregate(targets['rseqc']['infer_experiment'], patterns['rseqc']['infer_experiment'], agg['rseqc']['infer_experiment'], parser['rseqc']['infer_experiment'])
    aggregate(targets['rseqc']['bam_stat'], patterns['rseqc']['bam_stat'], agg['rseqc']['bam_stat'], parser['rseqc']['bam_stat'])
#     aggregate(targets['feature_counts']['first']['counts'], patterns['feature_counts']['first']['counts'], agg['feature_counts']['first']['counts'], parser['feature_counts']['counts'])
    aggregate(targets['feature_counts']['first']['summary'], patterns['feature_counts']['first']['summary'], agg['feature_counts']['first']['summary'], parser['feature_counts']['summary'])
#     aggregate(targets['feature_counts']['second']['counts'], patterns['feature_counts']['second']['counts'], agg['feature_counts']['second']['counts'], parser['feature_counts']['counts'])
    aggregate(targets['feature_counts']['second']['summary'], patterns['feature_counts']['second']['summary'], agg['feature_counts']['second']['summary'], parser['feature_counts']['summary'])
    aggregate(targets['picard']['collectrnaseqmetrics']['metrics']['first'], patterns['picard']['collectrnaseqmetrics']['metrics']['first'], agg['picard']['collectrnaseqmetrics']['first']['summary'], parser['picard']['collectrnaseqmetrics']['summary'])
    aggregate(targets['picard']['collectrnaseqmetrics']['metrics']['first'], patterns['picard']['collectrnaseqmetrics']['metrics']['first'], agg['picard']['collectrnaseqmetrics']['first']['hist'], parser['picard']['collectrnaseqmetrics']['hist'])
    aggregate(targets['picard']['collectrnaseqmetrics']['metrics']['second'], patterns['picard']['collectrnaseqmetrics']['metrics']['second'], agg['picard']['collectrnaseqmetrics']['second']['summary'], parser['picard']['collectrnaseqmetrics']['summary'])
    aggregate(targets['picard']['collectrnaseqmetrics']['metrics']['second'], patterns['picard']['collectrnaseqmetrics']['metrics']['second'], agg['picard']['collectrnaseqmetrics']['second']['hist'], parser['picard']['collectrnaseqmetrics']['hist'])
    aggregate(targets['picard']['markduplicates']['metrics'], patterns['picard']['markduplicates']['metrics'], agg['picard']['markduplicates'], parser['picard']['markduplicates'])
    aggregate(targets['samtools_stats'], patterns['samtools_stats'], agg['samtools_stats'], parser['samtools_stats'])
    aggregate(targets['samtools_idxstats'], patterns['samtools_idxstats'], agg['samtools_idxstats'], parser['samtools_idxstats'])
    aggregate(targets['bamtools_stats'], patterns['bamtools_stats'], agg['bamtools_stats'], parser['bamtools_stats'])
    aggregate(targets['fastqc']['zip'], patterns['fastqc']['zip'], agg['fastqc']['per_seq_quality'], parser['fastqc']['per_seq_quality'])
    aggregate(targets['fastqc']['zip'], patterns['fastqc']['zip'], agg['fastqc']['per_base_seq_quality'], parser['fastqc']['per_base_seq_quality'])
    aggregate(targets['fastqc']['zip'], patterns['fastqc']['zip'], agg['fastqc']['adapter_content'], parser['fastqc']['adapter_content'])
    aggregate(targets['fastqc']['zip'], patterns['fastqc']['zip'], agg['fastqc']['per_base_seq_content'], parser['fastqc']['per_base_seq_content'])
    aggregate(targets['fastqc']['zip'], patterns['fastqc']['zip'], agg['fastqc']['sequence_length'], parser['fastqc']['sequence_length'])
    aggregate(targets['fastqc']['zip'], patterns['fastqc']['zip'], agg['fastqc']['overrepresented_seq'], parser['fastqc']['overrepresented_seq'])
    aggregate(targets['fastqc']['zip'], patterns['fastqc']['zip'], agg['fastqc']['basic_stats'], parser['fastqc']['basic_stats'])
    aggregate(targets['fastqc']['zip'], patterns['fastqc']['zip'], agg['fastqc']['kmer_content'], parser['fastqc']['kmer_content'])
    aggregate(targets['fastqc']['zip'], patterns['fastqc']['zip'], agg['fastqc']['per_base_n_content'], parser['fastqc']['per_base_n_content'])
    aggregate(targets['fastqc']['zip'], patterns['fastqc']['zip'], agg['fastqc']['per_seq_gc_content'], parser['fastqc']['per_seq_gc_content'])
    aggregate(targets['fastqc']['zip'], patterns['fastqc']['zip'], agg['fastqc']['seq_dup_level'], parser['fastqc']['seq_dup_level'])



if __name__ == '__main__':
    main()
