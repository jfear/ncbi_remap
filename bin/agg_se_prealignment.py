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
    'fastq_clean': {
        'r1': '../output/prealignment/raw/{experiment}/{sample}/{sample}_1.clean.fastq.gz',
    },
    'fastq_md5': '../output/prealignment/raw/{experiment}/{sample}/{sample}.md5',
    'fastq_clean_count': {
        'r1': '../output/prealignment/raw/{experiment}/{sample}/{sample}_1.clean.fastq.gz.libsize',
    },
    'atropos': {
        'r1': '../output/prealignment/raw/{experiment}/{sample}/{sample}_1.trim.clean.fastq.gz',
    },
    'fastq_atropos_count': {
        'r1': '../output/prealignment/raw/{experiment}/{sample}/{sample}_1.trim.clean.fastq.gz.libsize',
    },
    'fastq_screen': '../output/prealignment/raw/{experiment}/{sample}/{sample}_1.clean_screen.txt',
    'fastqc': {
        'html': '../output/prealignment/raw/{experiment}/{sample}/{sample}_fastqc.html',
        'zip': '../output/prealignment/raw/{experiment}/{sample}/{sample}_fastqc.zip',
    },
    'bed12': '../output/dm6_r6-11.bed12',
    'hisat2': {
        'splice_sites': '../output/known_splice_sites_r6-11.txt',
        'bam': '../output/prealignment/raw/{experiment}/{sample}/{sample}.fq.bam',
    },
    'bai': '../output/prealignment/raw/{experiment}/{sample}/{sample}.fq.bam.bai',
    'rseqc': {
        'bam_stat': '../output/prealignment/raw/{experiment}/{sample}/{sample}.bam_stat.txt',
        'infer_experiment': '../output/prealignment/raw/{experiment}/{sample}/{sample}.infer_experiment.txt',
        'geneBodyCoverage': {
            'txt': '../output/prealignment/raw/{experiment}/{sample}/{sample}.geneBodyCoverage.txt',
            'r': '../output/prealignment/raw/{experiment}/{sample}/{sample}.geneBodyCoverage.r',
            'img': '../output/prealignment/raw/{experiment}/{sample}/{sample}.geneBodyCoverage.pdf',
        },
        'tin': {
            'table': '../output/prealignment/raw/{experiment}/{sample}/{sample}.tin.tsv',
            'summary': '../output/prealignment/raw/{experiment}/{sample}/{sample}.tin.txt',
        },
    },
    'dupradar': {
        'density_scatter': '../output/prealignment/raw/{experiment}/{sample}/{sample}.density_scatter.png',
        'expression_histogram': '../output/prealignment/raw/{experiment}/{sample}/{sample}.expression_histogram.png',
        'expression_boxplot': '../output/prealignment/raw/{experiment}/{sample}/{sample}.expression_boxplot.png',
        'expression_barplot': '../output/prealignment/raw/{experiment}/{sample}/{sample}.expression_barplot.png',
        'multimapping_histogram': '../output/prealignment/raw/{experiment}/{sample}/{sample}.multimapping_histogram.png',
        'dataframe': '../output/prealignment/raw/{experiment}/{sample}/{sample}.dupradar.tsv',
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
            'bam': '../output/prealignment/raw/{experiment}/{sample}/{sample}.picard.dups.fq.bam',
            'metrics': '../output/prealignment/raw/{experiment}/{sample}/{sample}.fq.bam.picard.markduplicatesmetrics',
        },
    },
    'samtools_stats': '../output/prealignment/raw/{experiment}/{sample}/{sample}.fq.bam.samtools.stats',
    'samtools_idxstats': '../output/prealignment/raw/{experiment}/{sample}/{sample}.fq.bam.samtools.idxstats',
    'bamtools_stats': '../output/prealignment/raw/{experiment}/{sample}/{sample}.fq.bam.bamtools.stats',
}

agg = {
    'fastq_md5': '../output/prealignment/agg/se_md5.tsv',
    'fastq_clean_count': '../output/prealignment/agg/se_fastq_clean_libsize.tsv',
    'atropos': '../output/prealignment/agg/se_atropos.tsv',
    'fastq_atropos_count': '../output/prealignment/agg/se_atropos_libsize.tsv',
    'fastq_screen': '../output/prealignment/agg/se_fastq_screen.tsv',
    'hisat2': '../output/prealignment/agg/se_hisat2.tsv',
    'rseqc': {
        'geneBodyCoverage': '../output/prealignment/agg/se_geneBodyCoverage.tsv',
        'infer_experiment': '../output/prealignment/agg/se_inferExperiment.tsv',
        'bam_stat': '../output/prealignment/agg/se_bamStat.tsv',
        'tin': '../output/prealignment/agg/se_tin.tsv',
    },
    'dupradar': '../output/prealignment/agg/se_dupradar.tsv',
    'feature_counts': {
        'first': {
            'counts': '../output/prealignment/agg/se_feature_counts_first.tsv',
            'summary': '../output/prealignment/agg/se_feature_counts_summary_first.tsv',
            },
        'second': {
            'counts': '../output/prealignment/agg/se_feature_counts_second.tsv',
            'summary': '../output/prealignment/agg/se_feature_counts_summary_second.tsv',
            },
    },
    'picard': {
        'markduplicates': '../output/prealignment/agg/se_picard_markduplicates.tsv',
        'collectrnaseqmetrics': {
            'first': {
                'summary': '../output/prealignment/agg/se_picard_collectrnaseqmetrics_first_summary.tsv',
                'hist': '../output/prealignment/agg/se_picard_collectrnaseqmetrics_first_hist.tsv',
            },
            'second': {
                'summary': '../output/prealignment/agg/se_picard_collectrnaseqmetrics_first_summary.tsv',
                'hist': '../output/prealignment/agg/se_picard_collectrnaseqmetrics_first_hist.tsv',
            },
        }
    },
    'samtools_stats': '../output/prealignment/agg/se_samtools_stats.tsv',
    'samtools_idxstats': '../output/prealignment/agg/se_samtools_idxstats.tsv',
    'bamtools_stats': '../output/prealignment/agg/se_bamtools_stats.tsv',
    'fastqc': {
        'per_seq_quality': '../output/prealignment/agg/se_fastqc_per_seq_quality.tsv',
        'per_base_seq_quality': '../output/prealignment/agg/se_fastqc_per_base_seq_quality.tsv',
        'adapter_content': '../output/prealignment/agg/se_fastqc_adapter_content.tsv',
        'per_base_seq_content': '../output/prealignment/agg/se_fastqc_per_base_seq_content.tsv',
        'sequence_length': '../output/prealignment/agg/se_fastqc_sequence_length.tsv',
        'overrepresented_seq': '../output/prealignment/agg/se_fastqc_overrepresented_seq.tsv',
        'basic_stats': '../output/prealignment/agg/se_fastqc_basic_stats.tsv',
        'kmer_content': '../output/prealignment/agg/se_fastqc_kmer_content.tsv',
        'per_base_n_content': '../output/prealignment/agg/se_fastqc_per_base_n_content.tsv',
        'per_seq_gc_content': '../output/prealignment/agg/se_fastqc_per_seq_gc_content.tsv',
        'seq_dup_level': '../output/prealignment/agg/se_fastqc_seq_dup_level.tsv',
    },
}

parser = {
    'fastq_md5': parser.parse_md5,
    'fastq_clean_count': parser.parse_libsize,
    'atropos': parser.parse_atropos,
    'fastq_atropos_count': parser.parse_libsize,
    'fastq_screen': parser.parse_fqscreen,
    'hisat2': parser.parse_hisat2,
    'rseqc': {
        'geneBodyCoverage': parser.parse_geneBodyCoverage,
        'infer_experiment': parser.parse_inferExperiment,
        'bam_stat': parser.parse_bamStat,
        'tin': parser.parse_tin,
    },
    'dupradar': parser.parse_dupradar,
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
                    {'runs.pre_aln_flags': 'SE'},
                ]
            }
        },
        {'$project': {'_id': 0, 'experiment': '$_id', 'sample': '$runs.srr'}},
        {'$sort': {'sample': 1}},
    ])

    sample_table = pd.DataFrame(list(samples))
    targets = helpers.fill_patterns(patterns, sample_table)

    os.makedirs('../output/prealignment/agg', exist_ok=True)        # Make agg direcotry if does not exists

#     aggregate(targets['fastq_md5'], patterns['fastq_md5'], agg['fastq_md5'], parser['fastq_md5'])
    aggregate(targets['fastq_clean_count']['r1'], patterns['fastq_clean_count']['r1'], agg['fastq_clean_count'], parser['fastq_clean_count'])
    aggregate([x + '.log' for x in targets['atropos']['r1']], patterns['atropos']['r1'] + '.log', agg['atropos'], parser['atropos'])
    aggregate(targets['fastq_atropos_count']['r1'], patterns['fastq_atropos_count']['r1'], agg['fastq_atropos_count'], parser['fastq_atropos_count'])
    aggregate(targets['fastq_screen'], patterns['fastq_screen'], agg['fastq_screen'], parser['fastq_screen'])
#     aggregate([x + '.log' for x in targets['hisat2']['bam']], patterns['hisat2']['bam'] + '.log', agg['hisat2'], parser['hisat2'])
#     aggregate(targets['rseqc']['geneBodyCoverage']['txt'], patterns['rseqc']['geneBodyCoverage']['txt'], agg['rseqc']['geneBodyCoverage'], parser['rseqc']['geneBodyCoverage'])
    aggregate(targets['rseqc']['infer_experiment'], patterns['rseqc']['infer_experiment'], agg['rseqc']['infer_experiment'], parser['rseqc']['infer_experiment'])
    aggregate(targets['rseqc']['bam_stat'], patterns['rseqc']['bam_stat'], agg['rseqc']['bam_stat'], parser['rseqc']['bam_stat'])
#     aggregate(targets['rseqc']['tin']['summary'], patterns['rseqc']['tin']['summary'], agg['rseqc']['tin'], parser['rseqc']['tin'])
#     aggregate(targets['dupradar']['dataframe'], patterns['dupradar']['dataframe'], agg['dupradar'], parser['dupradar'])
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
#     aggregate(targets['samtools_idxstats'], patterns['samtools_idxstats'], agg['samtools_idxstats'], parser['samtools_idxstats'])
    aggregate(targets['bamtools_stats'], patterns['bamtools_stats'], agg['bamtools_stats'], parser['bamtools_stats'])
#     aggregate(targets['fastqc']['zip'], patterns['fastqc']['zip'], agg['fastqc']['per_seq_quality'], parser['fastqc']['per_seq_quality'])
#     aggregate(targets['fastqc']['zip'], patterns['fastqc']['zip'], agg['fastqc']['per_base_seq_quality'], parser['fastqc']['per_base_seq_quality'])
#     aggregate(targets['fastqc']['zip'], patterns['fastqc']['zip'], agg['fastqc']['adapter_content'], parser['fastqc']['adapter_content'])
#     aggregate(targets['fastqc']['zip'], patterns['fastqc']['zip'], agg['fastqc']['per_base_seq_content'], parser['fastqc']['per_base_seq_content'])
    aggregate(targets['fastqc']['zip'], patterns['fastqc']['zip'], agg['fastqc']['sequence_length'], parser['fastqc']['sequence_length'])
    aggregate(targets['fastqc']['zip'], patterns['fastqc']['zip'], agg['fastqc']['overrepresented_seq'], parser['fastqc']['overrepresented_seq'])
    aggregate(targets['fastqc']['zip'], patterns['fastqc']['zip'], agg['fastqc']['basic_stats'], parser['fastqc']['basic_stats'])
#     aggregate(targets['fastqc']['zip'], patterns['fastqc']['zip'], agg['fastqc']['kmer_content'], parser['fastqc']['kmer_content'])
#     aggregate(targets['fastqc']['zip'], patterns['fastqc']['zip'], agg['fastqc']['per_base_n_content'], parser['fastqc']['per_base_n_content'])
#     aggregate(targets['fastqc']['zip'], patterns['fastqc']['zip'], agg['fastqc']['per_seq_gc_content'], parser['fastqc']['per_seq_gc_content'])
#     aggregate(targets['fastqc']['zip'], patterns['fastqc']['zip'], agg['fastqc']['seq_dup_level'], parser['fastqc']['seq_dup_level'])



if __name__ == '__main__':
    main()
