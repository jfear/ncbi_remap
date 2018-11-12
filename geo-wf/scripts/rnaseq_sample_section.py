import pandas as pd

OUTPUT = '../output/geo-wf/rnaseq_sample_section.tsv'
HEADER = ['title', 'organism', 'study', 'runs', 'GEO Experiment', 'GEO Sample', 'BioSample ID', 'BioProject',
          'pubmed', 'pubmed_title', 'pubmed_citation', 'pubmed_authors', 'contact',
          'characteristics: sex', 'characteristics: developmental stage', 'characteristics: tissue',
          'characteristics: cell type', 'molecule', 'description', 'raw file', ]


def main():
    df = pd.read_csv('../output/geo-wf/rnaseq_metadata.tsv', sep='\t', index_col=0)

    char_map = {
        'sex': 'characteristics: sex',
        'tissue': 'characteristics: tissue',
        'developmental stage': 'characteristics: developmental stage',
        'cell type': 'characteristics: cell type',
    }

    df.rename(char_map, axis=1, inplace=True)
    df['molecule'] = 'total RNA'
    df[HEADER].to_csv(OUTPUT, sep='\t')


if __name__ == '__main__':
    main()
