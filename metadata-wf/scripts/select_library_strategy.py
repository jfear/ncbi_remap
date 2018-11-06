"""Take all the different library startegies and select the best one."""
import pandas as pd


def annot_library_strategy(x):
    """Selection function"""
    # If all same
    if (x.sra == x.author) & (x.sra == x['data']):
        return x.sra

    # If there is author metadata
    if x.author:
        # Re-annotate OTHER samples if author and machine match
        if (x.sra == 'OTHER') & (x.author == x['data']):
            return x.author

        # Re-annotate OTHER samples if author has a single value
        if (x.sra == 'OTHER') & ~('|' in x.author):
            return x.author

        # Ambiguous so output evertyhing, but keep author upfront.
        if (x.sra != x['data']):
            string = x.author

            if x['data'] not in string:
                string = string + '|' + x['data']

            if x.sra not in string:
                string = string + '|' + x.sra

            return string

    # If there is no author metadata
    # If SRA/Machine same
    if (x.sra == x['data']):
        return x.sra

    # EST and RNA-Seq are really similar, so if the SRA uses EST then stick
    # with that.
    if (x.sra == 'EST') & (x['data'] == 'RNA-Seq'):
        return 'EST'

    # WGA and WGS are really similar, so if the SRA uses WGA then stick
    # with that.
    if (x.sra == 'WGA') & (x['data'] == 'WGS'):
        return 'WGA'

    # Ambiguous so output evertyhing, but keep machine up front
    return '|'.join([x['data'], x.sra])


# Grab the different versions of library strategy annotation
sra = pd.read_parquet(snakemake.input['sra'])
sra.columns = ['sra']

author = pd.read_parquet(snakemake.input['free_text'])
author.columns = ['author']

data = pd.read_parquet(snakemake.input['forest'])
data = data.reset_index().melt(id_vars='srx').groupby('srx').value.value_counts().unstack().idxmax(axis=1).to_frame()
data.columns = ['data']

other = pd.read_parquet(snakemake.input['other'])
other = other.reset_index().melt(id_vars='srx').groupby('srx').value.value_counts().unstack().idxmax(axis=1).to_frame()
other.columns = ['data']
data = pd.concat([data, other])

libstrat = pd.concat([sra, author, data], sort=True, axis=1)

# remove samples with no data, these have not completed the workflow
libstrat = libstrat[~libstrat.data.isnull()]

# Re-annotate samples
df = libstrat.apply(annot_library_strategy, axis=1)
df.name = 'Fear_et_al_library_strategy'
df.index.name = 'srx'
df.to_frame().to_parquet(snakemake.output[0])
