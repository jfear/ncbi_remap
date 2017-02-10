#!/bin/bash

# Inputs
FNAME1=$1
FNAME2=$2

# Outputs
ONAME1=$3
ONAME2=$4

# Decompress gz files
DECOMP1=$(mktemp -p /lscratch/$SLURM_JOBID)
DECOMP2=$(mktemp -p /lscratch/$SLURM_JOBID)

gunzip -c $FNAME1 | paste - - - - | awk 'BEGIN{FS=" "; OFS="\t"}{print $1"\t"$0}' | sort -k1,1n > $DECOMP1
gunzip -c $FNAME2 | paste - - - - | awk 'BEGIN{FS=" "; OFS="\t"}{print $1"\t"$0}' | sort -k1,1n > $DECOMP2

# Merge data by first column
MERGED=$(mktemp -p /lscratch/$SLURM_JOBID)
join -t $'\t' -j 1 -o 1.2,1.3,1.4,1.5,2.2,2.3,2.4,2.5 $DECOMP1 $DECOMP2 > $MERGED

# Split back out
OUT1=$(mktemp -p /lscratch/$SLURM_JOBID)
OUT2=$(mktemp -p /lscratch/$SLURM_JOBID)

awk -v out1=$OUT1 -v out2=$OUT2 'BEGIN{FS="\t"}{for (i=1; i<=4; i++){print $i >> out1}; for(i=5;i<=8; i++){print $i >> out2}}' $MERGED

# Zip it back up
gzip -c $OUT1 > $ONAME1
gzip -c $OUT2 > $ONAME2

# Clean up
rm $DECOMP1 $DECOMP2
rm $MERGED
rm $OUT1 $OUT2
