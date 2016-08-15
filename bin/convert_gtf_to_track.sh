# Load ucsc module
module load ucsc

# Create 2bit version of reference
#faToTwoBit /data/Oliverlab/references/genomes/Dmelanogaster/dm6/seq/dm6.fasta ../output/dm6.2bit

# Get chromsome lengths from twobit genome
#twoBitInfo ../output/dm6.2bit stdout | sort -k2rn > ../output/dm6.chrom.sizes

# copy over chr2L
awk '{print "chr"$0}' /data/Oliverlab/data/FlyBase/FB2015_04/dmel.chrom.sizes > ../output/dmel.chrom.sizes

# Convert gtf to genomepred format
# NOTE: FBtr0307759 and FBtr0084080 had exons on both strands
# NOTE: I was getting the following errors later in the bigBed step, it is easier to fix it here
#       End coordinate 23513713 bigger than chr2L size of 23513712 line 64791 of ../output/all.ucsc.bed
#       End coordinate 1348137 bigger than chr4 size of 1348131 line 320307 of ../output/all.ucsc.edit.bed
#       End coordinate 3667354 bigger than chrY size of 3667352 line 391216 of ../output/all.ucsc.edit.bed
# Because this is quick and dirty I am going to just edit it.
awk -F '\t' -v OFS='\t' '{
    if ($1 == "chr2L" && $5 > 23513712 ){
        $5="23513712";
    }; 
    if($1 == "chr4" && $5 > 1348131 ){
        $5="1348131"
    }; 
    if($1 == "chrY" && $5 > 3667352 ){
        $5="3667352"
    }; 

    print $0

    }' ../data/all.ucsc.gtf > ../output/all.ucsc.edit.gtf

gtfToGenePred -infoOut=../output/infoOut.txt -genePredExt ../output/all.ucsc.edit.gtf -allErrors ../output/all.ucsc.gp

# Check the genomepred file and make sure it is ok
genePredCheck ../output/all.ucsc.gp

# Convert genePred file to BED format
genePredToBed ../output/all.ucsc.gp stdout | sort -k1,1 -k2,2n > ../output/all.ucsc.bed

# Convert BED to bigBED
bedToBigBed -type=bed12 -extraIndex=name ../output/all.ucsc.bed ../output/dmel.chrom.sizes ../output/all.ucsc.bb

# Required for indexing step
#grep -v -e "^#" ../output/infoOut.txt | awk '{printf "%s\t%s,%s,%s,%s,%s\n", $1,$2,$3,$8,$9,$10}' > ../output/all.ucsc.nameIndex.txt 

# Create index for position search in browser
#ixIxx ../output/all.ucsc.nameIndex.txt ../output/all.ucsc.nameIndex.ix ../output/all.ucsc.nameIndex.ixx
