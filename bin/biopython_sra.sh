#!/bin/bash

if [ ! -e ~/.biopython/Bio/Entrez/XSDs ]; then mkdir -p ~/.biopython/Bio/Entrez/XSDs; fi
cd ~/.biopython/Bio/Entrez/XSDs

wget -O SRA.common.xsd "http://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/sra/doc/SRA_1-5/SRA.common.xsd?view=co"
wget -O SRA.submission.xsd "http://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/sra/doc/SRA_1-5/SRA.submission.xsd?view=co"
wget -O SRA.study.xsd "http://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/sra/doc/SRA_1-5/SRA.study.xsd?view=co"
wget -O SRA.sample.xsd "http://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/sra/doc/SRA_1-5/SRA.sample.xsd?view=co"
wget -O SRA.experiment.xsd "http://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/sra/doc/SRA_1-5/SRA.experiment.xsd?view=co"
wget -O SRA.run.xsd "http://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/sra/doc/SRA_1-5/SRA.run.xsd?view=co"
wget -O SRA.analysis.xsd "http://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/sra/doc/SRA_1-5/SRA.analysis.xsd?view=co"
