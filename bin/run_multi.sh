#!/bin/bash
source activate ncbi_remap

multiqc -c ../config/multiqc_config.yaml -o ../output/prealignment/ ../output/prealignment/raw
