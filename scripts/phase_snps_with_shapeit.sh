#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Phase SNPs with SHAPEIT
# ----------------------------------------------------------------------------------------

CHR=$1

module load shapeit

shapeit --input-vcf data/chr${CHR}.pass.snp.flt.vcf.gz \
        -O data/chr${CHR}.pass.snp.phased \
        --effective-size 1000000 \
        --thread 12 \
        --noped

exit
