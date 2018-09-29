#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Compute depth of coverage for each individual
# ----------------------------------------------------------------------------------------

module load vcftools/0.1.15
module load tabix

vcf-concat `ls data/chr*.pass.snp.flt.vcf.gz | grep -v "X"` | \
    bgzip -c > data/all.pass.snp.flt.vcf.gz

# --- Compute mean depth per individual

vcftools --gzvcf data/all.pass.snp.flt.vcf.gz --depth \
    --exclude-bed data/inversion_heterochromatin.bed \
    --out reports/all.pass.snp.flt

# --- Compute proportion and count of missing genotypes

vcftools --gzvcf data/all.pass.snp.flt.vcf.gz --missing-indv \
    --exclude-bed data/inversion_heterochromatin.bed \
    --out reports/all.pass.snp.flt

# --- Compute heterozygozity

vcftools --gzvcf data/all.pass.snp.flt.vcf.gz --het \
    --exclude-bed data/inversion_heterochromatin.bed \
    --out reports/all.pass.snp.flt
