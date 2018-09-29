#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Sample SNPs for computation of r2
# ----------------------------------------------------------------------------------------

module load vcftools

CHR=$1

# --- Figure out which SNPs meet our criteria for inclusion and
# --- make a list of them from which to sample

echo -e "3R\t1"

vcftools --vcf data/chr$CHR.pass.snp.phased.haps.vcf \
    --maf 0.1 \
    --recode \
    --exclude-bed data/inversion_heterochromatin.bed \
    --out results/chr$CHR.passing_SNPs_for_r2.full

# --- Create subsample of sites
shuf results/chr$CHR.passing_SNPs_for_r2.full.recode.vcf | \
    head -n 101 | \
    cut -f 1-2 > results/chr$CHR.passing_SNPs_for_r2.sample.txt
