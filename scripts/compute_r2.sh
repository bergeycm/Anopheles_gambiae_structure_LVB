#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Compute r2
# ----------------------------------------------------------------------------------------

module load vcftools

CHR=$1
SITE=$2

SITE_LIST=data/ssese.seqids.is-$SITE.txt

vcftools --vcf results/chr$CHR.passing_SNPs_for_r2.full.recode.vcf \
    --ld-window-bp 50000 \
    --maf 0.1 \
    --keep $SITE_LIST \
    --hap-r2-positions results/chr$CHR.passing_SNPs_for_r2.sample.txt \
    --min-r2 0.0 \
    --exclude-bed data/inversion_heterochromatin.bed \
    --out results/chr$CHR.$SITE
