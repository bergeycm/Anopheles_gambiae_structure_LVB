#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Convert Ssese + Ag1000G haps files to VCF format
# ----------------------------------------------------------------------------------------

chr=$1

SAMP_FILE=data/ssese_with_ag1000g/chr${chr}.pass.snp.phased.ag1000g.sample
cat data/ag1000g.phase1.ar3/ag1000g.phase1.ar3.1.haplotypes.${chr}.sample > $SAMP_FILE
grep "^ssese" data/chr${chr}.pass.snp.phased.sample | cut -d' ' -f 1-3   >> $SAMP_FILE

shapeit -convert \
    --input-haps data/ssese_with_ag1000g/chr${chr}.pass.snp.phased.ag1000g \
    --output-vcf data/ssese_with_ag1000g/chr${chr}.pass.snp.phased.ag1000g.vcf
