#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Count LVB SNPs after filtration (before LD filtering)
# ----------------------------------------------------------------------------------------

for VCF in `ls data/chr*.pass.snp.flt.vcf.gz`; do

    echo "Processing $VCF..."

    gunzip -c $VCF | wc -l

done > reports/snp_count_LVB_after_flt.txt

# 6066846 + 7370934 + 5338857 + 7354993 + 2438211 = 28,569,841 lines

# ----------------------------------------------------------------------------------------
# --- Count combined Ag+LVB SNPs
# ----------------------------------------------------------------------------------------

for VCF in `ls data/ssese_with_ag1000g/ssese_with_ag1000g.*.flt.strict.vcf`; do

    echo "Processing $VCF..."

    wc -l $VCF

done > reports/snp_count_AG+LVB_after_flt.txt

