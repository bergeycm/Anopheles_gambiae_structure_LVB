#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Make GENEPOP files for each population and run NeEstimator
# ----------------------------------------------------------------------------------------

wget https://raw.githubusercontent.com/z0on/2bRAD_denovo/master/vcf2genepop.pl

module load vcftools

IN_VCF=data/chr3L.pass.snp.flt.vcf.gz

OUT_DIR=data/genepop_for_NeEstimator
mkdir -p $OUT_DIR

# --- Make VCF with just subset of sites so we can use the same SNPs in all pops

vcftools --gzvcf $IN_VCF \
    --exclude-bed data/inversion_heterochromatin.bed \
    --chr 3L \
    --thin 1000 \
    --recode \
    --out $OUT_DIR/SNP_set

# Colinear part of 3L: 1,815,119 - 4,264,713

# "After filtering, kept 7654 out of a possible 5338813 Sites"

# --- Grab individuals for each site and conver to GENEPOP format

for SITE in `ls data/ssese.seqids.is-*.txt | cut -d"." -f 3 | sed -e "s/is-//"`; do

    OUT_PREFIX=$OUT_DIR/chr3L.$SITE

    vcftools \
        --vcf $OUT_DIR/SNP_set.recode.vcf \
        --keep data/ssese.seqids.is-$SITE.txt \
        --mac 1 \
        --recode \
        --out $OUT_PREFIX

    perl vcf2genepop.pl vcf=$OUT_PREFIX.recode.vcf > $OUT_PREFIX.sanspop.gen

    head -n2 $OUT_PREFIX.sanspop.gen > $OUT_PREFIX.gen
    echo "POP $SITE" >> $OUT_PREFIX.gen
    tail -n +3 $OUT_PREFIX.sanspop.gen >> $OUT_PREFIX.gen

    # --- Call NeEstimator
    sed -e "s/SITE/$SITE/" data/NeEstimator_input_file.txt > ${SITE}_info

    ~/work/bin/NeEstimator/Ne2-1L i:${SITE}_info > $OUT_PREFIX.genepop.out

done

rm vcf2genepop.pl

# ----------------------------------------------------------------------------------------

grep "Estimated Ne" data/genepop_for_NeEstimator/chr3L.*.txt | \
    sed -e "s/.*chr3L\.//" -e "s/\.txt.*=//" -e "s/ \+/\t/g" | \
    cut -f 1,2 | sed -e "s/LD//" > \
    results/genepop_Ne_estimates.txt

grep -A1 "Parametric" data/genepop_for_NeEstimator/chr3L.*.txt | \
    sed -e "s/.*chr3L\.//" -e "s/Parametric//" -e "s/\.txt[^ ]\+ */ /" -e "s/ \+/\t/g" | \
    cut -f 1,2 | grep -v "\-\-" | sed -e "s/LD//" > \
    results/genepop_bounds.txt
