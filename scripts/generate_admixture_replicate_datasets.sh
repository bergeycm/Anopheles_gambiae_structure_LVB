#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Make datasets for ADMIXTURE replicates
# ----------------------------------------------------------------------------------------

module load plink

SUFFIX=
if [[ "$1" == "missingtoref" ]]; then
    SUFFIX=.missingtoref
fi

SUBSET_DIR=data/ssese_with_ag1000g$SUFFIX/adm_subsets
mkdir -p $SUBSET_DIR
rm -fr $SUBSET_DIR/*

# Make file with range of regions we desire to include
TMP_RANGE_PRE=`mktemp -u tmp.included_ranges.XXXXX`
echo "3 1 37000000 chr3_included"        > ${TMP_RANGE_PRE}.3L
echo "4 15000000 41000000 chr3_included" > ${TMP_RANGE_PRE}.3R

for M_STR in "" "withM."; do

    if [[ "$M_STR" == "" ]]; then
        EXCLUDE="data/ag1000g.phase1.ar3/excluded.individuals.txt"
    else
        EXCLUDE="data/ag1000g.phase1.ar3/excluded.individuals.notM.txt"
    fi

    for CHR in 3L 3R; do

        # Make set with just region we want to include
        plink \
            --bfile data/ssese_with_ag1000g$SUFFIX/chr$CHR$SUFFIX.pass.snp.phased.ag1000g.strict \
            --make-set ${TMP_RANGE_PRE}.$CHR \
            --write-set --out ${TMP_RANGE_PRE}.$CHR.included

        # Filter for biallelic with MAF > 0.01, including only SNPs in non-het regions
        # And removing individuals we don't want (maybe including the M individuals)
        plink \
            --bfile data/ssese_with_ag1000g$SUFFIX/chr$CHR$SUFFIX.pass.snp.phased.ag1000g.strict \
            --extract ${TMP_RANGE_PRE}.$CHR.included.set \
            --remove-fam $EXCLUDE \
            --biallelic-only strict \
            --maf 0.01 \
            --make-bed \
            --out $SUBSET_DIR/chr$CHR.${M_STR}full
    done

    # Do replicates for 10 seeds...

    RANDOM=`date +%N | sed s/...$//`

    for i in `seq 1 10`; do

        SEED=$RANDOM

        for CHR in 3L 3R; do

            plink \
                --seed $SEED \
                --bfile $SUBSET_DIR/chr$CHR.${M_STR}full \
                --thin-count 100000 \
                --make-bed \
                --out $SUBSET_DIR/chr$CHR.${M_STR}replicate$i.seed$SEED

            # Prune replicates for LD - Randomly select from pairs with r^2 > 0.01
            # using window of size 500 SNPS (stepsize 250 SNPs)
            plink \
                --bfile $SUBSET_DIR/chr$CHR.${M_STR}replicate$i.seed$SEED \
                --indep-pairwise 500 250 0.1 \
                --out $SUBSET_DIR/chr$CHR.${M_STR}replicate$i.seed$SEED.LD

            plink \
                --bfile $SUBSET_DIR/chr$CHR.${M_STR}replicate$i.seed$SEED \
                --exclude $SUBSET_DIR/chr$CHR.${M_STR}replicate$i.seed$SEED.LD.prune.out \
                --make-bed \
                --out $SUBSET_DIR/chr$CHR.${M_STR}replicate$i.seed$SEED.LD
        done

        # Combine 3L and 3R into single BED
        TMP_BED_LIST=`mktemp tmp.chr_bed_list.XXXXX`
        ls $SUBSET_DIR/chr3L.${M_STR}replicate$i.seed$SEED.LD.bed | \
            sed -e "s/\(.*\).bed/\1.bed \1.bim \1.fam/" >  $TMP_BED_LIST
        ls $SUBSET_DIR/chr3R.${M_STR}replicate$i.seed$SEED.LD.bed | \
            sed -e "s/\(.*\).bed/\1.bed \1.bim \1.fam/" >> $TMP_BED_LIST
        # Write combined BED file
        plink --merge-list $TMP_BED_LIST \
            --make-bed \
            --out $SUBSET_DIR/chr3.${M_STR}replicate$i.LD
        rm $TMP_BED_LIST
    done
done

rm $TMP_RANGE_PRE.3L
rm $TMP_RANGE_PRE.3R
rm $TMP_RANGE_PRE.3L.included.*
rm $TMP_RANGE_PRE.3R.included.*
