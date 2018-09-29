#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Also create subsats of Ag1000G input file
# ----------------------------------------------------------------------------------------

module load parallel
module load plink
module load tabix

BS_DIR=data/dadi_bootstraps

ORIG_POP_FILE=data/dadi/dadi_pop_list.ag1000g.txt
RAND_BED=data/random_hunks_for_dadi_bootstrap.bed
AGAM_FA=$BS_DIR/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa

function convert_hunk {

    LINE=$1

    RANDOM_COUNT=`echo $LINE | cut -d ' ' -f 1`
    # Ag1000G bootstraps are 1-indexed, while the rest are 0-indexed.
    # Use this to fix them
    RANDOM_COUNT=$((RANDOM_COUNT - 1))
    p=`echo $LINE | cut -d ' ' -f 2-4`

    echo -n "--- Processing region number $RANDOM_COUNT for Ag1000G"

    TMP_BED=`mktemp`
    echo $p | cut -f1-3 -d' ' | tr ' ' '\t' | sed -e "s/$/\tbs_$RANDOM_COUNT.ag1000g/" | \
        sed -e "s/2L/1/" -e "s/2R/2/" -e "s/3L/3/" -e "s/3R/4/" > $TMP_BED

    TMP_VCF_PRE=$BS_DIR/bs_${RANDOM_COUNT}.ag1000g
    TMP_VCF=$TMP_VCF_PRE.vcf

    plink --file data/ag1000g.phase1.ar3/ag1000g.phase1.ar3.pass.biallelic.all.for.dadi \
        --extract range $TMP_BED \
        --recode vcf-iid \
        --out $TMP_VCF_PRE &> /dev/null

    # Fix chromosome numbering
    # [move from PLINK format (e.g. 1) to standard format (e.g. 2L)]
    sed -i'.bak' -e "s/^1\t/2L\t/" -e "s/^2\t/2R\t/" \
        -e "s/^3\t/3L\t/" -e "s/^4\t/3R\t/" $TMP_VCF
    rm ${TMP_VCF}.bak

    # Convert to dadi mode

    POP_FILE=$BS_DIR/dadi_pop_list_ag1000G.bs_${RANDOM_COUNT}.txt
    cp $ORIG_POP_FILE $POP_FILE

    perl scripts/convert_vcf_to_dadi_input.pl \
        $AGAM_FA \
        $TMP_VCF \
        $POP_FILE &> /dev/null

    # Rename dadi file and polarize
    head -n1 $POP_FILE.data > ${TMP_VCF_PRE}.data
    perl scripts/add_dadi_outgroup.pl $POP_FILE.data >> ${TMP_VCF_PRE}.data

    rm $TMP_BED
    rm $TMP_VCF
    rm ${TMP_VCF_PRE}.nosex
    rm ${TMP_VCF_PRE}.log
    rm $POP_FILE
    rm $POP_FILE.data

    RANDOM_COUNT=$(($RANDOM_COUNT + 1))

}

export BS_DIR
export ORIG_POP_FILE
export AGAM_FA

export -f convert_hunk

# Add iteration number to random hunk BED file
TMP_RAND_BED=`mktemp`
cat -n $RAND_BED | sed -e "s/^ \+//" | tr "\t" " " | head -n 250 > $TMP_RAND_BED

parallel \
    --progress \
    --env $BS_DIR --env $ORIG_POP_FILE --env $AGAM_FA \
    --jobs 8 \
    -a $TMP_RAND_BED \
    convert_hunk

rm $TMP_RAND_BED

exit
