#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Prepare to create bootstraps for dadi
# ----------------------------------------------------------------------------------------

module load samtools
module load bedtools

NUM_BS=1000

# --- Bring in three pop lists
# --- These filenames determines the output dadi filename

BS_DIR=data/dadi_bootstraps
mkdir -p $BS_DIR

cp data/dadi/dadi_pop_list.species.txt             $BS_DIR
cp data/dadi/dadi_pop_list.island.txt              $BS_DIR
cp data/dadi/dadi_pop_list.mainlandisland.txt      $BS_DIR
cp data/dadi/dadi_pop_list.mainlandsmallisland.txt $BS_DIR

# ----------------------------------------------------------------------------------------

RAND_BED=data/random_hunks_for_dadi_bootstrap.bed

# --- Temporarily download Agam genome
AGAM_URL=https://www.vectorbase.org/sites/default/files/ftp/downloads/
AGAM_URL=${AGAM_URL}Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa.gz

AGAM_FA=$BS_DIR/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa
wget $AGAM_URL -O $AGAM_FA.gz
gunzip -f ${AGAM_FA}.gz
touch $AGAM_FA

samtools faidx $AGAM_FA

# Figure out (chr3) hunks, excluding any that overlap heterochromatic, inversion regions
cat data/inversion_simple.bed data/heterochromatin.bed | \
    cut -f 1-3 > $BS_DIR/tmp_inv_het.bed
bedtools random -l 1000000 -n 1000 -g $AGAM_FA.fai | \
    grep "^3" | \
    cut -f 1-3 | \
    bedtools subtract -A -a stdin -b $BS_DIR/tmp_inv_het.bed | \
    head -n $NUM_BS > $RAND_BED
rm $BS_DIR/tmp_inv_het.bed
