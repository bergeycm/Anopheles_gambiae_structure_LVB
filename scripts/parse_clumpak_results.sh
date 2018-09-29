#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Parse CLUMPAK output
# ----------------------------------------------------------------------------------------

# CLUMPAK results (unzipped) are put into the following directory as
# CLUMPAK_sansM/ and CLUMPAK_withM/
# Rename to lose stupid equals sign in folder name

ADM_DIR=data/ssese_with_ag1000g/adm_subsets/

for K in `seq 2 10`; do

    for M_FLAG in sansM withM; do

        IN=$ADM_DIR/CLUMPAK_$M_FLAG/K$K/MajorCluster/CLUMPP.files/ClumppIndFile.output
        OUT=$ADM_DIR/CLUMPAK_$M_FLAG.$K.Q

        cut -d":" -f2 $IN | sed -e "s/^ //" > $OUT

    done
done
