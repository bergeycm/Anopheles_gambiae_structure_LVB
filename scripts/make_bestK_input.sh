#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Make input file for Best K by Evanno
# ----------------------------------------------------------------------------------------

ADM_DIR=data/ssese_with_ag1000g/adm_subsets

IN_SANS_M=$ADM_DIR/chr3.replicate*.LD.iter*.*.ADMIXTURE_log*.out
IN_WITH_M=$ADM_DIR/chr3.withM.replicate*.LD.iter*.*.ADMIXTURE_log*.out

# --- S-only

grep "^Loglikelihood" $IN_SANS_M | \
    sed -e "s/.*log\([0-9]\+\).* \(.*\)/\1\t\2/" | \
    sort -n -k 1 > $ADM_DIR/bestK.input.log_prob.txt

# --- Both M and S

grep "^Loglikelihood" $IN_WITH_M | \
    sed -e "s/.*log\([0-9]\+\).* \(.*\)/\1\t\2/" | \
    sort -n -k 1 > $ADM_DIR/bestK.input.withM.log_prob.txt
