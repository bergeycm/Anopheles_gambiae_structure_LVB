#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Find animals with high missingness
# ----------------------------------------------------------------------------------------

imiss=reports/all.pass.snp.imiss

awk '{if ($5 > 0.1 && $1 != "INDV") print $1 }' $imiss | \
    sort | uniq > reports/inds.high.missing.txt
