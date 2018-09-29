#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Make input file for CLUMPAK (zipped *.Q ADMIXTURE files)
# ----------------------------------------------------------------------------------------

TMP_DIR=data/ssese_with_ag1000g/adm_subsets/tmp.justQs

mkdir -p $TMP_DIR

cp data/ssese_with_ag1000g/adm_subsets/*.Q $TMP_DIR/

cd $TMP_DIR

zip chr3.allQs.zip       `ls chr3.*.Q | grep -v "withM"`
zip chr3.withM.allQs.zip `ls chr3.*.Q | grep    "withM"`

mv chr3.allQs.zip ..
mv chr3.withM.allQs.zip ..

cd ../../../..

rm -r $TMP_DIR
