#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Add dadi outgroup base in parallel (polarize)
# ----------------------------------------------------------------------------------------

IN_DADI=$1    # e.g. data/dadi/dadi.chr3.island.data

module load parallel
module load tabix

mkdir -p data/dadi_polarized

# Split input file into hunks of 200k lines
split -d -a 3 -l 200000 $IN_DADI ${IN_DADI}_tmp_pt

ls ${IN_DADI}_tmp_pt* > ${IN_DADI}_tmp_filelist.txt

parallel --jobs 20 -a ${IN_DADI}_tmp_filelist.txt \
    perl scripts/add_dadi_outgroup.pl {} ">" {}.polarized

OUT_DADI=`echo $IN_DADI | sed -e "s:data/dadi:data/dadi_polarized:"`

head -n1 $IN_DADI > $OUT_DADI
cat `ls -v ${IN_DADI}_tmp_pt*.polarized` >> $OUT_DADI

rm ${IN_DADI}_tmp*

exit
