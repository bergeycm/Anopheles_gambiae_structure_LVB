#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Parse Hudson Fst output
# ----------------------------------------------------------------------------------------

HUDSON_OUT=data/all.pass.snp.flt.eigen.fst.se.out

Z_HEAD_LINE=`grep -n "^Fst: Z" $HUDSON_OUT | cut -d":" -f1`

FST_HEAD_LINE=`grep -n "^fst \*1000" $HUDSON_OUT | cut -d":" -f1`

SD_HEAD_LINE=`grep -n "^s.dev \* 1000000:" $HUDSON_OUT | cut -d":" -f1`

LAST_LINE=`grep -n "^##end of smartpca run" $HUDSON_OUT | cut -d":" -f1`

POPS=`grep "population" $HUDSON_OUT | sed -e "s/ \+/\t/g" | cut -f3`

# Grab Z-scores
echo ${POPS[*]} | tr " " "\t" | sed -e "s/^/\t/" > $HUDSON_OUT.fstZ.txt
sed -n $((Z_HEAD_LINE + 2)),$((FST_HEAD_LINE - 2))p $HUDSON_OUT | \
    sed -e "s/ \+/\t/g" -e "s/^\t//" >> $HUDSON_OUT.fstZ.txt

# Grab Fst estimates
echo ${POPS[*]} | tr " " "\t" | sed -e "s/^/\t/" > $HUDSON_OUT.fst.txt
sed -n $((FST_HEAD_LINE + 2)),$((SD_HEAD_LINE - 2))p $HUDSON_OUT | \
    sed -e "s/ \+/\t/g" -e "s/^\t//" >> $HUDSON_OUT.fst.txt

# Grab SD
echo ${POPS[*]} | tr " " "\t" | sed -e "s/^/\t/" > $HUDSON_OUT.sd.txt
sed -n $((SD_HEAD_LINE + 2)),$((LAST_LINE - 2))p $HUDSON_OUT | \
    sed -e "s/ \+/\t/g" -e "s/^\t//" >> $HUDSON_OUT.sd.txt

exit
