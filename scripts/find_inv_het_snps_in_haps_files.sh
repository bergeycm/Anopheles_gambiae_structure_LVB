#!/bin/sh

# ========================================================================================
# --- Find SNPs to be removed from het or inversion regions in haps files
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Input files
# ----------------------------------------------------------------------------------------

CHR=$1
IN_HAP=data/chr${CHR}.pass.snp.phased.haps

OUT_SNP_LIST=data/chr${CHR}.pass.snp.phased.haps.excluded.snps
echo -n > $OUT_SNP_LIST

cat data/inversion_simple.bed data/heterochromatin.bed | grep "^$CHR" | \
    while read line; do

        START=`echo $line | cut -d' ' -f 2`
        END=`  echo $line | cut -d' ' -f 3`

        echo "Adding SNPs in range ${START}-$END to exclude list..."

        awk -v start="$START" -v end="$END" \
            '{ if ($3 >= start && $3 <= end) print $0}' $IN_HAP >> $OUT_SNP_LIST

    done
