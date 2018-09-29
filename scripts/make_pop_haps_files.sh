#!/bin/sh

# ========================================================================================
# --- Split hap file by population
# ========================================================================================

# Modified from script of same name in repo:
# https://github.com/bergeycm/batwa-bakiga-exomes

# ----------------------------------------------------------------------------------------
# --- Input files
# ----------------------------------------------------------------------------------------

chr=$1

# --- Lists of individuals to include
pop=$2
# Example: data/ssese.seqids.is-BANDA.txt

# --- Keys for naming pops in output filename
pop_ID=$3
# Example: BANDA

IN_HAP=data/nohetinv/chr${chr}.pass.snp.phased.nohetinv.haps
SAMP_LIST=data/nohetinv/chr${chr}.pass.snp.phased.nohetinv.sample

# Figure out line numbers for samples we want to include
POP_INDS=(`grep -n -f $pop $SAMP_LIST | cut -d':' -f1`)

# Subtract 3 to get to zero-indexed position in list of samples
# (so line 3 represents the 0th sample in the list)
# Then figure out nth pair of columns that goes with this nth individual, now 1-indexed
# But then add 5 to get to the 1-based column number in the hap file
# (E.g. 0th column is actually the sixth column counting starting with 1 in the haps file)
POP_COLS=()
for (( i = 0 ; i < ${#POP_INDS[@]} ; i++ )); do
	# Get zero-indexed position in sample list
    INDEX_IN_LIST=$((${POP_INDS[$i]} - 3))
    POP_COL1=$(($(($(($INDEX_IN_LIST + 1)) * 2)) - 1))
    POP_COL2=$(($POP_COL1 + 1))
    POP_COL1=$(($POP_COL1 + 5))
    POP_COL2=$(($POP_COL2 + 5))
    POP_COLS=( "${POP_COLS[@]}" $POP_COL1 $POP_COL2 )
done

# Smoosh values into comma-separated string
POP_COL_STR=$(printf ",%s" "${POP_COLS[@]}")
POP_COL_STR=${POP_COL_STR:1}

# Reduce haps file to only columns we want
OUT_HAPS=data/`basename ${IN_HAP/.nohetinv.haps/.}${pop_ID}.haps`
cut -d' ' -f 1-5,$POP_COL_STR $IN_HAP > $OUT_HAPS
