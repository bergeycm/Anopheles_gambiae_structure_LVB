#!/bin/sh

# ------------------------------------------------------------------------------
# --- Run ADMIXTURE in cross-validation mode
# ------------------------------------------------------------------------------

IN_BED=$1
K=$2
THREADS=${3:-1}

module load admixture/1.3.0

OUT_PREFIX=${IN_BED/\.bed/.ADMIXTURE_log}

OUT=$OUT_PREFIX${K}.out

RANDOM=`date +%N | sed s/...$//`
SEED=$RANDOM

admixture -s $SEED -j${THREADS} --cv ${IN_BED} ${K} | tee $OUT

exit;
