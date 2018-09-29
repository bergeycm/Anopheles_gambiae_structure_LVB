#!/bin/sh

# ------------------------------------------------------------------------------
# --- Run ADMIXTURE
# ------------------------------------------------------------------------------

IN_BED=$1
K=$2
THREADS=${3:-1}

module load admixture/1.3.0

admixture -B200 -j${THREADS} --cv ${IN_BED} ${K}

OUT_PREFIX=`basename $IN_BED | sed -e "s/chr3\(.*\)\.withM.bed/chr3.withM\\1/" | \
    sed -e "s/\.bed$//" -e "s/\.thinned_for_ADMIXTURE//"`

for END in P Q Q_se Q_bias; do
    mv `basename ${IN_BED%.bed}`.${K}.$END \
        data/ssese_with_ag1000g/$OUT_PREFIX.${K}.$END
done

exit;
