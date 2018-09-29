#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Submit all rehh jobs
# ----------------------------------------------------------------------------------------

POPS=($(find data -name "ssese.seqids.isAG-*" -type f -exec wc -l {} \; | \
    awk '$1 > 0 { print $2 }' | sed -e "s/.*isAG-\(.*\)\.txt/\1/"))

CHRS=(2L 2R 3L 3R X)

# 2.5 Mb * 2 is 5Mb analyzed
WID=2500000

for POP in ${POPS[*]}; do
    # GSTE:
    qsub -v CHR="3R",POP="${POP}",CENTER=28596971,WIDTH=${WID} pbs/compute_ehh_target.pbs
    sleep 1
    # VGSC:
    qsub -v CHR="2L",POP="${POP}",CENTER=2394887,WIDTH=${WID} pbs/compute_ehh_target.pbs
    sleep 1
    # CYP6:
    qsub -v CHR="2R",POP="${POP}",CENTER=28493196,WIDTH=${WID} pbs/compute_ehh_target.pbs
    sleep 1
    # CYP9K1:
    qsub -v CHR="X",POP="${POP}",CENTER=15241718,WIDTH=${WID} pbs/compute_ehh_target.pbs
    sleep 1
    # Peak in X around 9 Mb:
    qsub -v CHR="X",POP="${POP}",CENTER=9000000,WIDTH=${WID} pbs/compute_ehh_target.pbs
    sleep 1

    # Submit full chromosome jobs
    for CHR in ${CHRS[*]}; do
        qsub -v CHR="X",POP="${POP}" pbs/compute_ehh_all.pbs
        sleep 1
    done
done

exit
