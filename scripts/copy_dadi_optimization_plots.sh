#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Make copy of best dadi iteration's optimization plots
# ----------------------------------------------------------------------------------------

# --- 1-pop optmization plots

while read line; do
    PREFIX=`echo $line | sed -e "s/results/data/" \
                             -e "s/dadi.island-/dadi.chr3.island./" \
                             -e "s/.out.fix//"`
    ORIG_ITER=${PREFIX/iter/bestmodel.iter}.pdf
    BEST_ITER=`echo $ORIG_ITER | sed -e "s/iter.*/best-iter.pdf/"`

    cp $ORIG_ITER $BEST_ITER

done < results/dadi.best_iterations.txt

# --- 2-pop optmization plots

while read line; do
    PREFIX=`echo $line | sed -e "s/results/data/" \
                             -e "s/dadi.island.2pop/dadi.chr3.island/" \
                             -e "s/.out.fix//"`
    PREFIX=`echo $PREFIX | sed -e "s/island\.\(.*\)\.\(.*\).iter/island.\\1-\\2.iter/"`
    ORIG_ITER=$PREFIX.cmp.2d_fs.pdf
    BEST_ITER=`echo $ORIG_ITER | sed -e "s/iter.*/best-iter.pdf/"`

    cp $ORIG_ITER $BEST_ITER

done < results/dadi.best_iterations.2pop.txt

exit
