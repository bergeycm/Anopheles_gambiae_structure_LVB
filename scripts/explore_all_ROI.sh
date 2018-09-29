#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Prepare for local inference around regions of interest
# ----------------------------------------------------------------------------------------

# Peaks can be found in the result of:
# Rscript scripts/find_peaks_in_xp-ehh.R > results/xpehh_peaks.txt

#module load r/3.4
module load vcftools
module load plink
module load python/2.7.14-anaconda5.0.1
module load parallel

# ----------------------------------------------------------------------------------------

# Start file of commands
DEND_CMDS=tmp.dendro.cmds.sh
echo > $DEND_CMDS

# ----------------------------------------------------------------------------------------

echo "Bump around 34Mb on 2L..."
CHR=2L
CENTER=`grep "2L_near34Mb" results/xpehh_peaks.txt | cut -d":" -f2`
NAME=2L_34Mb
for WINDOW in 10000 100000; do
    for DIP_HAP in FALSE; do             # Only FALSE now
        for W_AG in TRUE; do             # Only TRUE now
            for DO_SNPEFF in FALSE; do   # Only FALSE now
                CMD="Rscript scripts/explore_around_ROI.R $CHR $CENTER $WINDOW $NAME \
                    $DIP_HAP $W_AG $DO_SNPEFF \
                    &> reports/dendro.logs.$NAME.$WINDOW.$DIP_HAP.$W_AG.log"
                echo $CMD >> $DEND_CMDS
            done
        done
    done
done

# ----------------------------------------------------------------------------------------

echo "CYP cluster on 2R..."
CHR=2R
CENTER=28501972
NAME=2R_CYP6P2
for WINDOW in 10000 100000; do
    for DIP_HAP in FALSE; do             # Only FALSE now
        for W_AG in TRUE; do             # Only TRUE now
            for DO_SNPEFF in FALSE; do   # Only FALSE now
                CMD="Rscript scripts/explore_around_ROI.R $CHR $CENTER $WINDOW $NAME \
                    $DIP_HAP $W_AG $DO_SNPEFF \
                    &> reports/dendro.logs.$NAME.$WINDOW.$DIP_HAP.$W_AG.log"
                echo $CMD >> $DEND_CMDS
            done
        done
    done
done

# ----------------------------------------------------------------------------------------

#   echo "CYP cluster on 3L..."
#   CHR=3L
#   CENTER=`grep "3L_CYP" results/xpehh_peaks.txt | cut -d":" -f2`
#   NAME=3L_CYP
#   for WINDOW in 10000 100000; do
#       for DIP_HAP in FALSE; do             # Only FALSE now
#           for W_AG in TRUE; do             # Only TRUE now
#               for DO_SNPEFF in FALSE; do   # Only FALSE now
#                   Rscript scripts/explore_around_ROI.R $CHR $CENTER $WINDOW $NAME \
#                       $DIP_HAP $W_AG $DO_SNPEFF \
#                       &> reports/dendro.logs.$NAME.$WINDOW.$DIP_HAP.$W_AG.log &
#               done
#           done
#       done
#   done

# ----------------------------------------------------------------------------------------

echo "Bump around 9Mb on X..."
CHR=X
CENTER=`grep "X_near9Mb" results/xpehh_peaks.txt | cut -d":" -f2`
NAME=X_9Mb
for WINDOW in 10000 100000; do
    for DIP_HAP in FALSE; do             # Only FALSE now
        for W_AG in TRUE; do             # Only TRUE now
            for DO_SNPEFF in FALSE; do   # Only FALSE now
                CMD="Rscript scripts/explore_around_ROI.R $CHR $CENTER $WINDOW $NAME \
                    $DIP_HAP $W_AG $DO_SNPEFF \
                    &> reports/dendro.logs.$NAME.$WINDOW.$DIP_HAP.$W_AG.log"
                echo $CMD >> $DEND_CMDS
            done
        done
    done
done

# ----------------------------------------------------------------------------------------

echo "CYP9K1 on X..."
CHR=X
CENTER=15241718
NAME=X_CYP9K1
for WINDOW in 100000; do     # For CYP9K1, only larger window is done
    for DIP_HAP in FALSE; do             # Only FALSE now
        for W_AG in TRUE; do             # Only TRUE now
            for DO_SNPEFF in FALSE; do   # Only FALSE now
                CMD="Rscript scripts/explore_around_ROI.R $CHR $CENTER $WINDOW $NAME \
                    $DIP_HAP $W_AG $DO_SNPEFF \
                    &> reports/dendro.logs.$NAME.$WINDOW.$DIP_HAP.$W_AG.log"
                echo $CMD >> $DEND_CMDS
            done
        done
    done
done

# Make empty files for CYP9K1 with the small window
CYP9K1_PREFIX=results/for_dendrograms
CYP9K1_PREFIX=$CYP9K1_PREFIX/chrX.X_CYP9K1.10000.haploid.with-ag1000g.sans-snpeff.haps
touch $CYP9K1_PREFIX.fancydendro.pdf
touch $CYP9K1_PREFIX.haplotypes.txt

# ----------------------------------------------------------------------------------------

echo "Random region at 4Mb on X..."
CHR=X
CENTER=4000000
NAME=X_4Mb
for WINDOW in 10000 100000; do
    for DIP_HAP in FALSE; do             # Only FALSE now
        for W_AG in TRUE; do             # Only TRUE now
            for DO_SNPEFF in FALSE; do   # Only FALSE now
                CMD="Rscript scripts/explore_around_ROI.R $CHR $CENTER $WINDOW $NAME \
                    $DIP_HAP $W_AG $DO_SNPEFF \
                    &> reports/dendro.logs.$NAME.$WINDOW.$DIP_HAP.$W_AG.log"
                echo $CMD >> $DEND_CMDS
            done
        done
    done
done

# ----------------------------------------------------------------------------------------

echo "GSTE on 3R..."
CHR=3R
CENTER=28598038
NAME=3R_GSTE
for WINDOW in 10000 100000; do
    for DIP_HAP in FALSE; do             # Only FALSE now
        for W_AG in TRUE; do             # Only TRUE now
            for DO_SNPEFF in FALSE; do   # Only FALSE now
                CMD="Rscript scripts/explore_around_ROI.R $CHR $CENTER $WINDOW $NAME \
                    $DIP_HAP $W_AG $DO_SNPEFF \
                    &> reports/dendro.logs.$NAME.$WINDOW.$DIP_HAP.$W_AG.log"
                echo $CMD >> $DEND_CMDS
            done
        done
    done
done

# ----------------------------------------------------------------------------------------

parallel --progress -j 12 -a $DEND_CMDS
