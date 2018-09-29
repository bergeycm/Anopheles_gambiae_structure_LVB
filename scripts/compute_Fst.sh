#!/bin/sh

# ========================================================================================
# --- Compute Fst between islands
# ========================================================================================

IN_VCF_GZ=$1
OUT_PREFIX=$2

POP_A=$3
POP_B=$4

# Compute site-based Fst?
DO_SITE_FST=$5

# Compute windowed Fst?
DO_WIN_FST=$6

# Window size and step
WINDOW_SIZE=$7

echo "Input VCF:       $IN_VCF_GZ"
echo "Output prefix:   $OUT_PREFIX"
echo "Pop. A:          $POP_A"
echo "Pop. B:          $POP_B"
echo "Do site Fst:     $DO_SITE_FST"
echo "Do windowed Fst: $DO_WIN_FST"
echo "Window size:     $WINDOW_SIZE"

# ----------------------------------------------------------------------------------------
# --- Compute Fst between sites/populations (within species)
# ----------------------------------------------------------------------------------------

CMD="vcftools --gzvcf $IN_VCF_GZ --maf 0.05 \
        --out ${OUT_PREFIX} \
        --weir-fst-pop data/ssese.seqids.is-${POP_A}.txt \
        --weir-fst-pop data/ssese.seqids.is-${POP_B}.txt"

if [[ $DO_SITE_FST == 1 ]]; then
    $CMD
fi

if [[ $DO_WIN_FST == 1 ]]; then
    $CMD --fst-window-size $WINDOW_SIZE --fst-window-step $WINDOW_SIZE
fi

exit
