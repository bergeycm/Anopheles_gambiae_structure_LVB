#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Compute XP-EHH
# ----------------------------------------------------------------------------------------

module load xpehh

XP_DIR=data/xpehh_input_jp
mkdir -p $XP_DIR

CHR=$1
POP_A=$2
POP_B=$3

xpehh -m $XP_DIR/chr$CHR.pass.snp.phased.jp.map \
      -h $XP_DIR/chr$CHR.pass.snp.phased.$POP_A.jp.haps.transpose \
         $XP_DIR/chr$CHR.pass.snp.phased.$POP_B.jp.haps.transpose > \
         results/xpehh_jp/xp-ehh.$POP_A-$POP_B.$CHR.xpehh.out

exit
