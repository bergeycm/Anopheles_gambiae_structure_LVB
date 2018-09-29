#!/bin/bash

# ========================================================================================
# --- Infer demography from IBS for two populations
# ========================================================================================

INFERRING_PATH=/storage/home/cxb585/work/bin/Inferring-demography-from-IBS/

POP_A=$1
POP_B=$2

module unload python
module load python/2.7.8

# ----------------------------------------------------------------------------------------
# --- Extract IBS tract length distributions
# ----------------------------------------------------------------------------------------

python $INFERRING_PATH/parse_between_pops_allpairs.py \
    IBS.$POP_A IBS.$POP_B 'None' 1

# ----------------------------------------------------------------------------------------

MIN_TRACT_LENGTH=5
NUM_REPS=10
GEN_TIME=1    # This can't be a float, so round up from 0.11 to 1.
              # Consider this during interpretation of results!

# Infer a simple population divergence (with no migration)

python $INFERRING_PATH/infer_onepop_adaptive_cumulative.py \
    results/IBS.${POP_A}_vs_IBS.${POP_B}_lengths.ibs \
    $MIN_TRACT_LENGTH $NUM_REPS $GEN_TIME \
    > results/IBS_inferred_size_${POP_A}_vs_${POP_B}.txt

DEMOG_HIST=results/IBS_inferred_size_${POP_A}_vs_${POP_B}.demographic_history.txt

grep -A1 "Parameters for plotting" results/IBS_inferred_size_${POP_A}_vs_${POP_B}.txt | \
    grep "^\[" > $DEMOG_HIST

python $INFERRING_PATH/plot_between.py \
    results/IBS.${POP_A}_vs_IBS.${POP_B}_lengths.ibs \
    $DEMOG_HIST $DEMOG_HIST.pdf
