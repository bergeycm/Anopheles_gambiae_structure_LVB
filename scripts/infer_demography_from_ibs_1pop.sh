#!/bin/bash

# ========================================================================================
# --- Infer demography from IBS for single population
# ========================================================================================

INFERRING_PATH=/storage/home/cxb585/work/bin/Inferring-demography-from-IBS/

POP=$1

module unload python
module load python/2.7.8

# ----------------------------------------------------------------------------------------
# --- Extract IBS tract length distributions
# ----------------------------------------------------------------------------------------

python $INFERRING_PATH/parse_within_pop_allpairs.py results/IBS.$POP 'None' 1

# ----------------------------------------------------------------------------------------
# --- Infer size changes within one population
# ----------------------------------------------------------------------------------------

MIN_TRACT_LENGTH=5
NUM_REPS=10
GEN_TIME=1    # This can't be a float, so round up from 0.11 to 1.
              # Consider this during interpretation of results!

python $INFERRING_PATH/infer_onepop_adaptive_cumulative.py \
    results/IBS.${POP}_lengths.ibs \
    $MIN_TRACT_LENGTH $NUM_REPS $GEN_TIME \
    > results/IBS_inferred_size_$POP.txt

DEMOG_HIST=results/IBS_inferred_size_${POP}.demographic_history.txt

grep -A1 "Parameters for plotting" results/IBS_inferred_size_$POP.txt | \
    grep "^\[" > $DEMOG_HIST

python $INFERRING_PATH/plot_onepop.py \
    results/IBS.${POP}_lengths.ibs $DEMOG_HIST $DEMOG_HIST.pdf
