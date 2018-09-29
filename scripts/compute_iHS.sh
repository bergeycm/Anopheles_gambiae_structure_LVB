#!/bin/sh

# ========================================================================================
# --- Compute iHS with selScan
# ========================================================================================

module load gcc

module load selscan

# ----------------------------------------------------------------------------------------
# --- Input files
# ----------------------------------------------------------------------------------------

chr=$1
pop_ID=$2

IN_HAP=data/chr${chr}.pass.snp.phased.${pop_ID}.haps

mkdir -p results/selscan
OUT_FILE=results/selscan/iHS.$pop_ID.$chr

# Make temporary directory to hold temporary files
TMP_DIR=`mktemp -d -p results/selscan/`

# ------------------------------------------------------------------------------------
# --- Reformat input data
# ------------------------------------------------------------------------------------

# Remove info columns from hap file
cut -d' ' -f 6- $IN_HAP > $TMP_DIR/`basename $IN_HAP`.fix

# ------------------------------------------------------------------------------------
# --- Transpose haps input file
# ------------------------------------------------------------------------------------

# awk script to transpose file modified from skeleton at:
# http://www.linuxquestions.org/questions/programming-9/awk-multiple-passes-794715/

transpose () {
    awk '
    NR == 1 {
        # Get number of lines in file
        ( "cat " FILENAME " | wc -l" ) | getline NL
    }
    { printf "%s%s", (FNR>1 ? OFS : ""), $ARGIND }
    FNR == NL {
        # Rerun for each column in file
        while ( ++count < NF ) {
            print ""
            ARGC++
            ARGV[ARGIND+1] = FILENAME
            nextfile
        }
    }
    END { print "" }' $1 > $1.transpose

}

echo "Transposing input file..."
time transpose $TMP_DIR/`basename $IN_HAP`.fix
echo "Finished transposing file."

# ----------------------------------------------------------------------------------------
# --- Fake map file
# ----------------------------------------------------------------------------------------

MAP_FILE=$TMP_DIR/`basename $OUT_FILE`.map

# Since we lack map for Anopheles, assume a constant recombination rate of 2.0 cM/Mb
awk '{ printf "%s %s %.8f %s\n", $1,$2,$3/500000,$3 }' $IN_HAP > \
    $MAP_FILE

# ------------------------------------------------------------------------------------
# --- Run selscan
# ------------------------------------------------------------------------------------

CMD="selscan --ihs \
    --threads 8 \
    --hap $TMP_DIR/`basename $IN_HAP`.fix.transpose \
    --map $MAP_FILE \
    --out $OUT_FILE"

echo "CMD: [$CMD]"

`$CMD`

# Delete temporary files
rm -r $TMP_DIR
