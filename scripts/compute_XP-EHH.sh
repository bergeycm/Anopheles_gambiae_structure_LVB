#!/bin/sh

# ========================================================================================
# --- Compute XP-EHH with selScan
# ========================================================================================

# Modified from script of same name in repo:
# https://github.com/bergeycm/batwa-bakiga-exomes

module load gcc
module load selscan

# ----------------------------------------------------------------------------------------
# --- Input files
# ----------------------------------------------------------------------------------------

chr=$1
popA_ID=$2
popB_ID=$3

IN_HAP_A=data/chr${chr}.pass.snp.phased.${popA_ID}.haps
IN_HAP_B=data/chr${chr}.pass.snp.phased.${popB_ID}.haps

mkdir -p results/selscan
OUT_FILE=results/selscan/xp-ehh.$popA_ID-$popB_ID.$chr

# Make temporary directory to hold temporary files
TMP_DIR=`mktemp -d -p results/selscan/`

# ------------------------------------------------------------------------------------
# --- Reformat input data
# ------------------------------------------------------------------------------------

# Remove info columns from hap file
cut -d' ' -f 6- $IN_HAP_A > $TMP_DIR/`basename $IN_HAP_A`.fix
cut -d' ' -f 6- $IN_HAP_B > $TMP_DIR/`basename $IN_HAP_B`.fix

# ------------------------------------------------------------------------------------
# --- Transpose haps input file
# ------------------------------------------------------------------------------------

# awk script to transpose file taken from:
# http://stackoverflow.com/a/1729980

# transpose () {
#     awk '{
#         for (i=1; i<=NF; i++)  {
#             a[NR,i] = $i
#         }
#     }
#     NF>p { p = NF }
#     END {
#         for(j=1; j<=p; j++) {
#             str=a[1,j]
#             for(i=2; i<=NR; i++){
#                 str=str" "a[i,j];
#             }
#             print str
#         }
#     }' $1 > $1.transpose
# }

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

echo "Transposing first file..."
time transpose $TMP_DIR/`basename $IN_HAP_A`.fix
echo "Transposing second file..."
time transpose $TMP_DIR/`basename $IN_HAP_B`.fix
echo "Finished transposing files."

# ----------------------------------------------------------------------------------------
# --- Fake map file
# ----------------------------------------------------------------------------------------

MAP_FILE=$TMP_DIR/`basename $OUT_FILE`.map

# Since we lack map for Anopheles, assume a constant recombination rate of 2.0 cM/Mb
awk '{ printf "%s %s %.8f %s\n", $1,$2,$3/500000,$3 }' $IN_HAP_A > \
    $MAP_FILE

# ------------------------------------------------------------------------------------
# --- Run selscan
# ------------------------------------------------------------------------------------

CMD="selscan --xpehh \
    --threads 8 \
    --hap $TMP_DIR/`basename $IN_HAP_A`.fix.transpose \
    --ref $TMP_DIR/`basename $IN_HAP_B`.fix.transpose \
    --map $MAP_FILE \
    --out $OUT_FILE"

echo "CMD: [$CMD]"

`$CMD`

# Delete temporary files
rm -r $TMP_DIR
