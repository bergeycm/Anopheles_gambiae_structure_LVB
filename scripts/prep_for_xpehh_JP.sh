xpehh WIP

# ========================================================================================
# --- Prep to call XP-EHH program
# ========================================================================================

CHR=$1

XP_DIR=data/xpehh_input_jp
mkdir -p $XP_DIR

POPS=( BANDA BUGALAIS BUGALAML BUKASA BUWAMA KAZZI KIYINDI MITYANA NSADZI SSERINYA )

# ----------------------------------------------------------------------------------------
# --- Make input haplotype files
# ----------------------------------------------------------------------------------------

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

for POP in ${POPS[*]}; do

    awk '{for (i=6; i<NF; i++) printf $i " "; print $NF}' \
        data/chr$CHR.pass.snp.phased.$POP.haps > \
        $XP_DIR/chr$CHR.pass.snp.phased.$POP.jp.haps

    echo "Transposing input file..."
    time transpose $XP_DIR/chr$CHR.pass.snp.phased.$POP.jp.haps
    echo "Finished transposing file."

done

# ----------------------------------------------------------------------------------------
# --- Make input map file
# ----------------------------------------------------------------------------------------

awk '{ print $1"_"$3,$3,$3 / 2000000,$4,$5 }' \
    data/chr$CHR.pass.snp.phased.BANDA.haps > \
    $XP_DIR/chr$CHR.pass.snp.phased.jp.map

exit
