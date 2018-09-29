#!/bin/bash

# ========================================================================================
# --- Make stairway plot
# ========================================================================================

# Before run
# module load python/2.7.8
# export PYTHONPATH=/storage/home/cxb585/local_python
# export PYTHONPATH=$PYTHONPATH:/storage/home/cxb585/bin/dadi:/usr/global/python/2.7.8/lib/python2.7/site-packages

POPID=$1   # E.g. island.SSERINYA or mainlandisland.ISLAND or mainlandsmallisland.MAINLAND

module load parallel/20170422

# ----------------------------------------------------------------------------------------
# --- Make population list
# ----------------------------------------------------------------------------------------

SS_IND_INFO=data/ssese_individual_info_simple_bugala_split.txt

mkdir -p data/stairway

POP_FILE=data/stairway/stairway_pop_list.ssese.$POPID.txt

# island mode
if [[ $POPID =~ ^island. ]]; then
    awk '{ if ($6 == "AG") print $2,$3 }' $SS_IND_INFO | \
        sed -e "s/ /\t/" > $POP_FILE
fi

# mainlandisland mode
if [[ $POPID =~ "mainlandisland." ]]; then
    awk '{ if ($6 == "AG") print $2,$3 }' $SS_IND_INFO | \
        sed -e "s/BUWAMA/_MAINLAND/"  -e "s/KAZZI/_MAINLAND/" \
            -e "s/KIYINDI/_MAINLAND/" -e "s/MITYANA/_MAINLAND/" | \
        sed -e "s/ [^_].*/ ISLAND/" | sed -e "s/_//" | \
        sed -e "s/ /\t/" > $POP_FILE
fi

# mainlandsmallisland mode
if [[ $POPID =~ "mainlandsmallisland." ]]; then
    awk '{ if ($6 == "AG") print $2,$3 }' $SS_IND_INFO | \
        sed -e "s/BUWAMA/_MAINLAND/"  -e "s/KAZZI/_MAINLAND/" \
            -e "s/KIYINDI/_MAINLAND/" -e "s/MITYANA/_MAINLAND/" \
            -e "s/BUGALA/_BUGALA/" | \
        sed -e "s/ [^_].*/ ISLAND/" | sed -e "s/_//" | \
        sed -e "s/ /\t/" > $POP_FILE
fi

# Remove excluded animals (just for simplicity's sake) and sort
grep -v "ssese99" $POP_FILE | \
    sort -k2 > \
    ${POP_FILE/.txt/.flt.txt}

# ----------------------------------------------------------------------------------------
# --- Make input file for Stairway plot program
# ----------------------------------------------------------------------------------------

POPID_SIMP=`echo $POPID | sed -e "s/.*\.//"`
POP_COUNT=`grep $POPID_SIMP ${POP_FILE/.txt/.flt.txt} | wc -l`

# Grab dadi-outputted SFS
POP_SFS=`cat data/dadi/dadi.chr3.$POPID.sfs.txt`

# Set random breakpoints to try
# Roughly (nseq-2)/4, (nseq-2)/2, (nseq-2)*3/4, nseq-2
NSEQ=$((POP_COUNT * 2))
BREAK1=`echo "($NSEQ - 2) / 4" | bc`
BREAK2=`echo "($NSEQ - 2) / 2" | bc`
BREAK3=`echo "($NSEQ - 2) * 3/4" | bc`
BREAK4=`echo "($NSEQ - 2)" | bc`

# Mutation rate per year from Miles et al
PARAM_mu="3.5e-9"

# Generation time, 11 per year
PARAM_g=`echo "1/11" | bc -l`

# Length of genome
#   3L length = 41,963,435
#   3R length = 53,200,684
# Heterochromatic regions excluded:
#   3L           1   1815119 = 1,815,119 bp
#   3R    52161877  53200684 = 1,038,808 bp
#   3L     4264713   5031692 =  766,980 bp
#   3R    38988757  41860198 = 2,871,442 bp
# Total length: 41963435 + 53200684 - (1815119 + 1038808 + 766980 + 2871442)
#               = 88,671,770
GENOME_LEN=88671770

mkdir -p data/stairway_plot_$POPID/plots

cat /storage/work/cxb585/bin/stairway_plot_v2beta2/two-epoch_fold.blueprint | \
    sed -e "s/popid: .*/popid: $POPID/" | \
    sed -e "s/nseq: .*/nseq: $NSEQ/" | \
    sed -e "s/L: .*/L: $GENOME_LEN/" | \
    sed -e "s/SFS: .*/SFS: $POP_SFS/" | \
    sed -e "s/nrand: .*/nrand: $BREAK1\t$BREAK2\t$BREAK3\t$BREAK4/" | \
    sed -e "s/project_dir: .*/project_dir: data\/stairway_plot_$POPID/" | \
    sed -e "s/project_dir: .*/project_dir: data\/stairway_plot_$POPID\/plots/" | \
    sed -e "s/mu: .*/mu: $PARAM_mu/" | \
    sed -e "s/year_per_generation: .*/year_per_generation: $PARAM_g/" | \
    sed -e "s/plot_title: .*/plot_title: Stairway plot $POPID/" | \
    sed -e "s/xrange: .*/xrange: 0,15/" > \
    tmp.stairway_plot.$POPID.config

# ----------------------------------------------------------------------------------------
# --- Make stairway plot
# ----------------------------------------------------------------------------------------

# Bring in stairway plot directory
if [ ! -d "stairway_plot_es" ]; then
    cp -r /storage/work/cxb585/bin/stairway_plot_v2beta2/stairway_plot_es/ .
else
    # Hack to wait in case another script is mid-copying
    sleep 10
fi

# Create bash file to do analysis
java -cp stairway_plot_es Stairbuilder tmp.stairway_plot.$POPID.config

# Split bash file into parallelizable and non-parallelizable parts
STEP2_LINE=`grep -n "# Step 2" tmp.stairway_plot.${POPID}.config.sh | cut -d":" -f1`

# Parallel portion:
head -n $((STEP2_LINE - 1)) tmp.stairway_plot.$POPID.config.sh > \
    tmp.stairway_plot.$POPID.config.parallel.sh

# Serial portion:
sed 1,${STEP2_LINE}s/^/#/ tmp.stairway_plot.$POPID.config.sh > \
    tmp.stairway_plot.$POPID.config.serial.sh

# Run analysis

# Step 1 in parallel:
parallel -a tmp.stairway_plot.$POPID.config.parallel.sh --jobs 12

# Later steps serially:
mkdir -p data/stairway_plot_$POPID/plots/$BREAK1
mkdir -p data/stairway_plot_$POPID/plots/$BREAK2
mkdir -p data/stairway_plot_$POPID/plots/$BREAK3
mkdir -p data/stairway_plot_$POPID/plots/$BREAK4

bash tmp.stairway_plot.$POPID.config.serial.sh

# ----------------------------------------------------------------------------------------
# --- Clean up
# ----------------------------------------------------------------------------------------

# Wait because of file system latency
sleep 180

cp data/stairway_plot_$POPID/plots/*$POPID.final.summary.pdf \
    data/stairway_plot_$POPID/$POPID.final.summary.pdf

rm $POP_FILE
