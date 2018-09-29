#!/bin/sh

# ========================================================================================
# --- Bootstrap Fst between islands
# ========================================================================================

module load vcftools

IN_VCF_GZ=$1

POP_A=$2
POP_B=$3

# Number of iterations
NUM_ITER=$4

# Distance between SNPs to include in thinned subset
THIN_DIST=$5

echo "Input VCF:         $IN_VCF_GZ"
echo "Pop. A:            $POP_A"
echo "Pop. B:            $POP_B"
echo "Num. iterations:   $NUM_ITER"
echo "SNP thinning dist: $THIN_DIST"

#    # Test:
#    IN_VCF_GZ=data/chr3L.pass.snp.flt.vcf.gz
#    POP_A=BANDA
#    POP_B=BUGALAIS
#    NUM_ITER=1000
#    THIN_DIST=5000

# ----------------------------------------------------------------------------------------
# --- Compute Fst between sites/populations (within species) after shuffling groups
# ----------------------------------------------------------------------------------------

# Make folder to house temporary population files
THIS_POPLIST_DIR=temp_shuffled_pop_lists/${POP_A}-${POP_B}
rm -rf $THIS_POPLIST_DIR
mkdir -p $THIS_POPLIST_DIR

# Make folder for Fst results
CHR=`basename $IN_VCF_GZ | cut -d'.' -f1 | sed -e "s/chr//"`
BOOT_FST_DIR=results/bootstrapped_Fst
mkdir -p $BOOT_FST_DIR
OUT_PREFIX=$BOOT_FST_DIR/${CHR}.${POP_A}_${POP_B}

# Downsample VCF and reduce to just individuals in these pops
vcftools --gzvcf $IN_VCF_GZ \
         --keep data/ssese.seqids.is-${POP_A}.txt \
         --keep data/ssese.seqids.is-${POP_B}.txt \
         --maf 0.05 \
         --thin $THIN_DIST \
         --recode \
         --out $OUT_PREFIX

# Compute Fst with real population data
vcftools \
    --vcf $OUT_PREFIX.recode.vcf \
    --out $OUT_PREFIX.real \
    --weir-fst-pop data/ssese.seqids.is-${POP_A}.txt \
    --weir-fst-pop data/ssese.seqids.is-${POP_B}.txt 2> $OUT_PREFIX.real.log

for ITER in `seq 1 $NUM_ITER`; do

    echo "Iteration $ITER..."

    # --- Shuffle groups
    POP_BOTH_TMP=$(mktemp $THIS_POPLIST_DIR/POP_BOTH.XXXXXX)
    POP_1_TMP=$(mktemp $THIS_POPLIST_DIR/POP_1.XXXXXX)
    POP_2_TMP=$(mktemp $THIS_POPLIST_DIR/POP_2.XXXXXX)

    cat data/ssese.seqids.is-${POP_A}.txt data/ssese.seqids.is-${POP_B}.txt | \
        shuf > $POP_BOTH_TMP

    head -n `wc -l data/ssese.seqids.is-${POP_A}.txt | cut -f1` $POP_BOTH_TMP > $POP_1_TMP
    tail -n `wc -l data/ssese.seqids.is-${POP_B}.txt | cut -f1` $POP_BOTH_TMP > $POP_2_TMP

    vcftools \
        --vcf $OUT_PREFIX.recode.vcf \
        --out $OUT_PREFIX.iter$ITER \
        --weir-fst-pop $POP_1_TMP \
        --weir-fst-pop $POP_2_TMP 2> $OUT_PREFIX.iter$ITER.log

done

# ----------------------------------------------------------------------------------------
# --- Compute p-value
# ----------------------------------------------------------------------------------------

REAL_FST=`grep "weighted" $OUT_PREFIX.real.log | \
    cut -d' ' -f7`
grep "weighted" $OUT_PREFIX.iter*log | cut -d' ' -f7 > $OUT_PREFIX.weighted.fst

# Get count of lines with higher Fst value
NUM_GREATER=`grep "weighted" $OUT_PREFIX.iter*log | cut -d' ' -f7 | \
    awk -v real="$REAL_FST" '{ if ($1 > real) print $0 }' | wc -l`

TOT_TRIALS=`ls $OUT_PREFIX.iter*log | wc -l`

echo "100 * $NUM_GREATER / $TOT_TRIALS" | bc -l

# ----------------------------------------------------------------------------------------
# --- Clean up
# ----------------------------------------------------------------------------------------

rm $OUT_PREFIX.iter*log
rm $OUT_PREFIX.iter*weir.fst
rm $OUT_PREFIX.real.weir.fst
rm $OUT_PREFIX.recode.vcf

rm -rf $THIS_POPLIST_DIR

# Output files:
#   $OUT_PREFIX.weighted.fst    # iteration Fst values
#   $OUT_PREFIX.real.log        # real Fst values (can be extracted from log)

exit
