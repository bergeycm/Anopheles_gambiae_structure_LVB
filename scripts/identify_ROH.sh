#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Identify runs of homozygosity (ROH)
# ----------------------------------------------------------------------------------------

module load plink/1.90b3.40

# Temporary files and prefix (needed to not conflict with other jobs when writing
# temporary BED, BIM, and FAM files)
TMP_BED_LIST=`mktemp tmp.chr_bed_list.XXXXX`
TMP_INV_HET_SIMP=`mktemp tmp.inversion_simple.XXXXX`
TMP_EXCLUDE=`mktemp tmp.exclude.XXXXX`
TMP_PREFIX=`mktemp -u data/tmp.all.pass.snp.flt.noinv.XXXXX`

# Create list of input BED files (by chr) and simple inversions BED
ls data/chr*.pass.snp.flt.bed | \
    sed -e "s/\(.*\).bed/\\1.bed \\1.bim \\1.fam/" > $TMP_BED_LIST

cat data/inversion_simple.bed data/heterochromatin.bed | \
    sed -e "s/^2L/1/" -e "s/^2R/2/" -e "s/^3L/3/" -e "s/^3R/4/" > $TMP_INV_HET_SIMP

# Make SNP set covering the 2La and 2Rb inversions and het regions
plink --merge-list $TMP_BED_LIST \
    --exclude data/all.pass.snp.flt.prune.out \
    --make-set $TMP_INV_HET_SIMP --write-set \
    --out $TMP_PREFIX.inversions

cat data/all.pass.snp.flt.prune.out $TMP_PREFIX.inversions.set > $TMP_EXCLUDE

plink \
    --merge-list $TMP_BED_LIST \
    --exclude $TMP_EXCLUDE \
    --homozyg group-verbose \
    --homozyg-snp 10 --homozyg-kb 100 \
    --homozyg-window-het 3 --homozyg-window-missing 5 \
    --out $TMP_PREFIX

mkdir -p results/ROH
mv $TMP_PREFIX.hom         results/ROH/all.pass.snp.flt.noinv.hom
mv $TMP_PREFIX.hom.indiv   results/ROH/all.pass.snp.flt.noinv.hom.indiv
mv $TMP_PREFIX.hom.summary results/ROH/all.pass.snp.flt.noinv.hom.summary
# ROH pool report
mv $TMP_PREFIX.hom.overlap results/ROH/all.pass.snp.flt.noinv.hom.overlap
# Per-pool reports
for file in $TMP_PREFIX.hom.overlap.S*.verbose; do
    mv $file ${file/$TMP_PREFIX/results/ROH/all.pass.snp.flt.noinv}
done

# --- Clean up
rm $TMP_BED_LIST
rm $TMP_INV_HET_SIMP
rm $TMP_EXCLUDE
rm $TMP_PREFIX*
