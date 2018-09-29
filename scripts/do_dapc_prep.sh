#!/bin/bash

# ========================================================================================
# --- Prepare to do DAPC on putatively selected SNPs
# ========================================================================================

module load plink

# ----------------------------------------------------------------------------------------
# --- Prep to do DAPC
# ----------------------------------------------------------------------------------------

# Temporary files and prefix (needed to not conflict with other jobs when writing
# temporary BED, BIM, and FAM files)
TMP_BED_LIST=`mktemp tmp.chr_bed_list.XXXXX`
TMP_INV_HET_SIMP=`mktemp tmp.inversion_het_simple.XXXXX`
TMP_EXCLUDE=`mktemp tmp.exclude.XXXXX`
TMP_PREFIX=`mktemp -u data/tmp.all.pass.snp.flt.noinv.XXXXX`

# Create list of input BED files (by chr) and simple inversions/het BED
ls data/chr*.pass.snp.flt.bed | \
    sed -e "s/\(.*\).bed/\\1.bed \\1.bim \\1.fam/" > $TMP_BED_LIST
cat data/inversion_simple.bed > $TMP_INV_HET_SIMP
awk 'BEGIN {{OFS="\t"}} {{ print $1,$2,$3,"het" }}' data/heterochromatin.bed \
    >> $TMP_INV_HET_SIMP
sed -e "s/^2L/1/" -e "s/^2R/2/" -e "s/^3L/3/" -e "s/^3R/4/" -iBACKUP $TMP_INV_HET_SIMP

# Make SNP set to exclude the inversions and het regions.
plink --merge-list $TMP_BED_LIST \
    --make-set $TMP_INV_HET_SIMP \
    --write-set --out data/for_selection.to_exclude.pass.snp.flt.pca.inversions

# Combine stuff to exclude: LD pruned stuff and inv/het regions
cat data/all.pass.snp.flt.prune.out \
    data/for_selection.to_exclude.pass.snp.flt.pca.inversions.set \
    > $TMP_EXCLUDE

# --- Write files of SNPs (with no Fst filtering)

# First command is just to make a MAP file by going through PED intermediary
plink --merge-list $TMP_BED_LIST \
    --exclude $TMP_EXCLUDE \
    --maf 0.01 \
    --recode \
    --out results/dapc_subset_snps.ALL

# Next makes raw format file
plink --file results/dapc_subset_snps.ALL \
    --maf 0.01 \
    --recode A \
    --out results/dapc_subset_snps.ALL

# --- Do the same (SNPs with no Fst filtering) but just chr3

# First command is just to make a MAP file by going through PED intermediary
plink --merge-list $TMP_BED_LIST \
    --exclude $TMP_EXCLUDE \
    --chr 3 \
    --maf 0.01 \
    --recode \
    --out results/dapc_subset_snps.chr3L

plink --merge-list $TMP_BED_LIST \
    --exclude $TMP_EXCLUDE \
    --chr 4 \
    --maf 0.01 \
    --recode \
    --out results/dapc_subset_snps.chr3R

# Next makes raw format file
plink --file results/dapc_subset_snps.chr3L \
    --maf 0.01 \
    --recode A \
    --out results/dapc_subset_snps.chr3L

plink --file results/dapc_subset_snps.chr3R \
    --maf 0.01 \
    --recode A \
    --out results/dapc_subset_snps.chr3R

# --- Clean up

rm $TMP_BED_LIST
rm $TMP_EXCLUDE
rm $TMP_INV_HET_SIMP
