#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Convert to EIGENSTRAT format
# ----------------------------------------------------------------------------------------

module load plink/1.9
module load eigensoft

PREFIX=data/all.pass.snp.flt.eigen

# Create list of input BED files (by chr) and combine
ls data/chr*.pass.snp.flt.bed | \
    sed -e "s/\(.*\).bed/\\1.bed \\1.bim \\1.fam/" > $PREFIX.bedlist.txt

cat data/inversion_simple.set  > data/invhet.set
awk 'BEGIN {OFS="\t"} { print $1,$2,$3,"het" }' \
    data/heterochromatin.set >> data/invhet.set

plink --merge-list $PREFIX.bedlist.txt \
    --exclude range data/invhet.set \
    --make-bed --out $PREFIX

PAR=${PREFIX}.par.BED.EIGENSTRAT

# ----------------------------------------------------------------------------------------

# Rename
cp ${PREFIX}.fam ${PREFIX}.pedind

echo "genotypename:    ${PREFIX}.bed"    >  $PAR
echo "snpname:         ${PREFIX}.bim"    >> $PAR
echo "indivname:       ${PREFIX}.pedind" >> $PAR
echo "outputformat:    EIGENSTRAT"       >> $PAR

echo "genotypeoutname: ${PREFIX}.eigenstratgeno" >> $PAR
echo "snpoutname:      ${PREFIX}.snp"            >> $PAR
echo "indivoutname:    ${PREFIX}.ind"            >> $PAR

echo "familynames:     NO" >> $PAR
echo "numchrom:        3"  >> $PAR

# chrom:       Only output SNPs on this chromosome.
# lopos:       Only output SNPs with physical position >= this value.
# hipos:       Only output SNPs with physical position <= this value.
# echo "chrom:           1"        >> $PAR
# echo "lopos:           1"        >> $PAR
# echo "hipos:           10000000" >> $PAR

# ----------------------------------------------------------------------------------------

convertf -p $PAR

# ----------------------------------------------------------------------------------------
# --- Add sex and population information
# ----------------------------------------------------------------------------------------

perl scripts/add_pop_sex_eigen.pl ${PREFIX}.ind
