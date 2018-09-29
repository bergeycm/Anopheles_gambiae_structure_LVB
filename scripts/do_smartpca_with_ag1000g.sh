#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Convert to EIGENSTRAT format for PCA with Ssese and Ag1000G
# ----------------------------------------------------------------------------------------

module load plink/1.9
###module load eigensoft
export PATH=$PATH:/storage/work/cxb585/bin/EIG/bin        
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/storage/work/cxb585/bin/OpenBLAS/build/include     
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/storage/work/cxb585/bin/OpenBLAS/build/lib/     

PREFIX=data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict

# Create list of input BED files (by chr) and combine
ls data/ssese_with_ag1000g/chr3*.pass.snp.phased.ag1000g.strict.bed | \
    sed -e "s/\(.*\).bed/\\1.bed \\1.bim \\1.fam/" > $PREFIX.bedlist.txt

cat data/inversion_simple.set  > $PREFIX.invhet.set
awk 'BEGIN {OFS="\t"} { print $1,$2,$3,"het" }' \
    data/heterochromatin.set >> $PREFIX.invhet.set

plink --merge-list $PREFIX.bedlist.txt \
    --exclude range $PREFIX.invhet.set \
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

# ----------------------------------------------------------------------------------------
# --- Call smartpca
# ----------------------------------------------------------------------------------------

sed -e "s/Burkina Faso/BF/g" -e "s/Guinea-Bissau/GW/g" -e "s/-//g" \
    data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.fixed.ind > \
    data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.simplepop.ind

sed -e "s/://g" data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.snp | \
    sed -e "s/^\( \+[^ ]\+ \+[^ ]\+ \+\)[^ ]\+\( \)/\10\2/" > \
    data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.fixed.snp

PAR=${PREFIX}.par.smartpca

echo "genotypename:    ${PREFIX}.eigenstratgeno" >  $PAR
echo "snpname:         ${PREFIX}.fixed.snp"            >> $PAR
echo "indivname:       ${PREFIX}.simplepop.ind"      >> $PAR

echo "evecoutname:     ${PREFIX}.evec"           >> $PAR
echo "evaloutname:     ${PREFIX}.eval"           >> $PAR

#echo "numthreads:      12"                       >> $PAR
#echo "fastmode:        YES"                      >> $PAR

# ----------------------------------------------------------------------------------------

echo "genotypename:    ${PREFIX}.eigenstratgeno" >  $PAR
echo "snpname:         ${PREFIX}.snp"            >> $PAR
echo "indivname:       ${PREFIX}.fixed.ind"      >> $PAR

echo "evecoutname:     ${PREFIX}.evec"           >> $PAR
echo "evaloutname:     ${PREFIX}.eval"           >> $PAR

echo "fstonly:         YES"                      >> $PAR
echo "fstz:            YES"                      >> $PAR

echo "phylipoutname:   ${PREFIX}.phylip"         >> $PAR
echo "numthreads:      12"                       >> $PAR

# ----------------------------------------------------------------------------------------

smartpca -p $PAR >> $PREFIX.smartpca.out 2>&1








genotypename:    ../CONVERTF/example.ped
snpname:         ../CONVERTF/example.map
indivname:       ../CONVERTF/example.ped
evecoutname:     example.evec
evaloutname:     example.eval
altnormstyle:    NO
numoutevec:      2
familynames:     NO
grmoutname:      grmjunk
# ----------------------------------------------------------------------------------------

smartpca -p $PAR >> $PREFIX.smartpca.out 2>&1





genotypename: data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.eigenstratgeno
5002192

head -n 10 data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.eigenstratgeno > data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.eigenstratgeno.mini

snpname: data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.snp
5002192

head -n 10 data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.snp > data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.snp.mini

indivname: data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.fixednodash.ind
881

inds  snp
   5    5 works
   6    6 fails
   7    7 fails
  10   10 fails
  50   50 fails

inds  snp
   5    5 works
   5    6 works
   5  600 works
   6    6 fails


NUM_INDS=5
NUM_SNPS=600

# One line per SNP, one colum per ind
sed -e "s/\([0-9]\)/\1\t/g" ../data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.eigenstratgeno | cut -f1-$NUM_INDS | sed -e "s/\t//g" | head -n $NUM_SNPS> tmp.eigenstratgeno
# One line per SNP
head -n $NUM_SNPS ../data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.snp > tmp.snp
head -n $NUM_INDS ../data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.strict.simplepop.ind > tmp.ind




