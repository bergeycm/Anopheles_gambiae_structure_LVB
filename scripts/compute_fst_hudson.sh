#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Compute Fst with EIGENSTRAT smartpca (Uses Hudson's estimator and computes S.E.)
# ----------------------------------------------------------------------------------------

module load eigensoft

PREFIX=data/all.pass.snp.flt.eigen
PAR=${PREFIX}.par.smartpca

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

smartpca -p $PAR >> data/all.pass.snp.flt.eigen.fst.se.out 2>&1
