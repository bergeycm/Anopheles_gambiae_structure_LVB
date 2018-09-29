#!/bin/bash

# ========================================================================================
# --- Make input files for dadi
# ========================================================================================

# Call with:
# sh scripts/make_dadi_input.sh
# It will run the conversion script for the full genome (sans inversions, etc.) in modes
# species, island, and mainlandisland,
# where species means do all for this species, island means split by island/site,
# and mainlandisland means split by whether animal comes from an island or the mainland.

# module load perl
module load vcftools
module load plink/1.9
module load parallel

IND_INFO=data/ssese_individual_info_simple_bugala_split.txt

# ----------------------------------------------------------------------------------------
# --- Create VCF with inversions and heterochromatic regions removed
# ----------------------------------------------------------------------------------------

# Create list of input BED files (by chr) - EXCLUDING the X-chromosome
TMP_BED_LIST=`mktemp tmp.chr_bed_list.XXXXX`
ls data/chr*.pass.snp.flt.bed | grep -v "X" | \
    sed -e "s/\(.*\).bed/\1.bed \1.bim \1.fam/" > $TMP_BED_LIST

# Add heterochromatic regions to inversion BED and switch to chromosome numbers
TMP_INVHET_SIMP=`mktemp tmp.inversion_het_simple.XXXXX`
cat data/inversion_simple.bed data/heterochromatin.bed | grep -v "X" | \
    sed -e "s/\([0-9]\)$/\1\thet/" | \
    sed -e "s/^2L/1/" -e "s/^2R/2/" -e "s/^3L/3/" -e "s/^3R/4/" > $TMP_INVHET_SIMP

# Make SNP set covering the 2La and 2Rb inversions and heterochromatic regions
plink --merge-list $TMP_BED_LIST \
    --make-set $TMP_INVHET_SIMP --write-set \
    --out data/all.pass.snp.phased.ag1000g.inversionshet

# Make list of all SNPs to exclude: now just inv/het
TMP_EXCLUDE=`mktemp tmp.exclude.XXXXX`
cat data/all.pass.snp.phased.ag1000g.inversionshet.set > $TMP_EXCLUDE

# Make temporary VCF sans inversions and het regions
TMP_VCF_PRE=`mktemp tmp.all.noinvhet.XXXXX`

plink --merge-list $TMP_BED_LIST --exclude $TMP_EXCLUDE \
    --recode vcf-iid --out $TMP_VCF_PRE

# Fix chromosome numbering [move from PLINK format (e.g. 1) to standard format (e.g. 2L)]
sed -i'.bak' -e "s/^1\t/2L\t/" -e "s/^2\t/2R\t/" \
             -e "s/^3\t/3L\t/" -e "s/^4\t/3R\t/" $TMP_VCF_PRE.vcf
rm $TMP_VCF_PRE.vcf.bak

# ----------------------------------------------------------------------------------------
# --- Create pop files that list individuals and their groups for dadi analyses
# ----------------------------------------------------------------------------------------

mkdir -p data/dadi

# --- Species mode

# This filename determines the output dadi filename
POP_FILE=data/dadi/dadi_pop_list.species.txt
cut -d' ' -f2,6 $IND_INFO | grep "AG" > $POP_FILE

# --- Island/site mode

POP_FILE=data/dadi/dadi_pop_list.island.txt
awk '{ if ($6 == "AG") print $2,$3 }' $IND_INFO > $POP_FILE

# --- Island vs. mainland mode

POP_FILE=data/dadi/dadi_pop_list.mainlandisland.txt
awk '{ if ($6 == "AG") print $2,$3 }' $IND_INFO | \
    sed -e "s/BUWAMA/_MAINLAND/"  -e "s/KAZZI/_MAINLAND/" \
        -e "s/KIYINDI/_MAINLAND/" -e "s/MITYANA/_MAINLAND/" | \
    sed -e "s/ [^_].*/ ISLAND/" | sed -e "s/_//" > $POP_FILE

# --- Small island vs. mainland mode

POP_FILE=data/dadi/dadi_pop_list.mainlandsmallisland.txt
awk '{ if ($6 == "AG") print $2,$3 }' $IND_INFO | \
    sed -e "s/BUWAMA/_MAINLAND/"  -e "s/KAZZI/_MAINLAND/" \
        -e "s/KIYINDI/_MAINLAND/" -e "s/MITYANA/_MAINLAND/" | \
    sed -e "s/BUGALA/_BUGALA/" | \
    sed -e "s/ [^_].*/ ISLAND/" | sed -e "s/_//" > $POP_FILE

# ----------------------------------------------------------------------------------------
# --- Do conversion of VCF to dadi format for all analyses to be run
# ----------------------------------------------------------------------------------------

# --- Temporarily download Anopheles genome FASTA file

AGAM_URL=https://www.vectorbase.org/sites/default/files/ftp/downloads/
AGAM_URL=${AGAM_URL}Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa.gz

AGAM_FA=data/dadi/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa
wget $AGAM_URL -O $AGAM_FA.gz
gunzip -f ${AGAM_FA}.gz

# --- Do conversion

mkdir -p parallel_tmp

#for mode in species island mainlandisland mainlandsmallisland; do
for mode in island mainlandsmallisland; do

    POP_FILE=data/dadi/dadi_pop_list.${mode}.txt
    VCF=$TMP_VCF_PRE.vcf

    NUM_LINES=`wc -l $VCF | cut -d' ' -f1`

    echo date > tmp_$mode.dadi.convert.commands.txt

    for START in `seq 0 1000000 $NUM_LINES`; do

        END=$((START + 999999))
        echo "Processing lines $START to $END..."

        echo "perl scripts/convert_vcf_to_dadi_input.pl \
            $AGAM_FA $VCF $POP_FILE \
            $START $END &> $POP_FILE.$START.$END.log.txt" \
            >> tmp_$mode.dadi.convert.commands.txt
    done

    echo date >> tmp_$mode.dadi.convert.commands.txt

    parallel -a tmp_$mode.dadi.convert.commands.txt --jobs 8 --tmpdir parallel_tmp

    # Combine output
    head -n 1 `ls -v data/dadi/dadi_pop_list.$mode.txt.*.data | head -n1` > $POP_FILE.data
    cat `ls -v data/dadi/dadi_pop_list.$mode.txt.*.data` | \
        grep -v "^NAME" >> $POP_FILE.data

done

# --- Rename output file

#for mode in species island mainlandisland mainlandsmallisland; do
for mode in island mainlandsmallisland; do

    POP_FILE=data/dadi/dadi_pop_list.${mode}.txt
    OUT_DADI=${POP_FILE%.txt}.data
    OUT_DADI=${OUT_DADI/_pop_list/}
    mv $POP_FILE.data $OUT_DADI

done

# ----------------------------------------------------------------------------------------
# --- Make smaller file of chr3 only
# ----------------------------------------------------------------------------------------

#for mode in species island mainlandisland mainlandsmallisland; do
for mode in island mainlandsmallisland; do

    FULL=data/dadi/dadi.$mode.data
    CHR3=${FULL/dadi./dadi.chr3.}

    # Chromosome is second from the last column.
    head -n1 $FULL > $CHR3
    awk '{ if ($(NF-1) == "3R" || $(NF-1) == "3L") print $0 }' $FULL >> $CHR3

done

rm $AGAM_FA
rm $TMP_VCF_PRE
rm $TMP_VCF_PRE.{log,vcf}
rm $TMP_VCF_PRE.{bak,bed,bim,fam,nosex}
rm $TMP_BED_LIST $TMP_INVHET_SIMP $TMP_EXCLUDE

exit
