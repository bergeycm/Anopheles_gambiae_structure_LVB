#!/bin/bash

# ========================================================================================
# --- Prepare to make input files for dadi from Ag1000G data
# ========================================================================================

module load perl
module load vcftools

# ----------------------------------------------------------------------------------------
# --- Create VCF with heterochromatic regions removed and LD-pruned
# ----------------------------------------------------------------------------------------

# --- Combine VCF files
export VCF=data/ag1000g.phase1.ar3/ag1000g.phase1.ar3.pass.biallelic.all.for.dadi.vcf

# Just chr3
vcf-concat data/ag1000g.phase1.ar3/ag1000g.phase1.ar3.pass.biallelic.3*.vcf.gz |
    sed -e "s/^3L\t/3\t/" -e "s/^3R\t/4\t/" > $VCF

PREFIX=${VCF/.vcf/}

# Remove heterochromatic regions and output PLINK format
vcftools --vcf $VCF --exclude-bed data/heterochromatin.bed --out $PREFIX \
    --plink --keep-INFO-all

# Fix chromosome numbering [move from PLINK format (e.g. 1) to standard format (e.g. 2L)]
sed -i'.bak' -e "s/^1\t/2L\t/" -e "s/^2\t/2R\t/" \
             -e "s/^3\t/3L\t/" -e "s/^4\t/3R\t/" $PREFIX.vcf
rm $PREFIX.vcf.bak

# ----------------------------------------------------------------------------------------
# --- Create pop files that list individuals and their groups for dadi analyses
# ----------------------------------------------------------------------------------------

export POP_FILE=data/dadi/dadi_pop_list.ag1000g.txt

grep "Uganda" data/ag1000g.phase1.ar3/samples.all.txt | \
    cut -f2 | \
    sed -e "s/$/ AG1000G_UG/" > $POP_FILE

# ----------------------------------------------------------------------------------------
# --- Do convesion of VCF to dadi format for all analyses to be run
# ----------------------------------------------------------------------------------------

AGAM_FA=data/dadi/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa

# --- Temporarily download Anopheles genome FASTA file

if [[ ! -e $AGAM_FA ]]; then

    AGAM_URL=https://www.vectorbase.org/sites/default/files/ftp/downloads/
    AGAM_URL=${AGAM_URL}Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa.gz

    wget $AGAM_URL -O $AGAM_FA.gz
    gunzip ${AGAM_FA}.gz
fi

# --- Do conversion

perl scripts/convert_vcf_to_dadi_input.pl \
    $AGAM_FA $VCF $POP_FILE &> $POP_FILE.log.txt

# --- Rename output file

OUT_DADI=${POP_FILE%.txt}.data
OUT_DADI=${OUT_DADI/_pop_list/}
mv $POP_FILE.data $OUT_DADI

exit
