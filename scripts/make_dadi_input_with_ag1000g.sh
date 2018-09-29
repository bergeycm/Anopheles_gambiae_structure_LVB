#!/bin/bash

# ========================================================================================
# --- Make input files for dadi from Ag1000G data
# ========================================================================================

# Call with:
# sh scripts/make_dadi_input_with_ag1000g.sh
# It will run the conversion script for the full genome (sans inversions, etc.)
# to create dataset with all Ag1000G Uganda animals.

module load parallel/20150122

# ----------------------------------------------------------------------------------------

export VCF=data/ag1000g.phase1.ar3/ag1000g.phase1.ar3.pass.biallelic.all.for.dadi.vcf

export POP_FILE=data/dadi/dadi_pop_list.ag1000g.txt

export HUNK_SIZE=100000

# --- Temporarily download Anopheles genome FASTA file
AGAM_URL=https://www.vectorbase.org/sites/default/files/ftp/downloads/
AGAM_URL=${AGAM_URL}Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa.gz
export AGAM_FA=data/dadi/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa
wget $AGAM_URL -O $AGAM_FA.gz
rm -f ${AGAM_FA}
gunzip ${AGAM_FA}.gz

# --- Figure out how many lines in VCF (overestimate SNP count)
VCF_LINES=`wc -l $VCF | cut -d' ' -f1`

# --- Do conversion
convert_to_dadi() {
    start=$1
    end=$((start + $HUNK_SIZE))
    echo "Processing $start to $end."
    perl scripts/convert_vcf_to_dadi_input.pl \
        $AGAM_FA $VCF $POP_FILE $start $end &> $POP_FILE.log.${start}-${end}.txt
}

export -f convert_to_dadi

seq 0 $HUNK_SIZE $((VCF_LINES + $HUNK_SIZE)) | parallel --jobs 20 convert_to_dadi

# --- Combine output

head -n1 data/dadi/dadi_pop_list.ag1000g.txt.0*.data \
    > data/dadi/dadi.ag1000g.data

cat `ls -v data/dadi/dadi_pop_list.ag1000g.txt.*.data` | \
    grep -v "^NAME" >> data/dadi/dadi.ag1000g.data
