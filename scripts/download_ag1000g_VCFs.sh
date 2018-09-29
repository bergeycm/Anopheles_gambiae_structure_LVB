#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Download Ag1000G VCFs
# ----------------------------------------------------------------------------------------

mkdir -p data/ag1000g.phase1.ar3/

FTP_PREFIX=ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3/variation/main/vcf/
OUT_PREFIX=data/ag1000g.phase1.ar3/ag1000g.phase1.ar3.pass.biallelic

for CHR in 2L 2R 3L 3R X; do

    wget $FTP_PREFIX/ag1000g.phase1.ar3.pass.biallelic.$CHR.vcf.gz \
        -O $OUT_PREFIX.$CHR.vcf.gz

    wget $FTP_PREFIX/ag1000g.phase1.ar3.pass.biallelic.$CHR.vcf.gz.tbi \
        -O $OUT_PREFIX.$CHR.vcf.gz.tbi

done

# ----------------------------------------------------------------------------------------
# --- Download Ag1000G’s "accessibility" metric files
# ----------------------------------------------------------------------------------------

ACC_DIR=data/accessibility
mkdir -p $ACC_DIR

for CHR in 2L 2R 3L 3R X; do

    ACC_URL=ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3/accessibility
    ACC_URL=$ACC_URL/accessibility.$CHR.vcf.gz

    wget $ACC_URL -O $ACC_DIR/accessibility.$CHR.vcf.gz

    gunzip -c $ACC_DIR/accessibility.$CHR.vcf.gz > $ACC_DIR/accessibility.$CHR.vcf

done
