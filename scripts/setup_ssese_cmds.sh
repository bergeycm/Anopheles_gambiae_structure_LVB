#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Set up project directory for Ssese WGS analysis
# ----------------------------------------------------------------------------------------

cd /scratch365/cbergey/ssese_wgs

git clone https://github.com/bergeycm/ssese_wgs.git

cd ssese_wgs

# ----------------------------------------------------------------------------------------
# --- Set variables for Snakemake
# ----------------------------------------------------------------------------------------

sh set_vars_for_snakemake.sh

# ----------------------------------------------------------------------------------------
# --- Compress and tabix index VCF files
# ----------------------------------------------------------------------------------------

export PATH=$PATH:$HOME/bin:$HOME/bin/vcftools_0.1.12b/bin:$HOME/bin/tabix-0.2.6

NGS_DIR=../NGS-map/AgamP4_AG_snps

for CHR in 2L 2R 3L 3R X; do
    bgzip -c $NGS_DIR/chr$CHR.pass.snp.vcf > $NGS_DIR/chr$CHR.pass.snp.vcf.gz
    tabix -p vcf $NGS_DIR/chr$CHR.pass.snp.vcf.gz
done

# ----------------------------------------------------------------------------------------
# --- Link to SNPs
# ----------------------------------------------------------------------------------------

NGS_FULL_PATH=/scratch365/cbergey/ssese_wgs/NGS-map/AgamP4_AG_snps

for CHR in 2L 2R 3L 3R X; do
    ln -s $NGS_FULL_PATH/chr$CHR.pass.snp.vcf        data/chr$CHR.pass.snp.vcf
    ln -s $NGS_FULL_PATH/chr$CHR.pass.snp.vcf.idx    data/chr$CHR.pass.snp.vcf.idx
    ln -s $NGS_FULL_PATH/chr$CHR.pass.snp.vcf.gz     data/chr$CHR.pass.snp.vcf.gz
    ln -s $NGS_FULL_PATH/chr$CHR.pass.snp.vcf.gz.tbi data/chr$CHR.pass.snp.vcf.gz.tbi
done

# ----------------------------------------------------------------------------------------
# --- Copy in other data
# ----------------------------------------------------------------------------------------

cp /afs/crc.nd.edu/group/BesanskyNGS/ssese_wgs_bergey/ssese_individual_info.csv data/
cp /afs/crc.nd.edu/group/BesanskyNGS/ssese_wgs_bergey/ssese_plate* data/

exit
