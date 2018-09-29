#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Set up project directory for Ssese WGS analysis - PSU version
# ----------------------------------------------------------------------------------------

cd /gpfs/scratch/cxb585/ssese-wgs

git clone https://bergeycm@github.com/bergeycm/ssese_wgs.git

cd ssese_wgs

# ----------------------------------------------------------------------------------------
# --- Set variables for Snakemake
# ----------------------------------------------------------------------------------------

sh set_vars_for_snakemake_psu.sh

# ----------------------------------------------------------------------------------------
# --- Bring in VCF files and indexes
# ----------------------------------------------------------------------------------------

scp cbergey@crcfe01.crc.nd.edu:/scratch365/cbergey/ssese_wgs/NGS-map/AgamP4_AG_snps/*.pass.*gz* data/

# ----------------------------------------------------------------------------------------
# --- Unzip VCF files and index
# ----------------------------------------------------------------------------------------

module load tabix/0.2.6

for CHR in 2L 2R 3L 3R X; do
    gunzip -c data/chr$CHR.pass.snp.vcf.gz > data/chr$CHR.pass.snp.vcf
done

# ----------------------------------------------------------------------------------------
# --- Copy in other data
# ----------------------------------------------------------------------------------------

scp cbergey@crcfe01.crc.nd.edu:/afs/crc.nd.edu/group/BesanskyNGS/ssese_wgs_bergey/ssese_individual_info.csv data/
scp cbergey@crcfe01.crc.nd.edu:/afs/crc.nd.edu/group/BesanskyNGS/ssese_wgs_bergey/ssese_plate* data/
scp cbergey@crcfe01.crc.nd.edu:/afs/crc.nd.edu/group/BesanskyNGS/ssese_wgs_bergey/heterochromatin.bed data/
scp cbergey@crcfe01.crc.nd.edu:/afs/crc.nd.edu/group/BesanskyNGS/ssese_wgs_bergey/inversion_sites.bed data/
scp cbergey@crcfe01.crc.nd.edu:/afs/crc.nd.edu/group/BesanskyNGS/ssese_wgs_bergey/insecticide_genes.bed data/

# ----------------------------------------------------------------------------------------
# --- Copy in phased data from Ag1000G
# ----------------------------------------------------------------------------------------

AG_DIR=data/ag1000g.phase1.ar3/
mkdir -p $AG_DIR

FTP=ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/haplotypes/main/shapeit/

for CHR in 2L 2R 3L 3R X; do
    wget $FTP/ag1000g.phase1.ar3.1.haplotypes.$CHR.haps.gz -O $AG_DIR/ag1000g.phase1.ar3.1.haplotypes.$CHR.haps.gz
    gunzip $AG_DIR/ag1000g.phase1.ar3.1.haplotypes.$CHR.haps.gz
    wget $FTP/ag1000g.phase1.ar3.1.haplotypes.$CHR.sample.gz -O $AG_DIR/ag1000g.phase1.ar3.1.haplotypes.$CHR.sample.gz
    gunzip $AG_DIR/ag1000g.phase1.ar3.1.haplotypes.$CHR.sample.gz
done

wget ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/samples/samples.all.txt \
 -O data/ag1000g.phase1.ar3/samples.all.txt
