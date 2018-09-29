#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Download phased Ag1000G data (shapeit)
# ----------------------------------------------------------------------------------------

mkdir -p data/ag1000g.phase1.ar3/

# Get sample info
wget ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3/samples/samples.all.txt \
    -O data/ag1000g.phase1.ar3/samples.all.txt

FTP_PREFIX=ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3.1/haplotypes/main/shapeit
OUT_PREFIX=data/ag1000g.phase1.ar3/ag1000g.phase1.ar3.1.haplotypes

for CHR in 2L 2R 3L 3R X; do

    for FILE in sample haps; do

        wget $FTP_PREFIX/ag1000g.phase1.ar3.1.haplotypes.$CHR.$FILE.gz \
            -O $OUT_PREFIX.$CHR.$FILE.gz

        gunzip -c $OUT_PREFIX.$CHR.$FILE.gz > $OUT_PREFIX.$CHR.$FILE
    done
done

exit
