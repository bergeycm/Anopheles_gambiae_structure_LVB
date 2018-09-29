#!/bin/sh

module load java

SNPEFF=${HOME}/bin/snpEff/

# Get Ensembl name of PEST genome
AGAM_ENS=`java -jar ${SNPEFF}/snpEff.jar databases | grep "Anopheles_gambiae" | tail -n 1 | cut -f1`

GENOME_NAME=Anopheles_gambiae

# ----------------------------------------------------------------------------------------

# Grab chromosome via command line
CHR=$1

VCF=data/chr${CHR}.pass.snp.vcf
OUT_HTML=${VCF/.snp.vcf/.snp.eff.snpEff_summary.html}
OUT_HTML=${OUT_HTML/data/reports}

echo "Running snpEff on chr $CHR";
java -Xmx4g -jar ${SNPEFF}/snpEff.jar \
    -v -v -v $AGAM_ENS \
    -dataDir $GRP/snpEff_DB \
    -stats $OUT_HTML \
    $VCF \
    > ${VCF/.snp.vcf/.snp.eff.vcf}

# Ultimately this should be piped to gzip.

exit
