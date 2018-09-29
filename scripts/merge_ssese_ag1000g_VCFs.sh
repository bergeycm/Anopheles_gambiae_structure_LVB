#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Merge (and filter) Ssese and Ag1000G VCFs
# ----------------------------------------------------------------------------------------

CHR=$1
SUFFIX=$2    # Empty string or ".missingtoref"

module load vcftools
module load bcftools/1.5
module load tabix/1.5

AG_VCF=data/ag1000g.phase1.ar3/ag1000g.phase1.ar3.pass.biallelic.$CHR.vcf.gz

SS_VCF=data/chr$CHR.pass.snp.vcf.gz

echo "--- Merging files --- `date`"

if [[ "$SUFFIX" == ".missingtoref" ]]; then

    mkdir -p data/ssese_with_ag1000g.missingtoref

    # Raw merge file without filtration and assuming missing sites are reference
    # -0 --missing-to-ref    assume genotypes at missing sites are 0/0

    bcftools merge -0 -O z $SS_VCF $AG_VCF > \
        data/ssese_with_ag1000g$SUFFIX/ssese_with_ag1000g.$CHR$SUFFIX.vcf.gz

else

    # Raw merged file without filtration
    bcftools merge -O z $SS_VCF $AG_VCF > \
        data/ssese_with_ag1000g/ssese_with_ag1000g.$CHR.vcf.gz
fi

echo "--- Finished merging files --- `date`"

# ----------------------------------------------------------------------------------------

# --- Finding high missingness sites to filter:

echo "--- Finding high missingness sites to filter --- `date`"

# Temporarily unzip
gunzip -c data/ssese_with_ag1000g$SUFFIX/ssese_with_ag1000g.$CHR$SUFFIX.vcf.gz > \
    data/ssese_with_ag1000g$SUFFIX/ssese_with_ag1000g.$CHR$SUFFIX.vcf

# Copy in header
head -n 1000 data/ssese_with_ag1000g$SUFFIX/ssese_with_ag1000g.$CHR$SUFFIX.vcf | \
    grep "^#"  > \
    data/ssese_with_ag1000g$SUFFIX/ssese_with_ag1000g.$CHR$SUFFIX.flt.vcf

# Filter VCF
perl scripts/filter_merged_vcfs.pl \
    data/ssese_with_ag1000g$SUFFIX/ssese_with_ag1000g.$CHR$SUFFIX.vcf \
    data/accessibility/accessibility.$CHR.vcf >> \
    data/ssese_with_ag1000g$SUFFIX/ssese_with_ag1000g.$CHR$SUFFIX.flt.vcf

echo "--- Finished finding high missingness sites to filter --- `date`"

# ----------------------------------------------------------------------------------------

# Just make link if we're in missingtoref mode (e.g. don't filter on missingness)

if [[ "$SUFFIX" == ".missingtoref" ]]; then

    ln -s \
        `pwd`/data/ssese_with_ag1000g$SUFFIX/ssese_with_ag1000g.$CHR$SUFFIX.flt.vcf \
        `pwd`/data/ssese_with_ag1000g$SUFFIX/ssese_with_ag1000g.$CHR$SUFFIX.flt.strict.vcf

else

    echo "--- Filtering merged files for missingness --- `date`"

    # --- Use missingness file created in above step to make file filtered for missingness

    MISSINGNESS=data/ssese_with_ag1000g/ssese_with_ag1000g.missingness.$CHR$SUFFIX.txt

    # STRICT filtering: anything with >10% missing in either dataset is cut.
    # Resultant list is positions to KEEP not to exclude.
    awk '{ if ($3 < 0.1 && $4 < 0.1) print $0 }' $MISSINGNESS > \
        ${MISSINGNESS/.txt/.strict.txt}

    vcftools --vcf data/ssese_with_ag1000g/ssese_with_ag1000g.$CHR$SUFFIX.flt.vcf \
        --positions ${MISSINGNESS/.txt/.strict.txt} \
        --recode \
        --out data/ssese_with_ag1000g/ssese_with_ag1000g.$CHR$SUFFIX.flt.strict.vcf

    # Rename to remove VCFtools's ".recode.vcf"
    mv data/ssese_with_ag1000g/ssese_with_ag1000g.$CHR$SUFFIX.flt.strict.vcf.recode.vcf \
        data/ssese_with_ag1000g/ssese_with_ag1000g.$CHR$SUFFIX.flt.strict.vcf

    # Clean-up temporarily unzipped file
    rm data/ssese_with_ag1000g/ssese_with_ag1000g.$CHR$SUFFIX.vcf

    echo "--- Finished filtering merged files for missingness --- `date`"

fi

# ----------------------------------------------------------------------------------------

# Compress and index filtered VCF

echo "--- Compressing and indexing filtered VCF --- `date`"

bgzip -c data/ssese_with_ag1000g$SUFFIX/ssese_with_ag1000g.$CHR$SUFFIX.flt.strict.vcf > \
         data/ssese_with_ag1000g$SUFFIX/ssese_with_ag1000g.$CHR$SUFFIX.flt.strict.vcf.gz
tabix -p vcf \
    data/ssese_with_ag1000g$SUFFIX/ssese_with_ag1000g.$CHR$SUFFIX.flt.strict.vcf.gz

echo "--- Finished compressing and indexing filtered VCF --- `date`"

exit
