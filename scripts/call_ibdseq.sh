#!/usr/bin/sh

# ----------------------------------------------------------------------------------------
# --- Infer IBD regions with IBDseq
# ----------------------------------------------------------------------------------------

CHR=$1

module load vcftools
module load ibdseq

IN_VCF=data/chr$CHR.pass.snp.flt.vcf.gz

OUT_PREFIX=results/`basename $IN_VCF | sed -e "s/vcf.gz/ibdseq/"`

# Temporarily decompress VCF
# gunzip -c $IN_VCF > ${IN_VCF/.vcf.gz/.forIBDseq.vcf}
vcftools --gzvcf $IN_VCF \
    --maf 0.05 \
    --recode --out ${IN_VCF/.vcf.gz/.forIBDseq}

if [[ "$CHR" == "3L" ]]; then
    REGION='3L:15000000-41000000'
fi
if [[ "$CHR" == "3R" ]]; then
    REGION='3R:1-37000000'
fi

java -Xmx24g -jar $IBDSEQ_JAR \
    gt=${IN_VCF/.vcf.gz/.forIBDseq.recode.vcf} \
    chrom=$REGION \
    out=$OUT_PREFIX \
    nthreads=20

rm ${IN_VCF/.vcf.gz/.forIBDseq.recode.vcf}
rm ${IN_VCF/.vcf.gz/.forIBDseq.log}

# Obsolete. Now Ag1000G map is used

#    # --- Fake map file
#    # Since we lack map for Anopheles, assume 1Mb = 1 cM
#
#    awk '{ printf "%s %s %.8f %s\n", $1,$2,$3/1000000,$3 }' \
#        data/chr$CHR.pass.snp.phased.haps > \
#        data/chr$CHR.pass.snp.phased.fake.map

exit
