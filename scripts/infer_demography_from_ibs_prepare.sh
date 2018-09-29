#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Prepare to infer demography from IBS
# ----------------------------------------------------------------------------------------

module load vcftools
module load parallel

IN_VCF=data/chr3L.pass.snp.phased.haps.vcf

SITES=$(ls data/ssese.seqids.is-*.txt | cut -d'-' -f 2 | \
    grep -v "BUGALA\.txt" | sed -e "s/\.txt//")

make_popdata () {

    site=$1

    SITE_INDS=`cat data/ssese.seqids.is-$site.txt | \
        grep -v "ssese99\." | tr "\n" "," | sed -e "s/,$//"`

    vcf-query $IN_VCF -f '%CHROM\t%POS\t%REF\t[%GT]\n' -c $SITE_INDS | \
        sed -e "s:[/\|]::g" | sed -e "s/^3L/3/" | \
        awk '{ if ($2 > 2000000 && $2 < 3000000) print $0 }' \
        > results/IBS.$site.popdata
}

export -f make_popdata
export IN_VCF=$IN_VCF

echo ${SITES[*]} | tr " " "\n" | parallel make_popdata

# ----------------------------------------------------------------------------------------

# --- Set symbolic links to scripts

ln -s ~/work/bin/Inferring-demography-from-IBS/condense_sorted_lengths.py .
ln -s ~/work/bin/Inferring-demography-from-IBS/infer_onepop_cumulative.py .
ln -s ~/work/bin/Inferring-demography-from-IBS/infer_onepop_cumulative_specify_start.py .

# ----------------------------------------------------------------------------------------

# Set the mutation rate to 5.5e-9
#     Since theta = 4Ne*mu, assuming Ne is a million, theta = 0.0055
# Set rho to 13.4 / kb (Genome-Wide Fine-Scale Recombination Rate Variation in
#                         Drosophila melanogaster, Chan et al 2012)
#     or 0.0134 / bp.

VARMU=~/work/bin/Inferring-demography-from-IBS/calc_ibs_backcoal_varmu.py
cp ${VARMU}{,_ORIGINAL}

sed -e "s/^theta=.*$/theta=0.0055/" ${VARMU}_ORIGINAL |
    sed -e "s/^rho=.*$/rho=0.0134/" > $VARMU
