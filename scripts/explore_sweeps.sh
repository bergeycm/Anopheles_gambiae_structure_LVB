#!/bin/bash

# ========================================================================================
# --- Explore around rgd and awh
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Annotate
# ----------------------------------------------------------------------------------------

for CHR in X 3L 3R 2L 2R; do
    sh scripts/annotate_SNPs_with_snpEff.sh data/ssese_with_ag1000g/ssese_with_ag1000g.$CHR.flt.strict.vcf
    sh scripts/annotate_SNPs_with_snpEff.sh data/chr$CHR.pass.snp.vcf
done

# ========================================================================================
# --- Pull out region around gene
# ========================================================================================

module load tabix

# --- rdg

echo -e "X\t9215505\t9266532" > data/interval.rdg.bed
echo -e "X\t8215505\t10266532" > data/interval.rdg.1mb.bed

# Strict, just gene
tabix -h data/ssese_with_ag1000g/ssese_with_ag1000g.X.flt.strict.snp.eff.vcf.gz \
    data/interval.rdg.bed -B > \
    data/ssese_with_ag1000g/ssese_with_ag1000g.X.flt.strict.rdg.vcf

# Strict, gene +/- 1Mb
tabix -h data/ssese_with_ag1000g/ssese_with_ag1000g.X.flt.strict.snp.eff.vcf.gz \
    data/interval.rdg.1mb.bed -B > \
    data/ssese_with_ag1000g/ssese_with_ag1000g.X.flt.strict.rdg.1mb.vcf

# Upstream, just gene
tabix -h data/chrX.pass.snp.snp.eff.vcf.gz \
    data/interval.rdg.bed -B > \
    data/chrX.pass.snp.rdg.vcf

# Upstream, gene +/- 1Mb
tabix -h data/chrX.pass.snp.snp.eff.vcf.gz \
    data/interval.rdg.1mb.bed -B > \
    data/chrX.pass.snp.rdg.1mb.vcf

# --- awh

echo -e "2L\t34022548\t34024344" > data/interval.awh.bed
echo -e "2L\t33022548\t35024344" > data/interval.awh.1mb.bed

# Strict, just gene
tabix -h data/ssese_with_ag1000g/ssese_with_ag1000g.2L.flt.strict.snp.eff.vcf.gz \
    data/interval.awh.bed -B > \
    data/ssese_with_ag1000g/ssese_with_ag1000g.2L.flt.strict.awh.vcf    # TO RUN    

# Strict, gene +/- 1Mb
tabix -h data/ssese_with_ag1000g/ssese_with_ag1000g.2L.flt.strict.snp.eff.vcf.gz \
    data/interval.awh.1mb.bed -B > \
    data/ssese_with_ag1000g/ssese_with_ag1000g.2L.flt.strict.awh.1mb.vcf   # TO RUN    

# Upstream, just gene
tabix -h data/chr2L.pass.snp.snp.eff.vcf.gz \
    data/interval.awh.bed -B > \
    data/chr2L.pass.snp.awh.vcf

# Upstream, gene +/- 1Mb
tabix -h data/chr2L.pass.snp.snp.eff.vcf.gz \
    data/interval.awh.1mb.bed -B > \
    data/chr2L.pass.snp.awh.1mb.vcf

# ========================================================================================
# --- Subset haps files to grab region +/- 1Mb around ROIs
# ========================================================================================

# X    8215505    10266532    rdg +/- 1 Mb

awk '{ if ($3 >= 8215505 && $3 <= 10266532) print $0}' \
    data/chrX.pass.snp.phased.haps > \
    data/chrX.pass.snp.phased.rdg.1mb.haps
cp data/chrX.pass.snp.phased.sample data/chrX.pass.snp.phased.rdg.1mb.sample

# 2L    33022548    35024344    awh +/- 1 Mb

awk '{ if ($3 >= 33022548 && $3 <= 35024344) print $0}' \
    data/chr2L.pass.snp.phased.haps > \
    data/chr2L.pass.snp.phased.awh.1mb.haps
cp data/chr2L.pass.snp.phased.sample data/chr2L.pass.snp.phased.awh.1mb.sample

# ========================================================================================
# --- Copy relevant data files into directory for transport
# ========================================================================================

CP_DIR=data/ssese_putative_sweep_data
mkdir $CP_DIR

cp data/chrX.pass.snp.rdg.1mb.vcf            $CP_DIR
cp data/chr2L.pass.snp.awh.1mb.vcf           $CP_DIR

cp data/chrX.pass.snp.phased.rdg.1mb.haps    $CP_DIR
cp data/chrX.pass.snp.phased.rdg.1mb.sample  $CP_DIR
cp data/chr2L.pass.snp.phased.awh.1mb.haps   $CP_DIR
cp data/chr2L.pass.snp.phased.awh.1mb.sample $CP_DIR

cp data/ssese_individual_info_simple.txt     $CP_DIR

tar -zcvf $CP_DIR.tar.gz $CP_DIR
