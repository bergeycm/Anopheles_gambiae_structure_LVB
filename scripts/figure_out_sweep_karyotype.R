#!/usr/bin/Rscript

options(stringsAsFactors = FALSE)

library(ggplot2)
library(grid)
library(gridExtra)
library(lattice)
library(scales)

# ----------------------------------------------------------------------------------------
# --- Figure out karyotype of individuals with sweep on 2L
# ----------------------------------------------------------------------------------------

sweep.inds = read.table("reports/individuals_with_uganda_specific_2L_haplotype.txt",
    header=TRUE, sep="\t")

mdspca.in = "data/ssese_with_ag1000g/chr2L.pass.snp.phased.ag1000g.mds"

mdspca = read.table(mdspca.in, header = T)

mdspca$has.sweep = FALSE

mdspca[mdspca$IID %in% sweep.inds$ind,]$has.sweep = TRUE

p = ggplot(mdspca, aes(C1, C2, col=has.sweep)) + geom_point()

ggsave(p, file="tmp_mds.pdf")
