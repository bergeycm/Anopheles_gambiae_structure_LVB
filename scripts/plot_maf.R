#!/usr/bin/env Rscript

# ========================================================================================
# --- Plot MAF
# ========================================================================================

library(ggplot2)
library(xtable)

options(stringsAsFactors=FALSE)

# Just 3 for now
maf.out = list.files(path="results/", pattern="chr3.*\\.*.sample.frq$", full.names=TRUE)

maf.all = do.call(rbind, lapply(maf.out, function(x) {
    this.maf = read.table(x, header=FALSE, skip=1, sep="\t")
    this.maf$pop = x
    this.maf
}))

names(maf.all) = c("CHROM", "POS", "N_ALLELES", "N_CHR", "A1.freq", "A2.freq", "pop")

maf.all$A1.freq = gsub(".*:", "", maf.all$A1.freq)
maf.all$A2.freq = gsub(".*:", "", maf.all$A2.freq)

maf.all$MAF = as.numeric(do.call(pmin, maf.all[,5:6]))

maf.all$pop = gsub(".*chr3[LR]\\.(.*)\\.sample.frq", "\\1", maf.all$pop)

# Remove unsplit Bugala
maf.all = maf.all[maf.all$pop != "BUGALA",]

maf.all$island_mainland = "Island"
is.sites = c("BANDA", "BUGALAIS", "BUKASA", "NSADZI", "SSERINYA")
ml.sites = c("BUGALAML", "BUWAMA", "KAZZI", "KIYINDI", "MITYANA")
maf.all[maf.all$pop %in% ml.sites,]$island_mainland = "Mainland"

maf.all[maf.all$pop == "BUGALAIS",]$pop = "BUGALA (I)"
maf.all[maf.all$pop == "BUGALAML",]$pop = "BUGALA (M)"
is.sites[is.sites == "BUGALAIS"] = "BUGALA (I)"
ml.sites[ml.sites == "BUGALAML"] = "BUGALA (M)"

p = ggplot(maf.all, aes(factor(pop, levels=rev(c(is.sites, ml.sites))), MAF,
        color=factor(island_mainland, levels=c("Mainland", "Island")))) +
    theme_bw() +
    xlab("") +
    ylab("Minor allele frequency") +
    scale_color_discrete(guide=FALSE) +
    coord_flip() +
    ylim(c(0,0.25))

ggsave(p + geom_boxplot(outlier.size=0.5), filename="reports/MAF_boxplot.pdf",
    height=3, width=4)
ggsave(p + geom_violin(outlier.size=0.5), filename="reports/MAF_violin.pdf",
    height=3, width=4)

simple.cap = function(x) {
    s = strsplit(x, " ")[[1]]
    lstr = paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
        sep="", collapse=" ")

    gsub("\\(i\\)", "(I)", gsub("\\(m\\)", "(M)", lstr))
}

maf.all$pop = as.vector(sapply(maf.all$pop, simple.cap))
is.sites = as.vector(sapply(is.sites, simple.cap))
ml.sites = as.vector(sapply(ml.sites, simple.cap))

# Plot MAF histogram
p = ggplot(maf.all, aes(MAF, ..ncount.., group = pop,
      col=factor(island_mainland, levels=c("Mainland", "Island")))) +
    geom_freqpoly(binwidth = 0.05) +
    xlim(c(0,0.5)) +
    guides(color=FALSE) +
    ylab("SNP density") +
    xlab("Minor allele frequency") +
    theme_bw(base_size=15) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black"))

ggsave(p, filename="reports/MAF_histogram_SFS.pdf",
    height=3, width=4)

# Also plot island vs. mainland
p = ggplot(maf.all, aes(island_mainland, MAF,
        color=factor(island_mainland, levels=rev(c("Mainland", "Island"))))) +
    geom_boxplot() +
    theme_bw(base_size=15) +
    xlab("") +
    ylab("Minor allele frequency") +
    scale_color_discrete(guide=FALSE) +
    coord_flip()

ggsave(p, filename="reports/MAF_boxplot.mainland-island.pdf", height=3, width=4)
