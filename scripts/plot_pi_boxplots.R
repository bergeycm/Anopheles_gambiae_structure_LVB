#!/usr/bin/env Rscript

# ========================================================================================
# --- Plot nucleotide diversity, pi
# ========================================================================================

library(ggplot2)
library(xtable)

options(stringsAsFactors=FALSE)

# Just 3 for now
pi.out = list.files(path="results/", pattern="chr3.*\\.*.windowed.pi$", full.names=TRUE)

pi.all = do.call(rbind, lapply(pi.out, function(x) {
    this.pi = read.table(x, header=TRUE, sep="\t")
    this.pi$pop = x
    this.pi
}))

pi.all$pop = gsub(".*chr3[LR]\\.(.*)\\.windowed.pi", "\\1", pi.all$pop)

# Remove unsplit Bugala
pi.all = pi.all[pi.all$pop != "BUGALA",]

pi.all$island_mainland = "Island"
is.sites = c("BANDA", "BUGALAIS", "BUKASA", "NSADZI", "SSERINYA")
ml.sites = c("BUGALAML", "BUWAMA", "KAZZI", "KIYINDI", "MITYANA")
pi.all[pi.all$pop %in% ml.sites,]$island_mainland = "Mainland"

pi.all[pi.all$pop == "BUGALAIS",]$pop = "BUGALA (I)"
pi.all[pi.all$pop == "BUGALAML",]$pop = "BUGALA (M)"
is.sites[is.sites == "BUGALAIS"] = "BUGALA (I)"
ml.sites[ml.sites == "BUGALAML"] = "BUGALA (M)"

pi.all[pi.all$pop == "KAZZI",]$pop = "KAAZI"
pi.all[pi.all$pop == "MITYANA",]$pop = "WAMALA"
ml.sites[ml.sites == "KAZZI"] = "KAAZI"
ml.sites[ml.sites == "MITYANA"] = "WAMALA"

simple.cap = function(x) {
    s = strsplit(x, " ")[[1]]
    lstr = paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
        sep="", collapse=" ")

    gsub("\\(i\\)", "(I)", gsub("\\(m\\)", "(M)", lstr))
}

pi.all$pop = as.vector(sapply(pi.all$pop, simple.cap))
is.sites = as.vector(sapply(is.sites, simple.cap))
ml.sites = as.vector(sapply(ml.sites, simple.cap))

p = ggplot(pi.all, aes(factor(pop, levels=rev(c(is.sites, ml.sites))), PI,
        #color=factor(island_mainland, levels=c("Mainland", "Island")),
        fill=factor(island_mainland, levels=c("Mainland", "Island")))) +
    geom_boxplot(outlier.size=0.5, outlier.shape = NA) +
    xlab("") +
    ylab(expression(atop("Nucleotide diversity"~pi))) +
    scale_color_discrete(guide=FALSE) +
    scale_fill_discrete(guide=FALSE) +
    coord_flip() +
    theme_bw(base_size=15) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black"))

ggsave(p, filename="reports/pi_boxplot.pdf", height=3, width=4)

# Also plot island vs. mainland
p = ggplot(pi.all, aes(island_mainland, PI,
        color=factor(island_mainland, levels=rev(c("Mainland", "Island"))))) +
    geom_boxplot() +
    theme_bw() +
    xlab("") +
    ylab(expression("Nucleotide diversity"~pi)) +
    scale_color_discrete(guide=FALSE) +
    coord_flip()

ggsave(p, filename="reports/pi_boxplot.mainland-island.pdf", height=3, width=4)
