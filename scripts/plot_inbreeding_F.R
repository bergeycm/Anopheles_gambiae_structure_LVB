#!/usr/bin/env Rscript

# ========================================================================================
# --- Plot inbreeding statistic, F
# ========================================================================================

library(ggplot2)
library(xtable)

options(stringsAsFactors=FALSE)

inb.out = list.files(path="results/", pattern="all.*\\.*.het$", full.names=TRUE)

inb.all = do.call(rbind, lapply(inb.out, function(x) {
    this.inb = read.table(x, header=TRUE, sep="\t")
    this.inb$pop = x
    this.inb
}))

inb.all$pop = gsub(".*all\\.(.*)\\.het", "\\1", inb.all$pop)

inb.all = inb.all[inb.all$pop != "BUGALA",]

inb.all$island_mainland = "Island"
is.sites = c("BANDA", "BUGALAIS", "BUKASA", "NSADZI", "SSERINYA")
ml.sites = c("BUGALAML", "BUWAMA", "KAZZI", "KIYINDI", "MITYANA")
inb.all[inb.all$pop %in% ml.sites,]$island_mainland = "Mainland"

inb.all[inb.all$pop == "BUGALAIS",]$pop = "BUGALA (I)"
inb.all[inb.all$pop == "BUGALAML",]$pop = "BUGALA (M)"
is.sites[is.sites == "BUGALAIS"] = "BUGALA (I)"
ml.sites[ml.sites == "BUGALAML"] = "BUGALA (M)"

inb.all[inb.all$pop == "KAZZI",]$pop   = "KAAZI"
inb.all[inb.all$pop == "MITYANA",]$pop = "WAMALA"
ml.sites[ml.sites == "KAZZI"] = "KAAZI"
ml.sites[ml.sites == "MITYANA"] = "WAMALA"

simple.cap = function(x) {
    s = strsplit(x, " ")[[1]]
    lstr = paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
        sep="", collapse=" ")

    gsub("\\(i\\)", "(I)", gsub("\\(m\\)", "(M)", lstr))
}

inb.all$pop = as.vector(sapply(inb.all$pop, simple.cap))
is.sites = as.vector(sapply(is.sites, simple.cap))
ml.sites = as.vector(sapply(ml.sites, simple.cap))

p = ggplot(inb.all, aes(factor(pop, levels=rev(c(is.sites, ml.sites))), F,
        fill=factor(island_mainland, levels=c("Mainland", "Island")))) +
    geom_boxplot(outlier.size=0.5, outlier.shape = NA) +
    xlab("") +
    ylab(expression(atop("Inbreeding statistic "~F))) +
    scale_fill_discrete(guide=FALSE) +
    coord_flip() +
    theme_bw(base_size=15) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black"))

ggsave(p, filename="reports/inbreeding_boxplot.pdf", height=3, width=4)

# Also plot island vs. mainland
p = ggplot(inb.all, aes(island_mainland, F,
        color=factor(island_mainland, levels=rev(c("Mainland", "Island"))))) +
    geom_boxplot() +
    theme_bw() +
    xlab("") +
    ylab(expression("Inbreeding statistic"~F)) +
    scale_color_discrete(guide=FALSE) +
    coord_flip()

ggsave(p, filename="reports/inbreeding_boxplot.mainland-island.pdf", height=3, width=4)
