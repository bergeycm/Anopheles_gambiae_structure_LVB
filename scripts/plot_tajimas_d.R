#!/usr/bin/env Rscript

# ========================================================================================
# --- Plot Tajima's D
# ========================================================================================

library(ggplot2)
library(xtable)

options(stringsAsFactors=FALSE)

# Just 3 for now
tajd.out = list.files(path="results/", pattern="chr3.*\\.*.Tajima.D$", full.names=TRUE)

tajd.all = do.call(rbind, lapply(tajd.out, function(x) {
    this.tajd = read.table(x, header=TRUE, sep="\t")
    this.tajd$pop = x
    this.tajd
}))

tajd.all$pop = gsub(".*chr3[LR]\\.(.*)\\.Tajima.D", "\\1", tajd.all$pop)

tajd.all = tajd.all[tajd.all$pop != "BUGALA",]

tajd.all$island_mainland = "Island"
is.sites = c("BANDA", "BUGALAIS", "BUKASA", "NSADZI", "SSERINYA")
ml.sites = c("BUGALAML", "BUWAMA", "KAZZI", "KIYINDI", "MITYANA")
tajd.all[tajd.all$pop %in% ml.sites,]$island_mainland = "Mainland"

tajd.all[tajd.all$pop == "BUGALAIS",]$pop = "BUGALA (I)"
tajd.all[tajd.all$pop == "BUGALAML",]$pop = "BUGALA (M)"
is.sites[is.sites == "BUGALAIS"] = "BUGALA (I)"
ml.sites[ml.sites == "BUGALAML"] = "BUGALA (M)"

tajd.all[tajd.all$pop == "KAZZI",]$pop = "KAAZI"
tajd.all[tajd.all$pop == "MITYANA",]$pop = "WAMALA"
ml.sites[ml.sites == "KAZZI"] = "KAAZI"
ml.sites[ml.sites == "MITYANA"] = "WAMALA"

simple.cap = function(x) {
    s = strsplit(x, " ")[[1]]
    lstr = paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
        sep="", collapse=" ")

    gsub("\\(i\\)", "(I)", gsub("\\(m\\)", "(M)", lstr))
}

tajd.all$pop = as.vector(sapply(tajd.all$pop, simple.cap))
is.sites = as.vector(sapply(is.sites, simple.cap))
ml.sites = as.vector(sapply(ml.sites, simple.cap))

p = ggplot(tajd.all, aes(factor(pop, levels=rev(c(is.sites, ml.sites))), TajimaD,
        fill=factor(island_mainland, levels=c("Mainland", "Island")))) +
    geom_boxplot(outlier.size=0.5, outlier.shape = NA) +
    xlab("") +
    ylab(expression(atop("Tajima's"~D))) +
    scale_fill_discrete(guide=FALSE) +
    coord_flip() +
    theme_bw(base_size=15) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black"))

ggsave(p, filename="reports/tajimas_d_boxplot.pdf", height=3, width=4)

# Also plot island vs. mainland
p = ggplot(tajd.all, aes(island_mainland, TajimaD,
        color=factor(island_mainland, levels=rev(c("Mainland", "Island"))))) +
    geom_boxplot() +
    theme_bw() +
    xlab("") +
    ylab(expression("Tajima's"~D)) +
    scale_color_discrete(guide=FALSE) +
    coord_flip()

ggsave(p, filename="reports/tajimas_d_boxplot.mainland-island.pdf", height=3, width=4)
