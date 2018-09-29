#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

library(ggplot2)

# ----------------------------------------------------------------------------------------
# --- Make ADMIXTURE CV plot
# ----------------------------------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
in.file = args[1]  # E.g. "data/ssese_with_ag1000g/adm_subsets/chr3.withM.CV.txt"

cv = read.table(in.file)

p = ggplot(cv, aes(V1, V2, col=V1==6)) +
    geom_point() +
    ylab("CV Error") +
    xlab("K") +
    scale_color_manual(values=c("black", "red"),
                       guide=FALSE) +
    theme_bw()

ggsave(p, file=gsub(".txt", ".pdf", in.file), width=4, height=4)
