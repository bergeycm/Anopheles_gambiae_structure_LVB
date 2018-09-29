#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

library(ggplot2)

# ----------------------------------------------------------------------------------------
# --- Plot XP-EHH
# ----------------------------------------------------------------------------------------

args = commandArgs(trailingOnly = TRUE)
xpehh.in = args[1]  # "results/selscan/xp-ehh.BUGALA-BUKASA.2R.xpehh.out.norm"

regex = ".*xp-ehh\\.(.*)-(.*)\\.(.*)\\.xpehh.out.norm"
pop1 = gsub(regex, "\\1", xpehh.in)
pop2 = gsub(regex, "\\2", xpehh.in)
chr  = gsub(regex, "\\3", xpehh.in)

xpehh = read.table(xpehh.in,
    colClasses = c("NULL", "numeric", rep("NULL", 5), rep("numeric", 2), "NULL"),
    sep="\t", row.names=NULL, header=TRUE)

xpehh$chr = chr

samp = xpehh[sample(1:nrow(xpehh), size=100000),]

insecticide = read.table("data/insecticide_genes.bed")

this.insecticide = rowMeans(insecticide[insecticide$V1 == chr, 2:3])

p = ggplot(samp, aes(pos, normxpehh, col=abs(normxpehh))) +
    geom_vline(xintercept=this.insecticide, lty=2, col='grey', alpha=0.5) +
    geom_point(alpha=0.5) +
    theme_bw() +
    ylab(paste("Normalized XP-EHH\n", pop1, "vs.", pop2)) +
    xlab(paste0("Position along ", chr)) +
    scale_colour_gradient(low='white', high='red', guide=FALSE) +
    scale_x_continuous(labels=function(x) paste(x/1000000, "Mb"), expand=c(0,0)) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

dir.create("reports/selscan", showWarnings=FALSE)
out.pdf = paste0(gsub("results", "reports", xpehh.in), ".pdf")

ggsave(out.pdf, p, height=3, width=8)

# Also make tile plot that's smaller in size

island.sites   = c("BANDA", "BUGALA", "BUGALAIS",
                   "BUKASA", "KOOME", "NSADZI", "SSERINYA")
mainland.sites = c("BUGALAML", "BUWAMA", "ENTEBBE",
                   "KAZZI", "KIYINDI", "MITYANA")

island.col = "#0A318C"
mainland.col = "#D28C00"

if (pop1 %in% island.sites) {
    low.val = island.col
} else {
    low.val = mainland.col
}

if (pop2 %in% island.sites) {
    high.val = island.col
} else {
    high.val = mainland.col
}

q1 = ggplot(samp, aes(pos, normxpehh)) +
    geom_vline(xintercept=this.insecticide, lty=2, col='pink', alpha=0.5) +
    stat_bin2d(aes(fill = ..density..), bins=200) +
    theme_bw() +
    ylab(paste("Normalized XP-EHH\n", pop1, "vs.", pop2)) +
    xlab(paste0("Position along ", chr)) +
    scale_fill_gradient(guide=FALSE) +
    scale_x_continuous(labels=function(x) paste(x/1000000, "Mb"), expand=c(0,0)) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ylim(c(-6,6))

# Not shaded by density
q2 = ggplot(samp, aes(pos, normxpehh, z=normxpehh)) +
    geom_vline(xintercept=this.insecticide, lty=2, col='grey', alpha=0.5) +
    stat_summary_2d(bins=100) +
    theme_bw() +
    ylab(paste("Normalized XP-EHH\n", pop1, "vs.", pop2)) +
    xlab(paste0("Position along ", chr)) +
    scale_fill_gradientn(
        colors=c(rep(low.val, 2), 'white', rep(high.val, 2)),
        values=c(-100, -2, 0, 2, 100),
        rescaler = function(x, ...) x,
        oob = identity,
        guide=FALSE) +
    scale_x_continuous(labels=function(x) paste(x/1000000, "Mb"), expand=c(0,0)) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ylim(c(-6,6))

out.pdf.sm = paste0(gsub("results", "reports", xpehh.in), ".sm.pdf")

ggsave(out.pdf.sm, q1, height=3, width=8)
