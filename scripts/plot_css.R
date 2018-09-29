#!/usr/bin/env Rscript

# module load R/3.3.0

# ----------------------------------------------------------------------------------------
# --- Plot CSS
# ----------------------------------------------------------------------------------------

options(stringsAsFactors=FALSE)

# Change settings to avoid use of exponential notation (e.g. 1.81e+08)
options("scipen"=100, "digits"=4)

library(ggplot2)

args = commandArgs(trailingOnly = TRUE)

chr = args[1]

css = read.table("results/css.all.txt", header=FALSE)

names(css) = c("chr", "pos", "pop1", "pop2", "css")

css = css[css$chr == chr,]

# Remove duplicate comparisons (BANDA-BUGALA stays, but BUGALA-BANDA goes)
css = css[css$pop1 < css$pop2,]

# ----------------------------------------------------------------------------------------
# --- Bring in location info for inversions, heterochromatic regions, and genes
# ----------------------------------------------------------------------------------------

# --- Genes

genes = read.table("data/insecticide_genes.bed")
names(genes) = c("chr", "start", "end", "gene")

this.genes = genes[genes$chr == chr,]

this.genes$midpoint = rowMeans(data.frame(this.genes$start, this.genes$end))

# --- Inversions

inversion_sites = read.table("data/inversion_simple.bed")
names(inversion_sites) = c("chr", "start", "end", "name")

this.inversion_sites = inversion_sites[inversion_sites$chr == chr,]

# --- Heterochromatic

het = read.table("data/heterochromatin.bed")
names(het) = c("chr", "start", "end")

this.het = het[het$chr == chr,]

# ----------------------------------------------------------------------------------------
# --- Do plot of all stats together
# ----------------------------------------------------------------------------------------

y.min = 0
y.max = 6

p = ggplot(css, aes(pos, css),
        environment = environment())

if (nrow(this.het) > 0) {
    p = p + geom_rect(data=this.het,
            aes(x=NULL, y=NULL, col=NULL, xmin=start, xmax=end, ymin=y.min, ymax=y.max),
            alpha=0.9, fill='grey')
}

if (nrow(this.inversion_sites) > 0) {
    p = p + geom_rect(data=this.inversion_sites,
            aes(x=NULL, y=NULL, col=NULL, xmin=start, xmax=end, ymin=y.min, ymax=y.max),
            alpha=0.2, fill='lightblue')
}

if (length(this.genes) > 0) {
    p = p + geom_vline(xintercept = c(this.genes$midpoint), lty=3)
}

p = p + geom_point(aes(color = css)) +
    scale_color_gradientn(colors=c("orange", "blue")) +
    facet_grid(pop1 ~ pop2) +
    coord_cartesian(ylim=c(y.min, y.max)) +
    ylab("CSS") + xlab("") +
    scale_x_continuous(labels=function(x) paste(x/1000000, "Mb"), expand=c(0,0)) +
    theme_bw() +
    theme(panel.margin = unit(2, "pt"),
        #axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.text.y = element_text(size = 8, angle = 270),
        strip.background = element_rect(fill = 'white', color='white'))

ggsave(p, file=paste0("reports/css_", chr, ".all.pdf"),
    width=20, height=20)
