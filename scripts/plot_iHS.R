#!/usr/bin/env Rscript

# module load R/3.3.0

# ----------------------------------------------------------------------------------------
# --- Plot iHS
# ----------------------------------------------------------------------------------------

options(stringsAsFactors=FALSE)

# Change settings to avoid use of exponential notation (e.g. 1.81e+08)
options("scipen"=100, "digits"=4)

library(ggplot2)

args = commandArgs(trailingOnly = TRUE)

chr = args[1]

write(paste0("Processing chr", chr, "."), stderr())

ihs.files = list.files(path="results/selscan",
        pattern="*.ihs.out.100bins.norm",
        full.names=TRUE)

this.chr.ihs.files = ihs.files[grepl(chr, ihs.files)]

read.stat = function (stat.file) {

    this.stat = read.table(stat.file, header=FALSE)

    names(this.stat) = c("locusID", "physicalPos", "freq1", "ihh1", "ihh0",
                         "ihs.raw", "ihs.norm", "extreme.flag")

    this.stat$pop = strsplit(stat.file, "\\.")[[1]][2]

    return(this.stat)
}

all.ihs = do.call(rbind, lapply(this.chr.ihs.files, read.stat))

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
y.max = 7.25

# Take absolute value of iHS and reduce dataset to just extreme values
all.ihs$ihs.norm.abs = abs(all.ihs$ihs.norm)
all.ihs.sm = all.ihs[all.ihs$ihs.norm.abs >= 2,]

p = ggplot(all.ihs.sm, aes(physicalPos, ihs.norm.abs),
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

p = p + geom_point(aes(color = ihs.norm.abs)) +
    scale_color_gradientn(colors=c("grey", "black")) +
    facet_grid(pop ~ .) +
    coord_cartesian(ylim=c(y.min, y.max)) +
    ylab("|iHS|") + xlab("") +
    scale_x_continuous(labels=function(x) paste(x/1000000, "Mb"), expand=c(0,0)) +
    theme_bw() +
    theme(panel.margin = unit(2, "pt"),
        #axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.text.y = element_text(size = 8, angle = 270),
        strip.background = element_rect(fill = 'white', color='white'))

ggsave(p, file=paste0("reports/ihs_", chr, ".all.pdf"),
    width=8, height=9)
