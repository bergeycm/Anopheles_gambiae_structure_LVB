#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Plot Tajima's D by site
# ----------------------------------------------------------------------------------------

library(ggplot2)

chrs = c("2R", "2L", "3R", "3L", "X")

tajd = do.call(rbind, lapply(chrs, function (chr) {

    this.chr.tajd.files = list.files(path="results/",
        pattern=paste0("chr", chr, "\\...*\\.Tajima.D"),
        full.names=TRUE)

    this.chr.tajd = do.call(rbind, lapply(this.chr.tajd.files, function (tajd.file) {
        x = read.table(tajd.file, sep="\t", header=TRUE)
        x$pop = gsub(".*\\.([^\\.]+)\\.Tajima.D", "\\1", tajd.file)
        x
    }))
    this.chr.tajd$chr = chr
    this.chr.tajd
}))

tajd.sm = unique(rbind(tajd[seq(from=1, to=nrow(tajd), by=10),],
                       tajd[tajd$TajimaD >= quantile(tajd$TajimaD, 0.99, na.rm=TRUE),],
                       tajd[tajd$TajimaD <= quantile(tajd$TajimaD, 0.01, na.rm=TRUE),]))

# Change Bugala labels
tajd.sm$pop[tajd.sm$pop == "BUGALAIS"] = "BUGALA (I)"
tajd.sm$pop[tajd.sm$pop == "BUGALAML"] = "BUGALA (M)"

# Remove unsplit Bugala
tajd.sm = tajd.sm[tajd.sm$pop != "BUGALA",]

tajd.sm = na.omit(tajd.sm[tajd.sm$TajimaD != "NaN",])

is.sites = c("BANDA", "BUGALA (I)", "BUKASA", "NSADZI", "SSERINYA")
ml.sites = c("BUGALA (M)", "BUWAMA", "KAZZI", "KIYINDI", "MITYANA")

tajd.sm[tajd.sm$pop == "KAZZI",]$pop = "KAAZI"
tajd.sm[tajd.sm$pop == "MITYANA",]$pop = "WAMALA"
ml.sites[ml.sites == "KAZZI"] = "KAAZI"
ml.sites[ml.sites == "MITYANA"] = "WAMALA"

tajd.sm$pop = factor(tajd.sm$pop, levels = c("", is.sites, ml.sites))

tajd.sm$chr = factor(tajd.sm$chr, levels = chrs)

# ----------------------------------------------------------------------------------------
# --- Bring in location info for inversions, heterochromatic regions, and genes
# ----------------------------------------------------------------------------------------

# --- Genes

genes = read.table("data/insecticide_genes_small.bed")
names(genes) = c("chr", "midpoint", "gene")

genes = do.call(rbind, lapply(c(is.sites, ml.sites),
    function (pop) { cbind(genes, pop) }))

genes$pop = factor(genes$pop, levels = c("", is.sites, ml.sites))
genes$chr = factor(genes$chr, levels = chrs)

# --- Inversions

inversion_sites = read.table("data/inversion_simple.bed")
names(inversion_sites) = c("chr", "start", "end", "name")

inversion_sites = do.call(rbind, lapply(c(is.sites, ml.sites),
    function (pop) { cbind(inversion_sites, pop) }))

inversion_sites$pop = factor(inversion_sites$pop, levels = c("", is.sites, ml.sites))
inversion_sites$chr = factor(inversion_sites$chr, levels = chrs)

# --- Heterochromatic

het = read.table("data/heterochromatin.bed")
names(het) = c("chr", "start", "end")

het = do.call(rbind, lapply(c(is.sites, ml.sites),
    function (pop) { cbind(het, pop) }))

het$pop = factor(het$pop, levels = c("", is.sites, ml.sites))
het$chr = factor(het$chr, levels = chrs)

# --- New sweeps

new.sweeps = data.frame(chr=c("X", "2L"), loc=c(9238942, 34044820))

new.sweeps = do.call(rbind, lapply(c(is.sites, ml.sites),
    function (pop) { cbind(new.sweeps, pop) }))

new.sweeps$pop = factor(new.sweeps$pop, levels = c("", is.sites, ml.sites))
new.sweeps$chr = factor(new.sweeps$chr, levels = chrs)

# ----------------------------------------------------------------------------------------
# --- Gene labels
# ----------------------------------------------------------------------------------------

gene.label.text = rbind(
    data.frame(BIN_START = 28497407, lab = "Cyp6p",
        chr = factor("2R", levels = levels(tajd.sm$chr))),
    data.frame(BIN_START =  2394888, lab = "Vgsc",
        chr = factor("2L", levels = levels(tajd.sm$chr))),
    data.frame(BIN_START = 25399104, lab = "Gaba",
        chr = factor("2L", levels = levels(tajd.sm$chr))),
    data.frame(BIN_START = 28598038, lab = "Gste",
        chr = factor("3R", levels = levels(tajd.sm$chr))),
    data.frame(BIN_START = 11204486, lab = "Tep1",
        chr = factor("3L", levels = levels(tajd.sm$chr))),
    data.frame(BIN_START = 15241718, lab = "Cyp9k1",
        chr = factor("X",  levels = levels(tajd.sm$chr)))
)

label.height = 50
label.expand.factor = 1

gene.label.text$TajimaD = label.height
gene.label.text$pop = factor("", levels = levels(tajd.sm$pop))

gene.label.fake.pts = gene.label.text
gene.label.fake.pts$BIN_START = 0
gene.label.fake.pts$TajimaD = seq(from = label.height - label.expand.factor,
                              to   = label.height + label.expand.factor,
                              length.out = 6)

novel.label.text = data.frame(BIN_START=c(9238942, 34044820), lab="",
    chr=factor(c("X", "2L"), levels = levels(tajd.sm$chr)),
    TajimaD = label.height, pop = factor("", levels = levels(tajd.sm$pop)))

# For putting a line at the origin
origin.df = data.frame(origin=0,
    pop = factor(c(is.sites, ml.sites), levels = levels(tajd.sm$pop)))

# For drawing in fake y-axis
axis.df = data.frame(x=0, y=-3, xend=0, yend=3,
    pop = factor(c(is.sites, ml.sites), levels = levels(tajd.sm$pop)),
    chr = factor("2R", levels = levels(tajd.sm$chr)))

# ----------------------------------------------------------------------------------------
# --- Do plot
# ----------------------------------------------------------------------------------------

y.min = -3; y.max = 4

p = ggplot(tajd.sm, aes(BIN_START, TajimaD, color=pop %in% is.sites)) +
    geom_hline(data=origin.df, aes(yintercept = origin), color = 'grey') +
    geom_rect(data=inversion_sites,
        aes(x=NULL, y=NULL, col=NULL,
            xmin=start, xmax=end, ymin=y.min, ymax=y.max + 1.4),
        alpha=0.2, fill='lightgrey') +
    geom_vline(data=genes, aes(xintercept = midpoint), lty=3, alpha=0.5) +
    geom_vline(data=new.sweeps, aes(xintercept = loc),
        lty=1, alpha=1, col="black", lwd=0.5) +
    geom_line() +
    geom_rect(data=het,
        aes(x=NULL, y=NULL, col=NULL,
            xmin=start, xmax=end, ymin=y.min, ymax=y.max + 1.4),
        alpha=1, fill='grey') +
    geom_segment(data=axis.df, aes(x=x, y=y, xend=xend, yend=yend), col="black", lwd=1) +
    facet_grid(pop ~ chr, scales = "free", space = "free", shrink=FALSE) +
    xlab("Genome position (Mb)") +
    ylab("Tajima's D") +
    guides(color=FALSE) +
    scale_x_continuous(labels = function(x) { x/1000000 },
                       breaks = seq(from=0e6, to=60e6, by=10e6),
                       expand = c(0,0)) +
    scale_y_continuous(breaks = c(-3,0,3), expand=c(0,0), limits=c(NA, NA)) +
    theme_bw() +
    theme(strip.text = element_text(size=7, lineheight=5.0),
          strip.background = element_rect(fill="white", color="black", size=0.5),
          panel.spacing.y = unit(0, "lines"),
          axis.line.x = element_line(color = "black"),
          axis.line.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = "grey", linetype = "dashed"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    geom_point(data = gene.label.fake.pts, mapping = aes(x = BIN_START, y = TajimaD),
               pch=NA, col=NA, fill=NA) +
    geom_point(data = gene.label.text, mapping = aes(x = BIN_START, y = TajimaD),
               pch=25, col='black', fill='black') +
    geom_text(data = gene.label.text,
              mapping = aes(x = BIN_START + 1000000, y = TajimaD, label = lab),
              color="black", fontface="italic", size=2, hjust = 'left') +
    geom_point(data = novel.label.text, mapping = aes(x = BIN_START, y = TajimaD),
               pch=25, col='black', fill='white')

ggsave("reports/tajimas_D_by_site.pdf", p, height=8, width=12)
