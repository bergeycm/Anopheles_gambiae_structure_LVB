#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Plot H12 by site
# ----------------------------------------------------------------------------------------

library(ggplot2)

chrs = c("2R", "2L", "3R", "3L", "X")

h12 = do.call(rbind, lapply(chrs, function (chr) {

    this.chr.h12.files = list.files(path="results/",
        pattern=paste0("chr", chr, "\\.pass\\.snp\\.phased\\..*\\.H12\\.out\\.txt"),
        full.names=TRUE)

    this.chr.h12 = do.call(rbind, lapply(this.chr.h12.files, function (h12.file) {
        x = read.table(h12.file, sep="\t")
        x$pop = gsub(".*phased\\.([^\\.]+)\\..*", "\\1", h12.file)
        x
    }))
    this.chr.h12$chr = chr
    this.chr.h12
}))

names(h12) = c("win.center", "win.left", "win.right", "k.haplotypes", "hap.freq.spec",
               "ind", "H1.het", "H2", "H12", "H2_H1", "pop", "chr")

h12.sm = unique(rbind(h12[seq(from=1, to=nrow(h12), by=100),],
                      h12[h12$H12 >= quantile(h12$H12, 0.99),]))

# Change Bugala labels
h12.sm$pop[h12.sm$pop == "BUGALAIS"] = "BUGALA (I)"
h12.sm$pop[h12.sm$pop == "BUGALAML"] = "BUGALA (M)"

# Fix names
h12.sm$pop[h12.sm$pop == "KAZZI"] = "KAAZI"
h12.sm$pop[h12.sm$pop == "MITYANA"] = "WAMALA"

is.sites = c("BANDA", "BUGALA (I)", "BUKASA", "NSADZI", "SSERINYA")
ml.sites = c("BUGALA (M)", "BUWAMA", "KAAZI", "KIYINDI", "WAMALA")
h12.sm$pop = factor(h12.sm$pop, levels = c("", is.sites, ml.sites))

h12.sm$chr = factor(h12.sm$chr, levels = chrs)

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
    data.frame(win.center = 28497407, lab = "Cyp6p",
        chr = factor("2R", levels = levels(h12.sm$chr))),
    data.frame(win.center =  2394888, lab = "Vgsc",
        chr = factor("2L", levels = levels(h12.sm$chr))),
    data.frame(win.center = 25399104, lab = "Gaba",
        chr = factor("2L", levels = levels(h12.sm$chr))),
    data.frame(win.center = 28598038, lab = "Gste",
        chr = factor("3R", levels = levels(h12.sm$chr))),
    data.frame(win.center = 11204486, lab = "Tep1",
        chr = factor("3L", levels = levels(h12.sm$chr))),
    data.frame(win.center = 15241718, lab = "Cyp9k1",
        chr = factor("X",  levels = levels(h12.sm$chr)))
)

label.height = 5
label.expand.factor = 0.1

gene.label.text$H12 = label.height
gene.label.text$pop = factor("", levels = levels(h12.sm$pop))

gene.label.fake.pts = gene.label.text
gene.label.fake.pts$win.center = 0
gene.label.fake.pts$H12 = seq(from = label.height - label.expand.factor,
                              to   = label.height + label.expand.factor,
                              length.out = 6)

novel.label.text = data.frame(win.center=c(9238942, 34044820), lab="",
    chr=factor(c("X", "2L"), levels = levels(h12.sm$chr)),
    H12=label.height, pop = factor("", levels = levels(h12.sm$pop)))

# For drawing in fake y-axis
axis.df = data.frame(x=0, y=0, xend=0, yend=1,
    pop = factor(c(is.sites, ml.sites), levels = levels(h12.sm$pop)),
    chr = factor("2R", levels = levels(h12.sm$chr)))

# ----------------------------------------------------------------------------------------
# --- Do plot
# ----------------------------------------------------------------------------------------

y.min = 0; y.max = 1

p = ggplot(h12.sm, aes(win.center, H12, color=pop %in% is.sites)) +
    geom_rect(data=inversion_sites,
        aes(x=NULL, y=NULL, col=NULL,
            xmin=start, xmax=end, ymin=y.min, ymax=y.max + 0.2),
        alpha=0.2, fill='lightgrey') +
    geom_vline(data=genes, aes(xintercept = midpoint), lty=3, alpha=0.5) +
    geom_vline(data=new.sweeps, aes(xintercept = loc),
        lty=1, alpha=1, col="black", lwd=0.25) +
    geom_line() +
    geom_rect(data=het,
        aes(x=NULL, y=NULL, col=NULL,
            xmin=start, xmax=end, ymin=y.min, ymax=y.max + 0.2, pop=pop),
        alpha=1, fill='grey') +
    geom_segment(data=axis.df, aes(x=x, y=y, xend=xend, yend=yend), col="black", lwd=1) +
    facet_grid(pop ~ chr, scales = "free", space = "free", shrink=FALSE) +
    xlab("Genome position (Mb)") +
    guides(color=FALSE) +
    #scale_x_continuous(labels = function(x) paste(x/1000000, "Mb"),
    #                   breaks = seq(from=10e6, to=60e6, by=10e6),
    #                   expand = c(0,0)) +
    scale_x_continuous(labels = function(x) { x/1000000 },
                       breaks = seq(from=0e6, to=60e6, by=10e6),
                       expand = c(0,0)) +
    scale_y_continuous(breaks = c(0,1), expand=c(0,0), limits=c(NA, NA)) +
    theme_bw() +
    theme(strip.text = element_text(size=7, lineheight=5.0),
          strip.background = element_rect(fill="white", color="black", size=0.5),
          panel.spacing.y = unit(0, "lines"),
          axis.line.x = element_line(color = "black"),
          axis.line.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    geom_point(data = gene.label.fake.pts, mapping = aes(x = win.center, y = H12),
               pch=NA, col=NA, fill=NA) +
    geom_point(data = gene.label.text, mapping = aes(x = win.center, y = H12),
               pch=25, col='black', fill='black') +
    geom_text(data = gene.label.text,
              mapping = aes(x = win.center + 1000000, y = H12, label = lab),
              color="black", fontface="italic", size=2, hjust = 'left') +
    geom_point(data = novel.label.text, mapping = aes(x = win.center, y = H12),
               pch=25, col='black', fill='white')

ggsave("reports/H12_by_site.pdf", p, height=8, width=12)
