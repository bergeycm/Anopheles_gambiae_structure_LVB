#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Plot H12 by site
# ----------------------------------------------------------------------------------------

library(ggplot2)

chrs = c("2L", "2R", "3L", "3R", "X")

pi = do.call(rbind, lapply(chrs, function (chr) {

    this.chr.pi.files = list.files(path="results/",
        pattern=paste0("chr", chr, "\\...*\\.windowed.pi"),
        full.names=TRUE)

    this.chr.pi = do.call(rbind, lapply(this.chr.pi.files, function (pi.file) {
        x = read.table(pi.file, sep="\t", header=TRUE)
        x$pop = gsub(".*\\.([^\\.]+)\\.windowed.pi", "\\1", pi.file)
        x
    }))
    this.chr.pi$chr = chr
    this.chr.pi
}))

pi.sm = unique(rbind(pi[seq(from=1, to=nrow(pi), by=10),],
                     pi[pi$PI >= quantile(pi$PI, 0.99, na.rm=TRUE),],
                     pi[pi$PI <= quantile(pi$PI, 0.01, na.rm=TRUE),]))

# Change Bugala labels
pi.sm$pop[pi.sm$pop == "BUGALAIS"] = "BUGALA (I)"
pi.sm$pop[pi.sm$pop == "BUGALAML"] = "BUGALA (M)"

# Remove unsplit Bugala
pi.sm = pi.sm[pi.sm$pop != "BUGALA",]

pi.sm = na.omit(pi.sm[pi.sm$PI != "NaN",])

is.sites = c("BANDA", "BUGALA (I)", "BUKASA", "NSADZI", "SSERINYA")
ml.sites = c("BUGALA (M)", "BUWAMA", "KAZZI", "KIYINDI", "MITYANA")
pi.sm$pop = factor(pi.sm$pop, levels = c(is.sites, ml.sites))

# ----------------------------------------------------------------------------------------
# --- Bring in location info for inversions, heterochromatic regions, and genes
# ----------------------------------------------------------------------------------------

# --- Genes

genes = read.table("data/insecticide_genes.bed")
names(genes) = c("chr", "start", "end", "gene")

genes$midpoint = rowMeans(data.frame(genes$start, genes$end))

# --- Inversions

inversion_sites = read.table("data/inversion_simple.bed")
names(inversion_sites) = c("chr", "start", "end", "name")

# --- Heterochromatic

het = read.table("data/heterochromatin.bed")
names(het) = c("chr", "start", "end")

# --- New sweeps

new.sweeps = data.frame(chr=c("X", "2L"), loc=c(9238942, 34044820))

# ----------------------------------------------------------------------------------------
# --- Do plot
# ----------------------------------------------------------------------------------------

y.min = 0; y.max = 0.03

p = ggplot(pi.sm, aes(BIN_START, PI, color=pop %in% is.sites)) +
    geom_rect(data=het,
        aes(x=NULL, y=NULL, col=NULL,
            xmin=start, xmax=end, ymin=y.min, ymax=y.max),
        alpha=0.9, fill='grey') +
    geom_rect(data=inversion_sites,
        aes(x=NULL, y=NULL, col=NULL,
            xmin=start, xmax=end, ymin=y.min, ymax=y.max),
        alpha=0.9, fill='grey') +
    geom_vline(data=genes, aes(xintercept = midpoint), lty=3, alpha=0.5) +
    geom_vline(data=new.sweeps, aes(xintercept = loc), lty=2, alpha=1, col="black") +
    geom_line() +
    facet_grid(pop ~ chr, scales = "free_x", space = "free_x") +
    xlab("Window center") +
    ylab(expression(pi)) +
    guides(color=FALSE) +
    scale_x_continuous(labels = function(x) paste(x/1000000, "Mb"),
                       breaks = seq(from=20e6, to=60e6, by=20e6),
                       expand = c(0,0)) +
    scale_y_continuous(breaks = c(y.min,y.max), expand=c(0,0.0001), limits=c(y.min,y.max)) +
    theme_bw() +
    theme(strip.text = element_text(size=7, lineheight=5.0),
          strip.background = element_rect(fill="white", color="black",
          size=0.5))

ggsave("reports/pi_by_site.pdf", p, height=8, width=12)
