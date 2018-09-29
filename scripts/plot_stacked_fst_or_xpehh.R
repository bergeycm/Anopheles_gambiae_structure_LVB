#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Make "stacked" plot of pairwise statistic
# ----------------------------------------------------------------------------------------

options(stringsAsFactors=FALSE)

# Change settings to avoid use of exponential notation (e.g. 1.81e+08)
options("scipen"=100, "digits"=4)

library(ggplot2)

chrs = c("2R", "2L", "3R", "3L", "X")

# ----------------------------------------------------------------------------------------

stat.files = list()

stat.files[["fst"]] = list.files(path="results/",
        pattern="*.windowed.*fst",
        full.names=TRUE)

stat.files[["xpehh"]] = list.files(path="results/selscan/",
    pattern="xp-ehh.*.xpehh.out.norm.avg",
    full.names=TRUE)

stat.files[["css"]] = "results/css.all.txt"

# Remove long distance Fst
stat.files[["fst"]] =
    stat.files[["fst"]][!grepl("long_distance", stat.files[["fst"]])]
# Remove unsplit Bugala
stat.files[["fst"]] =
    stat.files[["fst"]][!grepl("BUGALA[^MI]", stat.files[["fst"]])]
stat.files[["xpehh"]] =
    stat.files[["xpehh"]][!grepl("BUGALA[^MI]", stat.files[["xpehh"]])]

read.stat = function (stat.file) {

    if (grepl("fst", stat.file)) {
        pop.str = strsplit(stat.file, "\\.")[[1]][2]
        pops = strsplit(pop.str, "_")[[1]]
        chr = gsub("chr", "", strsplit(stat.file, "[\\./]")[[1]][3])
    } else if (grepl("xpehh", stat.file)) {
        pop.str = gsub(".*xp-ehh\\.(.*)\\..*\\.xpehh.*", "\\1", stat.file)
        pops = strsplit(pop.str, "-")[[1]]
        chr = gsub(".*xp-ehh\\..*\\.(.*)\\.xpehh.*", "\\1", stat.file)
    }

    if (grepl("css", stat.file)) {
        this.stat = read.table(stat.file, header=FALSE)
        names(this.stat) = c("chr", "pos", "pop1", "pop2", "statistic")

        # Double the data.frame, switching pop1 and pop2
        this.stat = rbind(this.stat,
                              data.frame(this.stat[c("chr", "pos", "statistic")],
                                         pop1 = this.stat$pop2,
                                         pop2 = this.stat$pop1))
    } else {
        this.stat = read.table(stat.file, header=TRUE)

        if (nrow(this.stat) > 0) {
            this.stat$chr = chr

            this.stat$pop1 = pops[1]
            this.stat$pop2 = pops[2]

            # CSS data.frame already doubled
        }
    }

    return(this.stat)
}

all.stat = list()
all.stat[["fst"]]   = do.call(rbind, lapply(stat.files[["fst"]],   read.stat))
all.stat[["xpehh"]] = do.call(rbind, lapply(stat.files[["xpehh"]], read.stat))
all.stat[["css"]]   = do.call(rbind, lapply(stat.files[["css"]],   read.stat))

# Standardize names
# Fst:
stat.col = which(names(all.stat[["fst"]]) == "WEIGHTED_FST")
pos.col  = which(names(all.stat[["fst"]]) == "BIN_START")
names(all.stat[["fst"]])[stat.col] = "statistic"
names(all.stat[["fst"]])[pos.col]  = "pos"

# XP-EHH:
stat.col = which(names(all.stat[["xpehh"]]) == "normxpehh")
pos.col  = which(names(all.stat[["xpehh"]]) == "pos")
names(all.stat[["xpehh"]])[stat.col] = "statistic"
names(all.stat[["xpehh"]])[pos.col]  = "pos"

# Set negative Fst values to zero
all.stat[["fst"]][all.stat[["fst"]]$statistic < 0,]$statistic = 0

all.stat[["fst"]]$stat.type   = "fst"
all.stat[["xpehh"]]$stat.type = "xpehh"
all.stat[["css"]]$stat.type   = "css"

# ----------------------------------------------------------------------------------------

# Reduce to just subset of populations:

to.include = c("KAZZI NSADZI", "BANDA KIYINDI", "BANDA NSADZI", "KAZZI KIYINDI")

all.stat[["fst"]] = all.stat[["fst"]][
    paste(all.stat[["fst"]]$pop1, all.stat[["fst"]]$pop2) %in% to.include,
]
all.stat[["xpehh"]] = all.stat[["xpehh"]][
    paste(all.stat[["xpehh"]]$pop1, all.stat[["xpehh"]]$pop2) %in% to.include,
]
all.stat[["css"]] = all.stat[["css"]][
    paste(all.stat[["css"]]$pop1, all.stat[["css"]]$pop2) %in% to.include,
]

# ----------------------------------------------------------------------------------------

pops = names(table(c(all.stat[["xpehh"]]$pop1, all.stat[["xpehh"]]$pop2)))
is.sites = c("BANDA", "BUGALA", "BUGALAIS", "BUKASA", "NSADZI", "SSERINYA")
ml.sites = pops[!pops %in% is.sites]

# ----------------------------------------------------------------------------------------

# Figure out which comparison class we're looking at
for (stat.type in c("fst", "xpehh", "css")) {
    island.island.boolean =     all.stat[[stat.type]]$pop1 %in% is.sites &
                                all.stat[[stat.type]]$pop2 %in% is.sites
    mainland.mainland.boolean = all.stat[[stat.type]]$pop1 %in% ml.sites &
                                all.stat[[stat.type]]$pop2 %in% ml.sites
    all.stat[[stat.type]]$comparison.type = "island.mainland"
    all.stat[[stat.type]][island.island.boolean,]$comparison.type     = "island.island"
    all.stat[[stat.type]][mainland.mainland.boolean,]$comparison.type = "mainland.mainland"
}

all.stat.to.plot = all.stat

###     # Negativize and scale Fst for plot
###     scaling.factor = max(all.stat[['xpehh']]$statistic, na.rm=TRUE) /
###                      max(all.stat[['fst']]$statistic,   na.rm=TRUE)
###     all.stat.to.plot[["fst"]]$statistic = -1 * all.stat[["fst"]]$statistic * scaling.factor

###     # Absolute value of XP-EHH for plot
###     all.stat.to.plot[["xpehh"]]$statistic = abs(all.stat[["xpehh"]]$statistic)

# Reduce Fst columns and combine datasets
all.stats.combined = rbind(all.stat.to.plot[["xpehh"]],
                           all.stat.to.plot[["fst"]][names(all.stat[["xpehh"]])],
                           all.stat.to.plot[["css"]][names(all.stat[["xpehh"]])])

# ----------------------------------------------------------------------------------------

all.stats = all.stats.combined[all.stats.combined$pop1 < all.stats.combined$pop2,]

all.stats.sm = na.omit(all.stats)
all.stats.sm$both.pops = paste(all.stats.sm$pop1, all.stats.sm$pop2, sep=" - ")

# Change Bugala labels
all.stats.sm$both.pops = gsub("BUGALAIS", "BUGALA (I)", all.stats.sm$both.pops)
all.stats.sm$both.pops = gsub("BUGALAML", "BUGALA (M)", all.stats.sm$both.pops)

all.stats.sm$both.pops = gsub("KAZZI",   "KAAZI",  all.stats.sm$both.pops)
all.stats.sm$both.pops = gsub("MITYANA", "WAMALA", all.stats.sm$both.pops)

comparisons = unique(all.stats.sm[c("pop1", "pop2", "both.pops", "comparison.type")])

all.stats.sm$both.pops = factor(all.stats.sm$both.pops,
    levels = comparisons$both.pops)

all.stats.sm$comparison.type = factor(all.stats.sm$comparison.type)

all.stats.sm$to.plot = FALSE

all.stats.sm = do.call(rbind, lapply(comparisons$both.pops, function (x) {

    this.pair = all.stats.sm[all.stats.sm$both.pops == x,]

    fst.cutoff = quantile(this.pair$statistic[this.pair$stat.type == "fst"],
        0.9, na.rm=TRUE)
    xpehh.cutoff = quantile(abs(this.pair$statistic[this.pair$stat.type == "xpehh"]),
        0.9, na.rm=TRUE)
    css.cutoff = quantile(this.pair$statistic[this.pair$stat.type == "css"],
        0.9, na.rm=TRUE)

    # Plot outliers
    this.pair[which(this.pair$stat.type == 'fst'   &
                    this.pair$statistic >= fst.cutoff),       ]$to.plot = TRUE
    this.pair[which(this.pair$stat.type == 'xpehh' &
                    abs(this.pair$statistic) >= xpehh.cutoff),]$to.plot = TRUE
    this.pair[which(this.pair$stat.type == 'css'   &
                    this.pair$statistic >= css.cutoff),       ]$to.plot = TRUE

    # Plot random points
    this.pair[sample(1:nrow(this.pair), size=5000),]$to.plot = TRUE

    this.pair
}))

# table(all.stats.sm$to.plot)

all.stats.sm[!all.stats.sm$to.plot,]$statistic = NA

# ----------------------------------------------------------------------------------------

# Remove unsplit Bugala
all.stats.sm = all.stats.sm[all.stats.sm$pop1 != "BUGALA",]
all.stats.sm = all.stats.sm[all.stats.sm$pop2 != "BUGALA",]

is.sites = c("BANDA", "BUGALA (I)", "BUKASA", "NSADZI", "SSERINYA")
ml.sites = c("BUGALA (M)", "BUWAMA", "KAZZI", "KIYINDI", "MITYANA")
all.stats.sm$pop1 = factor(all.stats.sm$pop1, levels = c("", is.sites, ml.sites))
all.stats.sm$pop2 = factor(all.stats.sm$pop2, levels = c("", is.sites, ml.sites))

all.stats.sm$chr = factor(all.stats.sm$chr, levels = chrs)

# ----------------------------------------------------------------------------------------
# --- Bring in location info for inversions, heterochromatic regions, and genes
# ----------------------------------------------------------------------------------------

y.max.xpehh = ceiling(max(abs(all.stat.to.plot[['xpehh']]$statistic), na.rm=TRUE))
y.min.xpehh = -1 * y.max.xpehh

y.max.fst = ceiling(max(abs(all.stat.to.plot[['fst']]$statistic) * 10,
    na.rm=TRUE)) / 10
y.min.fst = 0

y.max.css = ceiling(max(abs(all.stat.to.plot[['css']]$statistic), na.rm=TRUE))
y.min.css = 0

# --- Genes

genes = read.table("data/insecticide_genes_small.bed")
names(genes) = c("chr", "midpoint", "gene")

genes = do.call(rbind, lapply(1:nrow(comparisons),
    function (combo.idx) { cbind(comparisons[combo.idx,], genes) }))

genes$both.pops = factor(genes$both.pops,
    levels = levels(all.stats.sm$both.pops))
genes$chr = factor(genes$chr, levels = chrs)
genes$comparison.type = factor(genes$comparison.type,
    levels=levels(all.stats.sm$comparison.type))

# Duplicate for Fst, XP-EHH, and CSS
genes = rbind(
    cbind(genes, stat.type='fst',   ymin=y.min.fst,   ymax=y.max.fst),
    cbind(genes, stat.type='xpehh', ymin=y.min.xpehh, ymax=y.max.xpehh),
    cbind(genes, stat.type='css',   ymin=y.min.css,   ymax=y.max.css))

# --- Inversions

inversion_sites = read.table("data/inversion_simple.bed")
names(inversion_sites) = c("chr", "start", "end", "name")

inversion_sites = do.call(rbind, lapply(1:nrow(comparisons),
    function (combo.idx) { cbind(comparisons[combo.idx,], inversion_sites) }))

inversion_sites$both.pops = factor(inversion_sites$both.pops,
    levels = levels(all.stats.sm$both.pops))
inversion_sites$chr = factor(inversion_sites$chr, levels = chrs)

# Duplicate for Fst, XP-EHH, and CSS
inversion_sites = rbind(
    cbind(inversion_sites, stat.type='fst',   ymin=y.min.fst,   ymax=y.max.fst),
    cbind(inversion_sites, stat.type='xpehh', ymin=y.min.xpehh, ymax=y.max.xpehh),
    cbind(inversion_sites, stat.type='css',   ymin=y.min.css,   ymax=y.max.css))

# --- Heterochromatic

het = read.table("data/heterochromatin.bed")
names(het) = c("chr", "start", "end")

het = do.call(rbind, lapply(1:nrow(comparisons),
    function (combo.idx) { cbind(comparisons[combo.idx,], het) }))

het$both.pops = factor(het$both.pops,
    levels = levels(all.stats.sm$both.pops))
het$chr = factor(het$chr, levels = chrs)

# Duplicate for Fst, XP-EHH, and CSS
het = rbind(
    cbind(het, stat.type='fst',   ymin=y.min.fst,   ymax=y.max.fst),
    cbind(het, stat.type='xpehh', ymin=y.min.xpehh, ymax=y.max.xpehh),
    cbind(het, stat.type='css',   ymin=y.min.css,   ymax=y.max.css))

# --- New sweeps

new.sweeps = data.frame(chr=c("X", "2L"), loc=c(9238942, 34044820))

new.sweeps = do.call(rbind, lapply(1:nrow(comparisons),
    function (combo.idx) { cbind(comparisons[combo.idx,], new.sweeps) }))

new.sweeps$both.pops = factor(new.sweeps$both.pops,
    levels = levels(all.stats.sm$both.pops))
new.sweeps$chr = factor(new.sweeps$chr, levels = chrs)

# Duplicate for Fst, XP-EHH, and CSS
new.sweeps = rbind(
    cbind(new.sweeps, stat.type='fst',   ymin=y.min.fst,   ymax=y.max.fst),
    cbind(new.sweeps, stat.type='xpehh', ymin=y.min.xpehh, ymax=y.max.xpehh),
    cbind(new.sweeps, stat.type='css',   ymin=y.min.css,   ymax=y.max.css))

# --- Accessibility

#   #acc = read.table(paste0("data/accessibility/accessibility.", chr, ".vcf"),
#   #    colClasses = c("NULL", "integer", rep("NULL", 4), "character", "character"))
#
#   acc.cmd = "grep '^X' data/accessibility/accessibility.X.vcf | grep 'Accessible' | cut -f 2"
#   inacc.cmd = "grep '^X' data/accessibility/accessibility.X.vcf | grep -v 'Accessible' | cut -f 2"
#
#   acc = read.table(pipe(acc.cmd), colClasses = c("integer"))
#   inacc = read.table(pipe(inacc.cmd), colClasses = c("integer"))
#
#   acc$is.accessible = TRUE
#   inacc$is.accessible = TRUE
#
#   this.acc = rbind(acc, inacc)
#
#   names(this.acc)[1] = "pos"

# ----------------------------------------------------------------------------------------
# --- Gene labels
# ----------------------------------------------------------------------------------------

gene.label.text = rbind(
    data.frame(win.center = 28497407, lab = "Cyp6p",
        chr = factor("2R", levels = chrs)),
    data.frame(win.center =  2394888, lab = "Vgsc",
        chr = factor("2L", levels = chrs)),
    data.frame(win.center = 25399104, lab = "Gaba",
        chr = factor("2L", levels = chrs)),
    data.frame(win.center = 28598038, lab = "Gste",
        chr = factor("3R", levels = chrs)),
    data.frame(win.center = 11204486, lab = "Tep1",
        chr = factor("3L", levels = chrs)),
    data.frame(win.center = 15241718, lab = "Cyp9k1",
        chr = factor("X",  levels = chrs))
)

label.height = 4.5
label.expand.factor = 0.05

gene.label.text$stat.val = label.height
# Figure out top row in plot
top.row = head(all.stats.sm[order(all.stats.sm$comparison.type, all.stats.sm$both.pops),],
    n=1)
gene.label.text$both.pops = top.row$both.pops
gene.label.text$comparison.type = top.row$comparison.type
gene.label.text$stat.type = "css"

gene.label.fake.pts = gene.label.text
gene.label.fake.pts$win.center = 0
gene.label.fake.pts$stat.val = seq(from = label.height - label.expand.factor,
                                   to   = label.height + label.expand.factor,
                                   length.out = 6)

novel.label.text = data.frame(win.center=c(9238942, 34044820), lab="",
    chr=factor(c("X", "2L"), levels = chrs),
    stat.type="css",
    stat.val=label.height,
    both.pops = top.row$both.pops)
novel.label.text = merge(novel.label.text,
    comparisons[,c("both.pops", "comparison.type")])

# For drawing in fake y-axis
axis.df = data.frame(x=0, xend=0,
    both.pops = factor(levels(all.stats.sm$both.pops),
                       levels = levels(all.stats.sm$both.pops)),
    chr = factor("2R", levels = chrs))

axis.df = rbind(
    cbind(axis.df, stat.type='fst',   y=y.min.fst,   yend=y.max.fst),
    cbind(axis.df, stat.type='xpehh', y=y.min.xpehh, yend=y.max.xpehh),
    cbind(axis.df, stat.type='css',   y=y.min.css,   yend=y.max.css))

axis.df = merge(axis.df, comparisons[,c("both.pops", "comparison.type")])

# ----------------------------------------------------------------------------------------
# --- Make overall combined plot for text
# ----------------------------------------------------------------------------------------

facet.labeller.sites = function(value) {
    return(gsub(" - ", "\n", value))
}

facet.labeller.stat = as_labeller(c(xpehh="XP-EHH", fst="F[ST]", css="CSS"), label_parsed)

facet.labeller.isml = as_labeller(c(island.island     = "Between\nIslands",
                                    island.mainland   = "Island vs.\nMainland",
                                    mainland.mainland = "Between\nMainland Sites"))

p = ggplot(all.stats.sm, aes(pos, statistic, color=comparison.type, group=stat.type)) +
    geom_rect(data=inversion_sites,
        aes(x=NULL, y=NULL, col=NULL,
            xmin=start, xmax=end, ymin=ymin, ymax=ymax),
        alpha=0.2, fill='lightgrey') +
    geom_vline(data=genes, aes(xintercept = midpoint), lty=3, alpha=0.5) +
    geom_vline(data=new.sweeps, aes(xintercept = loc),
        lty=1, alpha=1, col="black", lwd=0.25) +
    geom_point(size=1e-10) +
    geom_rect(data=het,
        aes(x=NULL, y=NULL, col=NULL,
            xmin=start, xmax=end, ymin=ymin, ymax=ymax, both.pops=both.pops),
        alpha=1, fill='grey') +
    geom_segment(data=axis.df, aes(x=x, y=y, xend=xend, yend=yend), col="black", lwd=1) +
    scale_color_manual(values=c("#00BFC4", "#C060D6", "#F8766D"),
        drop=TRUE, limits = levels(all.stats.sm$comparison.type)) +
    geom_hline(yintercept=0, col="black", lwd=0.1) +
    facet_grid(comparison.type + both.pops + stat.type ~ chr,
        scales = "free", space = "free_x", shrink=FALSE, drop=TRUE,
        labeller = labeller(stat.type       = facet.labeller.stat,
                            both.pops       = facet.labeller.sites,
                            comparison.type = facet.labeller.isml)) +
    xlab("Genome position (Mb)") +
    ylab(expression("CSS,"~F[ST]~", and XP-EHH")) +
    guides(color=FALSE) +
    scale_x_continuous(labels = function(x) { x/1000000 },
                       breaks = seq(from=0e6, to=60e6, by=10e6),
                       expand = c(0,0)) +
    scale_y_continuous(expand=c(0,0), limits=c(NA, NA)) +
    theme_bw() +
    theme(strip.text = element_text(size=7, lineheight=1.0),
          strip.background = element_rect(fill="white", color="black", size=0.5),
          strip.text.y = element_text(angle = 0),
          panel.spacing.y = unit(0.75, "lines"),
          axis.line.x = element_line(color = "black"),
          axis.line.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    geom_point(data = gene.label.fake.pts, mapping = aes(x = win.center, y = stat.val),
               pch=NA, col=NA, fill=NA) +
    geom_point(data = gene.label.text, mapping = aes(x = win.center, y = stat.val),
               pch=25, col='black', fill='black') +
    geom_text(data = gene.label.text,
              mapping = aes(x = win.center + 1000000, y = stat.val, label = lab),
              color="black", fontface="italic", size=2, hjust = 'left') +
    geom_point(data = novel.label.text, mapping = aes(x = win.center, y = stat.val),
               pch=25, col='black', fill='white')

ggsave("reports/fst_xpehh_by_site.pdf", p, height=8, width=12)

#   p.ii = p %+% droplevels(subset(all.stats.sm,
#       all.stats.sm$comparison.type == "island.island"))
#   p.im = p %+% droplevels(subset(all.stats.sm,
#       all.stats.sm$comparison.type == "island.mainland"))
#   p.mm = p %+% droplevels(subset(all.stats.sm,
#       all.stats.sm$comparison.type == "mainland.mainland"))
#
#   ggsave("reports/fst_xpehh_by_site.island-island.pdf",     p.ii, height=8, width=12)
#   ggsave("reports/fst_xpehh_by_site.island-mainland.pdf",   p.im, height=8, width=12)
#   ggsave("reports/fst_xpehh_by_site.mainland-mainland.pdf", p.mm, height=8, width=12)
