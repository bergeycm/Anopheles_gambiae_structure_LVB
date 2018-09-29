#!/usr/bin/env Rscript

# ========================================================================================
# --- Plot Fst against pi
# ========================================================================================

options(stringsAsFactors=FALSE)

library(ggplot2)

# Swap Group1 and Group2 if they won't be in order after the re-ordering
islands = c("BANDA", "BUGALA", "BUGALAIS", "BUKASA", "NSADZI", "SSERINYA")
mainland = c("BUGALAML", "BUWAMA", "KAZZI", "KIYINDI", "MITYANA")

all.fst.files = list.files(path="results/", pattern="^chr.*.windowed.weir.fst$")

pop.pairs = unique(gsub(".*\\.(.*)\\.windowed.weir.fst", "\\1", all.fst.files))

all.fst = do.call(rbind, lapply(pop.pairs, function(pop.pair) {
    this.pair.fst.files = all.fst.files[grep(pop.pair, all.fst.files)]
    this.fst.df = do.call(rbind,
        lapply(this.pair.fst.files, function(f) {
            this.fst = read.table(paste0("results/", f), header=TRUE, nrow=1000)
            this.fst$pops = pop.pair
            this.fst
        })
    )
}))

all.fst$pop1 = do.call(rbind, strsplit(all.fst$pops, "_"))[,1]
all.fst$pop2 = do.call(rbind, strsplit(all.fst$pops, "_"))[,2]

# ----------------------------------------------------------------------------------------

all.pi.files = list.files(path="results/", pattern="^chr.*.windowed.pi$")

pops = unique(gsub(".*\\.(.*)\\.windowed.pi", "\\1", all.pi.files))

all.pi = do.call(rbind, lapply(pops, function(pop) {
    this.pop.pi.files = all.pi.files[grep(pop, all.pi.files)]
    this.pi.df = do.call(rbind,
        lapply(this.pop.pi.files, function(f) {
            this.pi = read.table(paste0("results/", f), header=TRUE, nrow=1000)
            this.pi$pop = pop
            this.pi
        })
    )
}))

all.fst.pi = merge(
    merge(all.fst, all.pi,
        by.x=c(names(all.fst)[1:3], "pop1"),
        by.y=c(names(all.fst)[1:3], "pop"),
        suffixes=c("", ".PI.pop1")), all.pi,
    by.x=c(names(all.fst)[1:3], "pop2"),
    by.y=c(names(all.fst)[1:3], "pop"),
    suffixes=c("", ".PI.pop2"))

names(all.fst.pi) = gsub("PI.PI.", "PI.", names(all.fst.pi))

# ----------------------------------------------------------------------------------------

all.fst.pi$is.island.pop1 = all.fst.pi$pop1 %in% islands
all.fst.pi$is.island.pop2 = all.fst.pi$pop2 %in% islands

# If pop2 is an island (listed first) and pop1 is a mainland site, switch
to.switch = all.fst.pi$is.island.pop2 & (!all.fst.pi$is.island.pop1)
pop1.backup = all.fst.pi$pop1
all.fst.pi[to.switch,]$pop1 = all.fst.pi$pop2[to.switch]
all.fst.pi[to.switch,]$pop2 = pop1.backup[to.switch]

all.fst.pi$comparison.type = "BOTH_MAINLAND"
all.fst.pi[all.fst.pi$is.island.pop1 | all.fst.pi$is.island.pop2,]$comparison.type =
    "ISLAND_MAINLAND"
all.fst.pi[all.fst.pi$is.island.pop1 & all.fst.pi$is.island.pop2,]$comparison.type =
    "BOTH_ISLAND"

# Order sites, islands first then mainland
all.fst.pi$pop1 = factor(all.fst.pi$pop1, levels=c(islands, mainland))
all.fst.pi$pop2 = factor(all.fst.pi$pop2, levels=c(islands, mainland))

# ----------------------------------------------------------------------------------------
# --- Remove inversions and heterochromatic regions
# ----------------------------------------------------------------------------------------

het = read.table("data/heterochromatin.bed")
names(het) = c("chr", "start", "end")

inv = read.table("data/inversion_sites.bed")
names(inv) = c("chr", "start", "end", "name")
inv$inversion = gsub("_.*", "", inv$name)

# Make a guess at where the distal breakpoint of 2Rb is
distal.2Rb.guess = c("2R", 18575300, 18575300, "2Rb_dist", "2Rb")
inv = rbind(inv, distal.2Rb.guess)

inv$start = as.numeric(inv$start)
inv$end   = as.numeric(inv$end)

full.coord = list()

for (inversion in unique(inv$inversion)) {

    this.chr   = inv$chr  [inv$inversion == inversion][1]
    this.start = inv$start[inv$inversion == inversion]
    this.end   = inv$end  [inv$inversion == inversion]

    full.coord[[inversion]] = c(min(this.start, this.end),
                                max(this.start, this.end))
}

# ----------------------------------------------------------------------------------------

# --- Mask out heterochromatic stuff

all.fst.pi$is.het = FALSE

for (this.het.idx in 1:nrow(het)) {

    this.het = het[this.het.idx,]

    het.bool = all.fst.pi$CHROM == this.het$chr &
        all.fst.pi$POS >= this.het$start &
        all.fst.pi$POS <= this.het$end

    if (sum(het.bool)) {
        all.fst.pi[het.bool,]$is.het = TRUE
    }
}

# --- Mask out inversion regions

all.fst.pi$inv.state = "collinear"

for (this.inv.idx in 1:length(full.coord)) {

    this.inv = full.coord[this.inv.idx]
    this.chr = substr(names(this.inv), start=0, stop=2)

    inv.bool = all.fst.pi$CHROM == this.chr &
        all.fst.pi$POS >= this.inv[[1]][1] &
        all.fst.pi$POS <= this.inv[[1]][2]

    if (sum(inv.bool)) {
        all.fst.pi[inv.bool,]$inv.state = names(this.inv)
    }
}

# Create boolean for Fst values we want: no heterochromatic stuff or inversion regions
all.fst.pi$to.include = all.fst.pi$inv.state == "collinear" & (!all.fst.pi$is.het)

all.fst.pi = all.fst.pi[all.fst.pi$to.include,]

# ----------------------------------------------------------------------------------------
# --- Reduce by sampling SNPs
# ----------------------------------------------------------------------------------------

sample.size = 1000

snps = unique(all.fst.pi[,1:2])
snps.samp = snps[sample(1:nrow(snps), size=sample.size),]

snps.samp$snp.str = paste(snps.samp$CHROM, snps.samp$BIN_START, sep="_")

# Backup copy
all.fst.pi.orig = all.fst.pi

all.fst.pi = all.fst.pi[paste(all.fst.pi$CHROM, all.fst.pi$BIN_START, sep="_") %in%
    snps.samp$snp.str,]

# ----------------------------------------------------------------------------------------
# --- Duplicate data, switching pop1 and pop2 to be able to plot against pi of both
# ----------------------------------------------------------------------------------------

all.fst.pi.part1 = all.fst.pi[,c(4,5,7,8,10,11,16,19)]
all.fst.pi.part2 = all.fst.pi[,c(5,4,7,8,12,13,16,19)]

names(all.fst.pi.part1) = names(all.fst.pi.part2) = c("pop1", "pop2",
    "WEIGHTED_FST", "MEAN_FST", "N_VARIANTS.PI", "PI",
    "comparison.type", "to.include")

all.fst.pi.full = rbind(all.fst.pi.part1, all.fst.pi.part2)

# ----------------------------------------------------------------------------------------
# --- Do scatter plot
# ----------------------------------------------------------------------------------------

p = ggplot(subset(all.fst.pi.full, to.include),
        aes(PI, WEIGHTED_FST, color=comparison.type)) +
    geom_point(size=0.2) +
    facet_grid(pop1 ~ pop2, drop=FALSE) +
    scale_color_discrete(guide=FALSE) +
    theme_bw() +
    theme(strip.text = element_text(size = 7),
          axis.text  = element_text(size=5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab(expression(pi)) + ylab(expression("F"["ST"]))

ggsave(p, filename="reports/fst_by_pi.pdf", height=6, width=8)
