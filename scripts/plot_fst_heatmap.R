#!/usr/bin/env Rscript

# ========================================================================================
# --- Do various plots of Fst
# ========================================================================================

options(stringsAsFactors=FALSE)

library(ggplot2)
library(reshape)

# Swap Group1 and Group2 if they won't be in order after the re-ordering
islands = c("BANDA", "BUGALA", "BUGALAIS", "BUKASA", "NSADZI", "SSERINYA")
mainland = c("BUGALAML", "BUWAMA", "KAZZI", "KIYINDI", "MITYANA")

# ========================================================================================
# --- Plot heatmap of pairwise Fst
# ========================================================================================

filters = c("no_het_or_inv", "no_het_2La", "no_het_2Rb", "no_het_2Rc")

fst = read.table("reports/fst_summary.txt", header=TRUE)

# Set negative values of Fst to zero
fst[which(fst$value < 0),]$value = 0

# Remove NaN's
fst = fst[which(fst$value != "NaN"),]

plot.fst.heatmap = function(filter, color) {

    # Mean Fst
    fst.mean = fst[which(fst$filter == filter & fst$stat == "mean"),]

    # Remove NA values
    fst.mean = fst.mean[!is.na(fst.mean$Group2),]

    # Fix labels for when Fst = 0
    labels = sprintf("%.3f", round(fst.mean$value, digits=3))

    return.idx = function(x) { which(x == c(islands, mainland)) }
    group1.idx = do.call(c, lapply(fst.mean$Group1, return.idx))
    group2.idx = do.call(c, lapply(fst.mean$Group2, return.idx))

    grp1.backup = fst.mean$Group1
    fst.mean[group2.idx < group1.idx,]$Group1 = fst.mean$Group2[group2.idx < group1.idx]
    fst.mean[group2.idx < group1.idx,]$Group2 = grp1.backup[group2.idx < group1.idx]

    # Change order to separate islands and mainland
    fst.mean$Group1 = factor(fst.mean$Group1, levels = c(islands, mainland))
    fst.mean$Group2 = factor(fst.mean$Group2, levels = c(islands, mainland))

    p = ggplot(fst.mean, aes(Group1, Group2)) +
        geom_tile(aes(fill = value), color = "black") +
        geom_text(data=fst.mean,
            aes(Group1, Group2,
               label = labels),
            color=c(color, "white")[1 + as.numeric(fst.mean$value > 0.01)],
            size=rel(5)) +
        scale_fill_gradient(low = "white", high = color) +
        coord_fixed() +
        labs(x='', y='') +
        theme_grey(base_size=9) +
        theme(panel.background = element_blank(),
            panel.grid.major=element_blank(),
            legend.position = "none",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.text=element_text(size=12, face="bold")) +
        scale_y_discrete(limits = rev(levels(factor(fst.mean$Group2))))

    # Add lines to separate islands and mainland
    p = p +
        geom_segment(x=length(islands) + 0.5,
                     xend=length(islands) + 0.5,
                     y=0.5,
                     yend=length(mainland) + 0.5,
                     col='black', lwd=2) +
        geom_segment(x=0.5,
                     xend=length(islands) + 0.5,
                     y=length(mainland) + 0.5,
                     yend=length(mainland) + 0.5,
                     col='black', lwd=2)

    ggsave(p, filename=paste0("reports/fst_heatmap_", filter, ".pdf"))

    return(p)
}

plot.fst.heatmap("no_het_or_inv", "black")
#plot.fst.heatmap("no_het_2La",    "darkred")
#plot.fst.heatmap("no_het_2Rb",    "darkgreen")

# ========================================================================================
# --- Plot Fst distributions
# ========================================================================================

all.fst.files = list.files(path="results/", pattern=".*3.*.fst$")
all.fst.files = all.fst.files[grep("window", all.fst.files, invert=TRUE)]

pop.pairs = unique(gsub(".*\\.(.*)\\.weir.fst", "\\1", all.fst.files))

all.fst = do.call(rbind, lapply(pop.pairs, function(pop.pair) {
    this.pair.fst.files = all.fst.files[grep(pop.pair, all.fst.files)]
    this.fst.df = do.call(rbind,
        lapply(this.pair.fst.files, function(f) {
            this.fst = read.table(paste0("results/", f), header=TRUE)
            this.fst$pops = pop.pair
            this.fst
        })
    )
}))

all.fst$pop1 = do.call(rbind, strsplit(all.fst$pops, "_"))[,1]
all.fst$pop2 = do.call(rbind, strsplit(all.fst$pops, "_"))[,2]

all.fst$is.island.pop1 = all.fst$pop1 %in% islands
all.fst$is.island.pop2 = all.fst$pop2 %in% islands

# If pop2 is an island (listed first) and pop1 is a mainland site, switch
to.switch = all.fst$is.island.pop2 & (!all.fst$is.island.pop1)
pop1.backup = all.fst$pop1
all.fst[to.switch,]$pop1 = all.fst$pop2[to.switch]
all.fst[to.switch,]$pop2 = pop1.backup[to.switch]

all.fst$comparison.type = "BOTH_MAINLAND"
all.fst[all.fst$is.island.pop1 | all.fst$is.island.pop2,]$comparison.type =
    "ISLAND_MAINLAND"
all.fst[all.fst$is.island.pop1 & all.fst$is.island.pop2,]$comparison.type =
    "BOTH_ISLAND"

# Order sites, islands first then mainland
all.fst$pop1 = factor(all.fst$pop1, levels=c(islands, mainland))
all.fst$pop2 = factor(all.fst$pop2, levels=c(islands, mainland))

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

all.fst$is.het = FALSE

for (this.het.idx in 1:nrow(het)) {

    this.het = het[this.het.idx,]

    het.bool = all.fst$CHROM == this.het$chr &
        all.fst$POS >= this.het$start &
        all.fst$POS <= this.het$end

    if (sum(het.bool)) {
        all.fst[het.bool,]$is.het = TRUE
    }
}

# --- Mask out inversion regions

all.fst$inv.state = "collinear"

for (this.inv.idx in 1:length(full.coord)) {

    this.inv = full.coord[this.inv.idx]
    this.chr = substr(names(this.inv), start=0, stop=2)

    inv.bool = all.fst$CHROM == this.chr &
        all.fst$POS >= this.inv[[1]][1] &
        all.fst$POS <= this.inv[[1]][2]

    if (sum(inv.bool)) {
        all.fst[inv.bool,]$inv.state = names(this.inv)
    }
}

# Create boolean for Fst values we want: no heterochromatic stuff or inversion regions
all.fst$to.include = all.fst$inv.state == "collinear" & (!all.fst$is.het)

# ----------------------------------------------------------------------------------------
# --- Reduce by sampling SNPs
# ----------------------------------------------------------------------------------------

sample.size = 1000

snps = unique(all.fst[,1:2])
snps.samp = snps[sample(1:nrow(snps), size=sample.size),]

snps.samp$snp.str = paste(snps.samp$CHROM, snps.samp$POS, sep="_")

# Backup copy
###all.fst.orig = all.fst

all.fst = all.fst[paste(all.fst$CHROM, all.fst$POS, sep="_") %in%
    snps.samp$snp.str,]

save.image("fst_heatmap.Rdata")

# ----------------------------------------------------------------------------------------
# --- Bring in medians
# ----------------------------------------------------------------------------------------

# Get medians together for text annotations

# Weir and Cockerham method
###fst.mean = fst[which(fst$filter == "no_het_or_inv" & fst$stat == "mean"),]
###names(fst.mean) = c("pop1", "pop2", "species", "filer", "stat", "value")

# Hudson method - Fst
hud.fst = read.table("data/all.pass.snp.flt.eigen.fst.se.out.fst.txt", row.names=NULL)
hud.fst$row.names = names(hud.fst)[-c(1)]
hud.fst.l = melt(hud.fst, id.vars="row.names")
names(hud.fst.l) = c("pop1", "pop2", "value")
hud.fst.l$value = hud.fst.l$value / 1000
hud.fst.l$stat = "mean"

# Hudson method - SD
hud.sd = read.table("data/all.pass.snp.flt.eigen.fst.se.out.sd.txt", row.names=NULL)
hud.sd$row.names = names(hud.sd)[-c(1)]
hud.sd.l = melt(hud.sd, id.vars="row.names")
names(hud.sd.l) = c("pop1", "pop2", "value")
hud.sd.l$value = hud.sd.l$value / 1000000
hud.sd.l$stat = "mean"

fst.mean = rbind(hud.fst.l, hud.sd.l)

# Fix labels for when Fst = 0
fst.mean$labels = ""
fst.mean[fst.mean$stat == "mean",]$labels =
    sprintf("%.3f",  round(fst.mean[fst.mean$stat == "mean",]$value, digits=3))
fst.mean[fst.mean$stat == "sd",  ]$labels =
    sprintf("%.1e", round(fst.mean[fst.mean$stat == "sd",  ]$value, digits=3))

# Add fake data for plotting text in correct place
fst.mean$WEIR_AND_COCKERHAM_FST = (max(all.fst$WEIR_AND_COCKERHAM_FST) +
    min(all.fst$WEIR_AND_COCKERHAM_FST)) / 2
fst.mean$comparison.type="ISLAND_MAINLAND"
fst.mean$scaled = 0.50

fst.mean$is.island.pop1 = fst.mean$pop1 %in% islands
fst.mean$is.island.pop2 = fst.mean$pop2 %in% islands

# If pop2 is an island (listed first) and pop1 is a mainland site, switch
to.switch = fst.mean$is.island.pop2 & (!fst.mean$is.island.pop1)
pop1.backup = fst.mean$pop1
fst.mean[to.switch,]$pop1 = fst.mean$pop2[to.switch]
fst.mean[to.switch,]$pop2 = pop1.backup[to.switch]

fst.mean$pop1 = factor(fst.mean$pop1, levels=c(islands, mainland))
fst.mean$pop2 = factor(fst.mean$pop2, levels=c(islands, mainland))

fst.mean$comparison.type = "BOTH_MAINLAND"
fst.mean[fst.mean$is.island.pop1 | fst.mean$is.island.pop2,]$comparison.type =
    "ISLAND_MAINLAND"
fst.mean[fst.mean$is.island.pop1 & fst.mean$is.island.pop2,]$comparison.type =
    "BOTH_ISLAND"

# Switch pop1 and pop2 to plot text in the lower diagonal
names(fst.mean)[1:2] = paste0("pop", 2:1)

# ----------------------------------------------------------------------------------------
# --- Add vertical line at origin
# ----------------------------------------------------------------------------------------

# Create additional dataframe just to put vertical line at upper diagonal zero
origin.line.df = unique(all.fst[,c(4:7)])

# ----------------------------------------------------------------------------------------
# --- Do distribution plot
# ----------------------------------------------------------------------------------------

p = ggplot(subset(all.fst, to.include),
        aes(WEIR_AND_COCKERHAM_FST, ..scaled.., color=comparison.type)) +
    geom_vline(data=origin.line.df, aes(xintercept=0), color="black", lty=3, lwd=0.2) +
    geom_density() +
    facet_grid(pop1 ~ pop2, drop=FALSE) +
    theme_bw() +
    theme(strip.text = element_text(size = 7),
          axis.text  = element_text(size=5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab(expression("F"["ST"])) + ylab("Scaled Density") +
    geom_text(data = fst.mean, aes(label = labels, color=comparison.type, scaled=0.5)) +
    scale_color_manual(guide=FALSE, values=c("#00B6EB", "#F8766D", "#53B400"))

ggsave(p, filename="reports/fst_distributions.pdf", height=6, width=8)
