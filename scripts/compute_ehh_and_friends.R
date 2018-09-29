#!/usr/bin/env Rscript

# ========================================================================================
# --- Compute various stats with rehh
# ========================================================================================

options(stringsAsFactors=FALSE)

library(rehh)
library(dplyr)
library(MASS)

source("scripts/tweaked_rehh_functions.R")
source("scripts/bifurcation_plot_andrewparkermorgan.R")

environment(calc_ehh)            = asNamespace('rehh')
environment(ihsplot)             = asNamespace('rehh')
environment(bifurcation.diagram) = asNamespace('rehh')

write("Finished loading packages", stderr())

args = commandArgs(trailingOnly = TRUE)

chr                = args[1]
center             = args[2]
width              = args[3]    # Actually 0.5x width
nickname           = args[4]    # Short name for writing duplicate copy of EHH plot
do.reverse.anc.der = args[5]    # Should we flip ancestral and derived?
# Test: chr="X"; center=9250808; width=10000; nickname="X_9Mb"; do.reverse.anc.der=FALSE
#       chr="2L"; center=34044820; width=100000; nickname="2L_34Mb"; do.reverse.anc.der=FALSE

do.checkpointing = FALSE

# Do full chromosome if no window info given
if (length(args) == 1) {
    center = 0
    width = 0
    nickname = "sweep"
    do.reverse.anc.der = FALSE
} else if (length(args) == 5) {
    center = as.numeric(center)
    width = as.numeric(width)
    do.reverse.anc.der = as.logical(do.reverse.anc.der)
} else {
    stop("ERROR: Must pass either 1 or 5 arguments.")
}

min.bp = center - width
max.bp = center + width

in.hap = paste0("data/chr", chr, ".pass.snp.phased.haps")

tmpdir = "tmp_ehh"
dir.create("tmp_ehh", showWarnings=FALSE)
in.map = tempfile(pattern = paste0("chr", chr, "_"), tmpdir = tmpdir, fileext = ".map")
in.hap.tmp = gsub(".map", ".tmp.hap", in.map)

# --- Fake map file
# Since we lack map for Anopheles, assume 1Mb = 1 cM
if (center != 0) {
    awk.cmd = paste0("awk '{ if ($3 > ", min.bp, " && $3 < ", max.bp,
        ")  printf \"%s %s %d %s %s\\n\", $3,$1,$3,$4,$5 }' ",
        in.hap, " > ", in.map)
} else {
    awk.cmd = paste0("awk '{ printf \"%s %s %d %s %s\\n\", $3,$1,$3,$4,$5 }' ",
        in.hap, " > ", in.map)
}
system(awk.cmd)

write("Finished writing map file", stderr())

# ----------------------------------------------------------------------------------------

samp = read.table(paste0("data/chr", chr, ".pass.snp.phased.sample"),
    header = TRUE, stringsAsFactors = FALSE)

pop.inds = read.table("data/ssese_individual_info_simple.txt")

pop.inds.ordered = samp$ID_1[samp$ID_1 %in% pop.inds$V2]

# --- Code to transform hap file modified from:
# --- https://github.com/ksamuk/rehh_helper/blob/master/run.R

# Reduce to just region of interest
if (center != 0) {
    awk.cmd = paste0("awk '{ if ($3 > ", min.bp, " && $3 < ", max.bp,
        ")  print $0 }' ", in.hap, " > ", in.hap.tmp)
} else {
    # Lazy and inefficient way to copy a file while avoiding renaming variable
    awk.cmd = paste0("awk '{ print $0 }' ", in.hap, " > ", in.hap.tmp)
}
system(awk.cmd)

hap = read.table(in.hap.tmp)

write("Finished reading in hap file", stderr())

in.hap.m = hap[,-(1:5)] %>% as.matrix

row.names(in.hap.m) = hap[,3]

hap.names = rep(pop.inds.ordered, each = 2)

colnames(in.hap.m) = hap.names
in.hap.m = in.hap.m + 1
in.hap.m = t(in.hap.m)

# Annotate rows with row.names
haps.1.indices = seq(from = 1, to = length(hap.names), by = 2)
haps.2.indices = seq(from = 2, to = length(hap.names), by = 2)

hap.names[haps.1.indices] = paste0(hap.names[haps.1.indices],"_1")
hap.names[haps.2.indices] = paste0(hap.names[haps.2.indices],"_2")

write("Finished rejiggering data frame", stderr())

in.hap.m.names = as.matrix(cbind(row.names(in.hap.m), as.data.frame(in.hap.m)))
colnames(in.hap.m.names) = NULL

if (do.checkpointing) {
    if (center != 0) {
        save.image(paste0("after_transform_chr", chr, ".",
            min.bp, ".", max.bp, ".Rdata"))
    } else {
        save.image(paste0("after_transform_chr", chr, ".",
            "full.Rdata"))
    }
}

write("Finished transform", stderr())

# Write to file
if (center != 0) {
    in.hap.fix = gsub(".haps$", paste0(".", min.bp, ".", max.bp, ".fixed.haps"),
        gsub("data", "results", in.hap))
} else {
    in.hap.fix = gsub(".haps$", ".fixed.haps",
        gsub("data", "results", in.hap))
}

# Reverse ancestral and derived
if (do.reverse.anc.der) {
    in.hap.m.names[in.hap.m.names == "1"] = "orig_1"
    in.hap.m.names[in.hap.m.names == "2"] = "1"
    in.hap.m.names[in.hap.m.names == "orig_1"] = "2"
}

write.matrix(in.hap.m.names, file = in.hap.fix, sep="\t")

if (do.checkpointing) {
    if (center != 0) {
        save.image(paste0("after_write_chr", chr, ".",
            min.bp, ".", max.bp, ".Rdata"))
    } else {
        save.image(paste0("after_write_chr", chr, ".",
            "full.Rdata"))
    }
}

# ----------------------------------------------------------------------------------------

# --- Read in haps file

hap = data2haplohh(hap_file = in.hap.fix,
                   map_file = in.map,
                   haplotype.in.columns = FALSE,
                   recode.allele = FALSE,
                   min_perc_geno.hap = 85,
                   chr.name = chr)

write("Finished reading in hap file", stderr())

# Caution: Parallelized
ehh.scan = scan_hh(hap, threads=8)

write("Finished EHH scan", stderr())

ihs = ihh2ihs(ehh.scan, freqbin=0.05)

# To avoid trying to plot things with color "X", etc.
chrs = c("2L", "2R", "3L", "3R", "X")
ihs$iHS$CHR = which (chr == chrs)

write("Finished iHS computation", stderr())

dir.create("results/ehh/", showWarnings=FALSE)
results.dir = paste0("results/ehh/chr", chr, ".pass.snp.phased.figs")
dir.create(results.dir, showWarnings=FALSE)

setwd(results.dir)

write("Moved into results directory", stderr())

pdf(file=paste0("iHS_chr", chr, ".", min.bp, "-", max.bp, "_%01d.pdf"), onefile=FALSE)
    ihsplot(ihs, plot.pval=TRUE, ylim.scan=2, main="iHS")
dev.off()

write("Finished iHS plot", stderr())

# Calculate EHH at a given SNP

#   # Previous method grabbed SNP with highest value of iHS
#   most.extreme = which(abs(ihs$iHS$iHS) == max(abs(ihs$iHS$iHS), na.rm=TRUE))
#   most.extreme.val = ihs$iHS[most.extreme,]
#   most.extreme.hap.idx = which(most.extreme.val$POSITION == hap@snp.name)

# Now grab the one closest to passed target
dist.to.center = abs(center - ihs$iHS$POSITION)
# Grab first if tie
most.extreme = which(dist.to.center == min(dist.to.center))[1]
most.extreme.val = ihs$iHS[most.extreme,]
most.extreme.hap.idx = which(most.extreme.val$POSITION == hap@snp.name)

pdf(file=paste0("ehh_chr", chr, ".", min.bp, "-", max.bp, "_%01d.pdf"), onefile=FALSE,
        width=3, height=3)
    res.ehh = calc_ehh(hap, most.extreme.hap.idx, col=c("grey", "darkred"),
        cex.lab=1.5, cex.main=1.5)
dev.off()

# Also write to file without positioning, so Snakefile can recognize it if peaks change
pdf(file=paste0("ehh_chr", chr, ".", nickname, "_%01d.pdf"), onefile=FALSE,
        width=3, height=3)
    res.ehh = calc_ehh(hap, most.extreme.hap.idx, col=c("grey", "darkred"),
        cex.lab=1.5, cex.main=1.5)
dev.off()

write("Finished EHH plot", stderr())

# Plot bifurcation diagram for both ancestral and derived allele of most extreme SNP
pdf(file=paste0("bifurcation_chr", chr, ".", min.bp, "-", max.bp, "_%01d.pdf"),
        onefile=FALSE, width=3, height=3)
    # Ancestral allele, hypothetically
    bifurcation.diagram(hap, most.extreme.hap.idx,
        limhapcount=1, all_foc=1, nmrk_l=200, nmrk_r=200,
        linecol="grey", main_leg="Ancestral Allele")
    # Derived allele, hypothetically
    bifurcation.diagram(hap, most.extreme.hap.idx,
        limhapcount=1, all_foc=2, nmrk_l=200, nmrk_r=200,
        linecol="darkred", main_leg="Derived Allele")

dev.off()

# Also write to file without positioning, so Snakefile can recognize it if peaks change
pdf(file=paste0("bifurcation_chr", chr, ".", nickname, "_%01d.pdf"),
        onefile=FALSE, width=3, height=3)
    # Ancestral allele, hypothetically
    bifurcation.diagram(hap, most.extreme.hap.idx,
        limhapcount=1, all_foc=1, nmrk_l=200, nmrk_r=200,
        linecol="grey", main_leg="Ancestral Allele")
    # Derived allele, hypothetically
    bifurcation.diagram(hap, most.extreme.hap.idx,
        limhapcount=1, all_foc=2, nmrk_l=200, nmrk_r=200,
        linecol="darkred", main_leg="Derived Allele")
dev.off()

#   # Version with ggplot
#   p = spiderplot(hap, most.extreme.hap.idx,
#           left = 200, right = 200, max.haps = 1, palette = "RdBu", reverse = FALSE)
#   ggsave(p, file=paste0("bifurcation_chr", chr, ".", min.bp, "-", max.bp, "_ggplot.pdf"),
#       width=5, height=5)

write("Finished bifurcation plot", stderr())

setwd("../../..")

if (do.checkpointing) {
    if (center != 0) {
        save.image(paste0("after_rehh_chr", chr, ".",
            min.bp, ".", max.bp, ".END.Rdata"))
    } else {
        save.image(paste0("after_rehh_chr", chr, ".END.Rdata"))
    }
}
