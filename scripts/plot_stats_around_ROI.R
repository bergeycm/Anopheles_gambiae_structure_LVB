#!/usr/bin/env Rscript

# module load R/3.3.0

options(stringsAsFactors=FALSE)

# Change settings to avoid use of exponential notation (e.g. 1.81e+08)
options("scipen"=100, "digits"=4)

library(ggplot2)
library(gtools)

# ========================================================================================
# --- Plot pop gen stats around putative sweeps
# ========================================================================================

args = commandArgs(trailingOnly=TRUE)
chr       = args[1]                # chr       = "X"
start     = as.numeric(args[2])    # start     = 9238942 - 1000000
end       = as.numeric(args[3])    # end       = 9238942 + 1000000
focal.pop = args[4]                # focal.pop = "NSADZI"
nickname  = args[5]                # nickname  = "X_9Mb"        # Short name for duplicate
                                                                # file naming

# ----------------------------------------------------------------------------------------
# --- Grab FST and XP-EHH (pairwise stats)
# ----------------------------------------------------------------------------------------

stat.files = list()

stat.files[["fst"]] = list.files(path="results/",
        pattern="*.windowed.*fst",
        full.names=TRUE)
stat.files[["fst"]] = stat.files[["fst"]][grepl("chr", stat.files[["fst"]])]
stat.files[["xpehh"]] = list.files(path="results/selscan/",
    pattern="xp-ehh.*.xpehh.out.norm.avg",
    full.names=TRUE)

this.chr.stat.files = lapply(stat.files, function(x) { x[grepl(chr, x)] })

read.stat = function (stat.file) {

    if (grepl("fst", stat.file)) {
        pop.str = strsplit(stat.file, "\\.")[[1]][2]
        pops = strsplit(pop.str, "_")[[1]]
    } else {
        pop.str = gsub(".*xp-ehh\\.(.*)\\..*\\.xpehh.*", "\\1", stat.file)
        pops = strsplit(pop.str, "-")[[1]]
    }

    this.stat = read.table(stat.file, header=TRUE)

    if (nrow(this.stat) > 0) {
        this.stat$pop1 = pops[1]
        this.stat$pop2 = pops[2]

        # Double the data.frame, switching pop1 and pop2
        num.cols = ncol(this.stat) - 2
        this.stat = rbind(this.stat,
                          data.frame(this.stat[c(1:num.cols)],
                                     pop1 = this.stat$pop2,
                                     pop2 = this.stat$pop1))
    }

    return(this.stat)
}

all.stat = list()
all.stat[["fst"]]   = do.call(rbind, lapply(this.chr.stat.files[["fst"]],   read.stat))
all.stat[["xpehh"]] = do.call(rbind, lapply(this.chr.stat.files[["xpehh"]], read.stat))

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

all.stat[["fst"]]$stat_name   = "FST"
all.stat[["xpehh"]]$stat_name = "XP-EHH"

# ----------------------------------------------------------------------------------------
# --- Grab H12 for each population
# ----------------------------------------------------------------------------------------

stat.files[["h12"]] = list.files(path="results/",
        pattern="*.phased\\..*\\.H12.out.txt",
        full.names=TRUE)

this.chr.stat.files = lapply(stat.files, function(x) { x[grepl(chr, x)] })

read.stat = function (stat.file) {

    pop = gsub(".*phased\\.(.*)\\.H12.out.txt", "\\1", stat.file)
    chr = gsub("results/*chr(.*)\\.pass.*", "\\1", stat.file)

    this.stat = read.table(stat.file, header=FALSE)

    if (nrow(this.stat) > 0) {
        this.stat$pop1 = pop
        this.stat$pop2 = NA
        this.stat$chr  = chr
    }

    names(this.stat) = c("win.center", "win.left", "win.right",
            "k.haplotypes", "hap.freq.spec",
            "ind", "H1.het", "H2", "H12", "H2_H1", "pop1", "pop2", "chr")

    return(this.stat)
}

all.stat[["h12"]] = do.call(rbind, lapply(this.chr.stat.files[["h12"]], read.stat))

# Standardize names
stat.col = which(names(all.stat[["h12"]]) == "H12")
pos.col  = which(names(all.stat[["h12"]]) == "win.left")
names(all.stat[["h12"]])[stat.col] = "statistic"
names(all.stat[["h12"]])[pos.col]  = "pos"

all.stat[["h12"]]$stat_name = "H12"

# ----------------------------------------------------------------------------------------
# --- Grab nucleotide diversity (pi) for each population
# ----------------------------------------------------------------------------------------

stat.files[["pi"]] = list.files(path="results/",
        pattern="*.windowed.pi",
        full.names=TRUE)

this.chr.stat.files = lapply(stat.files, function(x) { x[grepl(chr, x)] })

read.stat = function (stat.file) {

    pop = gsub(".*chr.*\\.(.*)\\.windowed.pi", "\\1", stat.file)
    chr = gsub(".*chr(.*)\\..*\\.windowed.pi", "\\1", stat.file)

    this.stat = read.table(stat.file, header=TRUE)

    if (nrow(this.stat) > 0) {
        this.stat$pop1 = pop
        this.stat$pop2 = NA
        this.stat$chr  = chr
    }

    return(this.stat)
}

all.stat[["pi"]] = do.call(rbind, lapply(this.chr.stat.files[["pi"]], read.stat))

# Standardize names
stat.col = which(names(all.stat[["pi"]]) == "PI")
pos.col  = which(names(all.stat[["pi"]]) == "BIN_START")
names(all.stat[["pi"]])[stat.col] = "statistic"
names(all.stat[["pi"]])[pos.col]  = "pos"

all.stat[["pi"]]$stat_name = "pi"

# ----------------------------------------------------------------------------------------
# --- Grab Tajima's D for each population
# ----------------------------------------------------------------------------------------

stat.files[["tajd"]] = list.files(path="results/",
        pattern="*.Tajima.D",
        full.names=TRUE)

this.chr.stat.files = lapply(stat.files, function(x) { x[grepl(chr, x)] })

read.stat = function (stat.file) {

    pop = gsub(".*chr.*\\.(.*)\\.Tajima.D", "\\1", stat.file)
    chr = gsub(".*chr(.*)\\..*\\.Tajima.D", "\\1", stat.file)

    this.stat = read.table(stat.file, header=TRUE)

    if (nrow(this.stat) > 0) {
        this.stat$pop1 = pop
        this.stat$pop2 = NA
        this.stat$chr  = chr
    }

    return(this.stat)
}

all.stat[["tajd"]] = do.call(rbind, lapply(this.chr.stat.files[["tajd"]], read.stat))

# Standardize names
stat.col = which(names(all.stat[["tajd"]]) == "TajimaD")
pos.col  = which(names(all.stat[["tajd"]]) == "BIN_START")
names(all.stat[["tajd"]])[stat.col] = "statistic"
names(all.stat[["tajd"]])[pos.col]  = "pos"

all.stat[["tajd"]]$stat_name = "Tajima's D"

# ----------------------------------------------------------------------------------------
# --- Grab SNP info (location, missingness, accessibility)
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# --- Grab gene info
# ----------------------------------------------------------------------------------------

# --- Write region to temporary BED file
tmp.bed = tempfile(pattern = "tmp_region_to_plot.", tmpdir = ".")
sink(tmp.bed)
cat(paste0(chr, "\t", start, "\t", end, "\n"))
sink()

genes = "data/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.9.gff3"
bedtools.cmd = paste("bedtools intersect -a", genes,
                     "-b", tmp.bed, " | awk '{ if ($3 == \"gene\") print $4,$5,$9 }'")

genes.df = read.table(text=system(bedtools.cmd, intern=TRUE))

names(genes.df) = c("gene.start", "gene.end", "gene.info")

genes.df$stat_name = "Genes"
genes.df$pop1 = focal.pop

genes.df$gene.id = gsub(".*ID=([^;]*);.*", "\\1", genes.df$gene.info)

genes.to.highlight = c()
if (chr == "2L") {
    genes.to.highlight = c("AGAP006549", "AGAP006550", "AGAP006551", "AGAP006553",
                           "AGAP006554", "AGAP006555", "AGAP006556")
} else if (chr == "X") {
    genes.to.highlight = c("AGAP000519", "AGAP000818")
} else if (chr == "2R") {
    genes.to.highlight = c("AGAP002862", "AGAP013128", "AGAP002863", "AGAP002865",
                           "AGAP002866", "AGAP002867", "AGAP002868", "AGAP002869",
                           "AGAP002870", "AGAP002894", "AGAP003343", "AGAP013511")
} else if (chr == "3R") {
    genes.to.highlight = c("AGAP009195", "AGAP009194", "AGAP009197", "AGAP009193",
                           "AGAP009192", "AGAP009191", "AGAP009196")
}

# devtools::install_github("wilkox/gggenes")

# ========================================================================================
# --- Combine stats into massive data.frame
# ========================================================================================

all.stat$h12$pop2  = as.character(all.stat$h12$pop2)
all.stat$pi$pop2   = as.character(all.stat$pi$pop2)
all.stat$tajd$pop2 = as.character(all.stat$tajd$pop2)

all.stat.df = do.call(smartbind, all.stat)

all.stat.df.sm = all.stat.df[all.stat.df$pos >= start & all.stat.df$pos <= end,]

all.stat.df.sm =  merge(all.stat.df.sm, genes.df, all=TRUE)

# ========================================================================================
# --- Do plot
# ========================================================================================

all.stat.df.sm$stat_name = factor(all.stat.df.sm$stat_name,
    levels = c("pi", "Tajima's D", "H12", "FST", "XP-EHH", "Genes"))

levels(all.stat.df.sm$stat_name) = c(expression(pi), expression("Tajima's"~D), "H12",
                                     expression("F[ST]"), "XP-EHH", "Genes")

# Change Bugala labels
all.stat.df.sm$pop2[all.stat.df.sm$pop2 == "BUGALAIS"] = "BUGALA (I)"
all.stat.df.sm$pop2[all.stat.df.sm$pop2 == "BUGALAML"] = "BUGALA (M)"

all.stat.df.sm$pop2[all.stat.df.sm$pop2 == "KAZZI"] = "KAAZI"
all.stat.df.sm$pop2[all.stat.df.sm$pop2 == "MITYANA"] = "WAMALA"

is.sites = c("BANDA", "BUGALA (I)", "BUKASA", "NSADZI", "SSERINYA")
ml.sites = c("BUGALA (M)", "BUWAMA", "KAAZI", "KIYINDI", "WAMALA")

simple.cap = function(x) {
    s = strsplit(x, " ")[[1]]
    lstr = paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
        sep="", collapse=" ")

    gsub("\\(i\\)", "(I)", gsub("\\(m\\)", "(M)", lstr))
}

all.stat.df.sm$pop1 = as.vector(sapply(all.stat.df.sm$pop1, simple.cap))
all.stat.df.sm$pop2 = as.vector(sapply(all.stat.df.sm$pop2, simple.cap))
is.sites = as.vector(sapply(is.sites, simple.cap))
ml.sites = as.vector(sapply(ml.sites, simple.cap))
focal.pop = simple.cap(focal.pop)

if (focal.pop == "Kazzi") {
    focal.pop = "Kaazi"
}
if (focal.pop == "Mityana") {
    focal.pop = "Wamala"
}
if (focal.pop == "Bugalaml") {
    focal.pop = "Bugala (M)"
}
if (focal.pop == "Bugalais") {
    focal.pop = "Bugala (I)"
}

all.stat.df.sm$pop2 = factor(all.stat.df.sm$pop2, levels=c(is.sites, ml.sites))

p = ggplot(subset(all.stat.df.sm, all.stat.df.sm$pop1 == focal.pop),
            aes(pos, statistic, col=pop2)) +
        geom_line() +
        geom_point(size=0.5, pch=1) +
        geom_rect(data=subset(all.stat.df.sm, all.stat.df.sm$stat_name == "Genes"),
            aes(xmin=gene.start, xmax=gene.end, ymin=0, ymax=1,
                fill=gene.id %in% genes.to.highlight), color=NA, size=0, alpha=1) +
        facet_grid(stat_name ~ pop1, scales="free_y", labeller=label_parsed) +
        theme_bw(base_size=15) +
#        expand_limits(x = 0, y = 0) +
        scale_x_continuous(labels=function(x) paste(x/1000000, "Mb"), expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        xlab(paste("Position along", chr)) + ylab("") +
        scale_color_discrete(name = "Comparison\nPopulation",
                             breaks = levels(all.stat.df.sm$pop2)) +
        scale_fill_manual(values=c("grey", "black"), guide=FALSE) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              #panel.border = element_blank(),
              axis.line = element_line(color = "black"),
              strip.background = element_rect(fill="white"),
              strip.text = element_text(color = 'black'),
              legend.position = "bottom",
              legend.background = element_rect())

out.file = paste0("reports/all_stats.",
    "chr", chr, "-", start, "-", end, ".", toupper(focal.pop), ".pdf")

ggsave(out.file, p, height=12, width=12)

# Save also as other file that doesn't use variable position info
out.file.nickname = paste0("reports/all_stats.", nickname, ".", toupper(focal.pop), ".pdf")

ggsave(out.file.nickname, p, height=12, width=12)
