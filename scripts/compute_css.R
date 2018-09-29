#!/usr/bin/env Rscript

# ========================================================================================
# --- Compute Composite Selection Score (CSS)
# ========================================================================================

# Examples:
#  Rscript scripts/compute_css.R

options(stringsAsFactors=FALSE)

library(ggplot2)

# ----------------------------------------------------------------------------------------
# --- Read in windowed Fst and XP-EHH files
# ----------------------------------------------------------------------------------------

stat.files = list()

stat.files[["fst"]] = list.files(path="results/",
        pattern="*.windowed.*fst",
        full.names=TRUE)

stat.files[["xpehh"]] = list.files(path="results/selscan/",
    pattern="xp-ehh.*.xpehh.out.norm.avg",
    full.names=TRUE)

read.stat = function (stat.file) {

    if (grepl("fst", stat.file)) {
        pop.str = strsplit(stat.file, "\\.")[[1]][2]
        pops = strsplit(pop.str, "_")[[1]]
        chr = gsub(".*chr([^\\.]+)\\..*", "\\1", stat.file)
    } else {
        pop.str = gsub(".*xp-ehh\\.(.*)\\..*\\.xpehh.*", "\\1", stat.file)
        pops = strsplit(pop.str, "-")[[1]]
        chr = strsplit(stat.file, "\\.")[[1]][3]
    }

    this.stat = read.table(stat.file, header=TRUE)

    if (nrow(this.stat) > 0) {
        this.stat$chr = chr
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
all.stat[["fst"]]   = do.call(rbind, lapply(stat.files[["fst"]],   read.stat))
all.stat[["xpehh"]] = do.call(rbind, lapply(stat.files[["xpehh"]], read.stat))

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

# Use absolute value of XP-EHH
all.stat[["xpehh"]]$statistic = abs(all.stat[["xpehh"]]$statistic)

# Remove NA values for PBS
xp.not.na.vals = which(is.na(all.stat[["xpehh"]]$statistic) == FALSE)
all.stat[["xpehh"]] = all.stat[["xpehh"]][xp.not.na.vals,]

# ----------------------------------------------------------------------------------------
# --- Compute empirical p-value
# ----------------------------------------------------------------------------------------

stat.ecdf.fst   = ecdf(all.stat[["fst"]]$statistic)
stat.ecdf.xpehh = ecdf(all.stat[["xpehh"]]$statistic)

all.stat[["fst"]]$p.val   = 1 - sapply(all.stat[["fst"]]$statistic,   stat.ecdf.fst)
all.stat[["xpehh"]]$p.val = 1 - sapply(all.stat[["xpehh"]]$statistic, stat.ecdf.xpehh)

# Convert p-values to z-scores
all.stat[["fst"]]$z   = qnorm(all.stat[["fst"]]$p.val)
all.stat[["xpehh"]]$z = qnorm(all.stat[["xpehh"]]$p.val)

# Set infinite z-scores to some real number
min.z.fst   = min(all.stat[["fst"]][abs(all.stat[["fst"]]$z)     != Inf,]$z, na.rm=TRUE)
min.z.xpehh = min(all.stat[["xpehh"]][abs(all.stat[["xpehh"]]$z) != Inf,]$z, na.rm=TRUE)

all.stat[["fst"]][abs(all.stat[["fst"]]$z)     == Inf,]$z = min.z.fst
all.stat[["xpehh"]][abs(all.stat[["xpehh"]]$z) == Inf,]$z = min.z.xpehh

# ----------------------------------------------------------------------------------------
# --- Merge Fst and XP-EHH data
# ----------------------------------------------------------------------------------------

both = merge(all.stat[["fst"]], all.stat[["xpehh"]],
    by=c("pos", "chr", "pop1", "pop2"), suffixes=c(".fst", ".xpehh"))

# ----------------------------------------------------------------------------------------
# --- Average z scores and covert back to p-value, which is then log transformed
# ----------------------------------------------------------------------------------------

both$z.both = rowMeans(cbind(both$z.fst, both$z.xpehh))
both$p.both = pnorm(-abs(both$z.both))

both$css = -1 * log(both$p.both, base=10)

# ----------------------------------------------------------------------------------------
# --- Write output
# ----------------------------------------------------------------------------------------

both.clean = both[,c(2,1,3,4,17)]

write.table(both.clean, file="results/css.all.txt",
    sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
