#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ========================================================================================
# --- Find H12 outliers
# ========================================================================================

library(xtable)
library(reshape2)

# ----------------------------------------------------------------------------------------
# --- Helper functions for formatting
# ----------------------------------------------------------------------------------------

simple.cap = function(x) {
    s = strsplit(x, " ")[[1]]
    lstr = paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
        sep="", collapse=" ")

    gsub("\\(i\\)", "(I)", gsub("\\(m\\)", "(M)", lstr))
}

make.italic <- function(x) {
    do.call(c, lapply(as.vector(x), function(y) {
        y = gsub("gambiae", "\\\\emph{gambiae}", gsub("coluzzii", "\\\\emph{coluzzii}", y))
    }))
}

# ----------------------------------------------------------------------------------------
# --- Bring in Ag1000G data for context
# ----------------------------------------------------------------------------------------

ag = read.table("data/h12.txt", header=TRUE, fill=TRUE)
ag.long = melt(ag, id.vars=names(ag)[1:3])

ag.cutoffs = lapply(names(ag)[-c(1:3)], function (site) {
    quantile(ag.long[ag.long$variable == site,]$value, 0.95, na.rm=TRUE)
})
names(ag.cutoffs) = names(ag)[-c(1:3)]
ag.cutoffs.df = cbind(site=names(ag.cutoffs), cutoff=ag.cutoffs)

ag.dict = list(
    AOM = 'Angola [coluzzii]',
    BFM = 'Burkina Faso [coluzzii]',
    BFS = 'Burkina Faso [gambiae]',
    CMS = 'Cameroon [gambiae]',
    GAS = 'Gabon [gambiae]',
    GNS = 'Guinea [gambiae]',
    GWA = 'Guinea-Bissau',
    UGS = 'Uganda [gambiae]'
)

# ----------------------------------------------------------------------------------------
# --- Read in H12 data
# ----------------------------------------------------------------------------------------

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

# Change Bugala labels
h12$pop[h12$pop == "BUGALAIS"] = "BUGALA (I)"
h12$pop[h12$pop == "BUGALAML"] = "BUGALA (M)"

h12$pop[h12$pop == "KAZZI"]   = "KAAZI"
h12$pop[h12$pop == "MITYANA"] = "WAMALA"

is.sites = c("BANDA", "BUGALA (I)", "BUKASA", "NSADZI", "SSERINYA")
ml.sites = c("BUGALA (M)", "BUWAMA", "KAAZI", "KIYINDI", "WAMALA")
h12$pop = factor(h12$pop, levels = c("", is.sites, ml.sites))

h12$chr = factor(h12$chr, levels = chrs)

cutoffs = lapply(c(is.sites, ml.sites), function (site) {
    quantile(h12[h12$pop == site,]$H12, 0.99)
})
names(cutoffs) = c(is.sites, ml.sites)
cutoffs.df = cbind(site=names(cutoffs), cutoff=cutoffs)

h12$is.island = FALSE
h12[h12$pop %in% is.sites,]$is.island = TRUE

# ----------------------------------------------------------------------------------------
# --- ID outliers
# ----------------------------------------------------------------------------------------

wins = do.call(rbind, lapply(chrs, function (chr) {
    chr.max = max(as.numeric(h12[h12$chr == chr,]$win.right))
    data.frame(
        chr=chr,
        win.start=seq(from=0, to=chr.max, by=100000),
        win.end=100000+seq(from=0, to=chr.max, by=100000)
    )
}))

assess.bin.for.is.ml = function (chr, bin.start, bin.end) {

    win = h12[h12$chr == chr & h12$win.center >= bin.start & h12$win.center <= bin.end,]

    # Remove Bugala mainland-like population
    win = win[win$pop != "BUGALA (M)",]

    win = merge(win, cutoffs.df, by.x="pop", by.y="site")

    win.high = win[win$H12 > win$cutoff,]

    if (nrow(win.high) == 0) {
        return(NA)
    }

    # Count of island sites with high H12
    is.high.ct = length(unique(win.high[win.high$is.island,]$pop))
    is.ct = length(unique(win[win$is.island,]$pop))

    # Count of mainland sites with high H12
    ml.high.ct = length(unique(win.high[!win.high$is.island,]$pop))
    ml.ct = length(unique(win[!win$is.island,]$pop))

    conting.tbl = matrix(c(is.high.ct, ml.high.ct, is.ct , ml.ct), nrow=2)
    chi = chisq.test(conting.tbl)

    # Any evidence of site-specific sweeps? (E.g., is only one site high?)
    high.pop.ct = length(unique(win.high$pop))
    high.pops = paste(unique(win.high$pop), collapse=';')

    res        = c(chr, bin.start, bin.end,
                   is.high.ct, is.ct, ml.high.ct, ml.ct,
                   chi$p.value,
                   high.pop.ct, high.pops)
    names(res) = c("chr", "bin.start", "bin.end",
                   "is.high.ct", "is.ct", "ml.high.ct", "ml.ct",
                   "chi.p",
                   "high.pop.ct", "high.pops")
    return(res)
}

isml.test.res = data.frame(do.call(rbind, lapply(1:nrow(wins), function (bin.idx) {
    assess.bin.for.is.ml(wins[bin.idx,1], wins[bin.idx,2], wins[bin.idx,3])
})))

save.image('h12.outlier.finding.Rdata')

isml.test.res$chi.p = as.numeric(isml.test.res$chi.p)

isml.test.res$bin.start = as.numeric(isml.test.res$bin.start)
isml.test.res$bin.end   = as.numeric(isml.test.res$bin.end)

is.specific = isml.test.res[which(isml.test.res$ml.high.ct <= 1 &
                                  isml.test.res$is.high.ct >= 4),]

ml.specific = isml.test.res[which(isml.test.res$ml.high.ct >= 3 &
                                  isml.test.res$is.high.ct <= 1),]

site.specific = isml.test.res[which(isml.test.res$high.pop.ct == 1),]

# ----------------------------------------------------------------------------------------
# --- Finalize island-specific results (merge contiguous regions)
# ----------------------------------------------------------------------------------------

is.specific.m = is.specific

is.specific.m$window.ct = 1

for (line in nrow(is.specific.m):2) {

    if (abs(is.specific.m[line-1,]$bin.end - is.specific.m[line,]$bin.start) < 200001 &
        is.specific.m[line-1,]$chr == is.specific.m[line,]$chr) {

        # Merge with above line
        is.specific.m[line-1,]$bin.end = is.specific.m[line,]$bin.end
        is.specific.m[line-1,]$window.ct = is.specific.m[line-1,]$window.ct + 1
        is.specific.m[line,]$window.ct = 0

        if (nchar(is.specific.m[line,]$high.pops) >
            nchar(is.specific.m[line-1,]$high.pops)) {

            is.specific.m[line-1,]$high.pops = is.specific.m[line,]$high.pops
            is.specific.m[line-1,]$high.pop.ct = is.specific.m[line,]$high.pop.ct
        }
    }
}

is.specific.m = is.specific.m[is.specific.m$window.ct != 0,]

is.specific.m$is.str = paste(is.specific.m$is.high.ct, "/", is.specific.m$is.ct)
is.specific.m$ml.str = paste(is.specific.m$ml.high.ct, "/", is.specific.m$ml.ct)

# Split sites into island and mainland

is.specific.m$high.pops = unlist(is.specific.m$high.pops)

is.specific.m$high.is.pops = do.call(c, lapply(strsplit(is.specific.m$high.pops, ";"),
    function (pops) { paste(is.sites[is.sites %in% pops], collapse=";") }))
is.specific.m$high.ml.pops = do.call(c, lapply(strsplit(is.specific.m$high.pops, ";"),
    function (pops) { paste(ml.sites[ml.sites %in% pops], collapse=";") }))

is.specific.m$high.is.pops[is.specific.m$high.is.pops == ""] = "None"
is.specific.m$high.ml.pops[is.specific.m$high.ml.pops == ""] = "None"

# --- Bring Ag1000G info into island-specific sweep exploration

is.specific.m$ag1kg = do.call(rbind, lapply(1:nrow(is.specific.m), function (row.idx) {

    this.roi = is.specific.m[row.idx,]
    ag.overlap = ag.long[ag.long$chrom == this.roi$chr &
                         ag.long$start >= this.roi$bin.start &
                         ag.long$stop <= this.roi$bin.end &
                         is.na(ag.long$value) == FALSE,]
    ag.overlap = merge(ag.overlap, ag.cutoffs.df, by.x="variable", "site")
    ag.overlap$is.outlier = ag.overlap$value >= ag.overlap$cutoff

    outlier.counts = aggregate(ag.overlap$is.outlier, by=list(ag.overlap$variable),
        FUN=sum)
    outlier.list = paste0(outlier.counts[outlier.counts$x > 0,]$Group.1, collapse=", ")
    if (nchar(outlier.list) == 0) {
        outlier.list = "None"
    }
    outlier.list
}))

# ----------------------------------------------------------------------------------------
# --- Finalize mainland-specific results (merge contiguous regions)
# ----------------------------------------------------------------------------------------

ml.specific.m = ml.specific

ml.specific.m$window.ct = 1

for (line in nrow(ml.specific.m):2) {

    if (abs(ml.specific.m[line-1,]$bin.end - ml.specific.m[line,]$bin.start) < 200001 &
        ml.specific.m[line-1,]$chr == ml.specific.m[line,]$chr) {

        # Merge with above line
        ml.specific.m[line-1,]$bin.end = ml.specific.m[line,]$bin.end
        ml.specific.m[line-1,]$window.ct = ml.specific.m[line-1,]$window.ct + 1
        ml.specific.m[line,]$window.ct = 0

        if (nchar(ml.specific.m[line,]$high.pops) >
            nchar(ml.specific.m[line-1,]$high.pops)) {

            ml.specific.m[line-1,]$high.pops = ml.specific.m[line,]$high.pops
            ml.specific.m[line-1,]$high.pop.ct = ml.specific.m[line,]$high.pop.ct
        }
    }
}

ml.specific.m = ml.specific.m[ml.specific.m$window.ct != 0,]

ml.specific.m$is.str = paste(ml.specific.m$is.high.ct, "/", ml.specific.m$is.ct)
ml.specific.m$ml.str = paste(ml.specific.m$ml.high.ct, "/", ml.specific.m$ml.ct)

# Split sites into island and mainland

ml.specific.m$high.pops = unlist(ml.specific.m$high.pops)

ml.specific.m$high.is.pops = do.call(c, lapply(strsplit(ml.specific.m$high.pops, ";"),
    function (pops) { paste(is.sites[is.sites %in% pops], collapse=";") }))
ml.specific.m$high.ml.pops = do.call(c, lapply(strsplit(ml.specific.m$high.pops, ";"),
    function (pops) { paste(ml.sites[ml.sites %in% pops], collapse=";") }))

ml.specific.m$high.is.pops[ml.specific.m$high.is.pops == ""] = "None"
ml.specific.m$high.ml.pops[ml.specific.m$high.ml.pops == ""] = "None"

# --- Bring Ag1000G info into mainland-specific sweep exploration

ml.specific.m$ag1kg = do.call(rbind, lapply(1:nrow(ml.specific.m), function (row.idx) {

    this.roi = ml.specific.m[row.idx,]
    ag.overlap = ag.long[ag.long$chrom == this.roi$chr &
                         ag.long$start >= this.roi$bin.start &
                         ag.long$stop <= this.roi$bin.end &
                         is.na(ag.long$value) == FALSE,]
    ag.overlap = merge(ag.overlap, ag.cutoffs.df, by.x="variable", "site")
    ag.overlap$is.outlier = ag.overlap$value >= ag.overlap$cutoff

    outlier.counts = aggregate(ag.overlap$is.outlier, by=list(ag.overlap$variable),
        FUN=sum)
    outlier.list = paste0(outlier.counts[outlier.counts$x > 0,]$Group.1, collapse=", ")
    if (nchar(outlier.list) == 0) {
        outlier.list = "None"
    }
    outlier.list
}))

# ----------------------------------------------------------------------------------------
# --- Do same for site-specific (within each island though)
# ----------------------------------------------------------------------------------------

site.specific.m = site.specific[order(site.specific$high.pops,
                                      site.specific$chr,
                                      site.specific$bin.start),]

site.specific.m$window.ct = 1

for (line in nrow(site.specific.m):2) {

    if (abs(site.specific.m[line-1,]$bin.end - site.specific.m[line,]$bin.start) < 200001 &
        site.specific.m[line-1,]$chr == site.specific.m[line,]$chr &
        site.specific.m[line-1,]$high.pops == site.specific.m[line,]$high.pops &
        sum(site.specific.m[line-1, 4:7] == site.specific.m[line, 4:7]) == 4) {

        # Merge with above line
        site.specific.m[line-1,]$bin.end = site.specific.m[line,]$bin.end
        site.specific.m[line-1,]$window.ct = site.specific.m[line-1,]$window.ct + 1
        site.specific.m[line,]$window.ct = 0

    }

}

site.specific.m = site.specific.m[site.specific.m$window.ct != 0,]

site.specific.ct = table(site.specific.m$high.pops)

# --- Bring Ag1000G info into site-specific sweep exploration

site.specific.m$ag1kg = do.call(rbind, lapply(1:nrow(site.specific.m), function (row.idx) {

    this.roi = site.specific.m[row.idx,]
    ag.overlap = ag.long[ag.long$chrom == this.roi$chr &
                         ag.long$start >= this.roi$bin.start &
                         ag.long$stop <= this.roi$bin.end &
                         is.na(ag.long$value) == FALSE,]
    ag.overlap = merge(ag.overlap, ag.cutoffs.df, by.x="variable", "site")
    ag.overlap$is.outlier = ag.overlap$value >= ag.overlap$cutoff

    outlier.counts = aggregate(ag.overlap$is.outlier, by=list(ag.overlap$variable),
        FUN=sum)
    outlier.list = paste0(outlier.counts[outlier.counts$x > 0,]$Group.1, collapse=", ")
    if (nchar(outlier.list) == 0) {
        outlier.list = "None"
    }
    outlier.list
}))

# ----------------------------------------------------------------------------------------
# --- Write tables - island-specific sweeps
# ----------------------------------------------------------------------------------------

is.toprint = is.specific.m[,c("chr", "bin.start", "bin.end",
    "is.str", "ml.str", "high.is.pops", "high.ml.pops", "ag1kg")]

is.toprint$bin.start = prettyNum(is.toprint$bin.start, big.mark=",", scientific=FALSE)
is.toprint$bin.end   = prettyNum(is.toprint$bin.end,   big.mark=",", scientific=FALSE)
is.toprint$high.is.pops = gsub(";", "; ", is.toprint$high.is.pops)
is.toprint$high.ml.pops = gsub(";", "; ", is.toprint$high.ml.pops)

# Capitalize
is.toprint[,6] = do.call(c, lapply(is.toprint[,6], simple.cap))
is.toprint[,7] = do.call(c, lapply(is.toprint[,7], simple.cap))

# Expand Ag1000G abbreviations
for (i in 1:length(ag.dict)) {
    is.toprint$ag1kg = gsub(names(ag.dict)[i], ag.dict[i], is.toprint$ag1kg)
}

# Italicize
is.toprint[,8] = make.italic(is.toprint[,8])

names(is.toprint) = c("", "", "",
                      "Island Sites with",
                      "Mainland Sites with",
                      "Outlier", "Outlier",
                      "Ag1000G Populations")
is.toprint = rbind(c("Chr.", "Region Start", "Region End",
                     "Putative Sweep", "Putative Sweep",
                     "Island Localities", "Mainland Localities",
                     "with Putative Sweep"), is.toprint)

display = c('s','s','s','s','s','s','s','s','s')
digits  = c( 0,  0,  0,  0,  0,  0,  0,  0,  0 )

caption = paste("Putative sweeps based on H12 statistic present on islands",
    "but rare or absent on LVB mainland.")

xt = xtable(is.toprint, display=display, digits=digits,
    caption=c(caption, "Putative sweeps present on islands, but rare or absent on LVB mainland"),
    label="table:is_sweeps")

align(xt) = c('l', 'r', 'r', 'r||', rep('c', 2), '|p{1.5in}|', 'l', '||p{1.5in}')

tex.out = "reports/island-specific-sweeps.tex"

sink(tex.out)

    cat("\\documentclass{article}",
        "\\usepackage{graphicx}",
        "\\usepackage{longtable}",
        "\\DeclareGraphicsExtensions{.pdf}",
        "\\usepackage[top=0.5in, bottom=0.5in, left=0.5in, right=0.5in]{geometry}",
        "\\usepackage{caption}",
        "\\captionsetup[table]{labelformat=empty}",
        "\\begin{document}", sep="\n")

    print.xtable(xt,
                include.rownames = FALSE,
                size="\\fontsize{9pt}{10pt}\\selectfont",
                tabular.environment = 'longtable', floating = FALSE,
                hline.after=sort(c(1, 1,
                    which(c(is.toprint[,1], NA) != c(NA, is.toprint[,1]))-1,
                    1:nrow(is.toprint),
                    nrow(is.toprint))),
                caption.placement = "top",
                sanitize.text.function = function(x){x})

    cat("\\end{document}", sep="\n")

sink()

# ----------------------------------------------------------------------------------------
# --- Write tables - mainland-specific sweeps
# ----------------------------------------------------------------------------------------

ml.toprint = ml.specific.m[,c("chr", "bin.start", "bin.end",
    "is.str", "ml.str", "high.is.pops", "high.ml.pops", "ag1kg")]

ml.toprint$bin.start = prettyNum(ml.toprint$bin.start, big.mark=",", scientific=FALSE)
ml.toprint$bin.end   = prettyNum(ml.toprint$bin.end,   big.mark=",", scientific=FALSE)
ml.toprint$high.is.pops = gsub(";", "; ", ml.toprint$high.is.pops)
ml.toprint$high.ml.pops = gsub(";", "; ", ml.toprint$high.ml.pops)

# Capitalize
ml.toprint[,6] = do.call(c, lapply(ml.toprint[,6], simple.cap))
ml.toprint[,7] = do.call(c, lapply(ml.toprint[,7], simple.cap))

# Expand Ag1000G abbreviations
for (i in 1:length(ag.dict)) {
    ml.toprint$ag1kg = gsub(names(ag.dict)[i], ag.dict[i], ml.toprint$ag1kg)
}

# Italicize
ml.toprint[,8] = make.italic(ml.toprint[,8])

names(ml.toprint) = c("", "", "",
                      "Island Sites with",
                      "Mainland Sites with",
                      "Outlier", "Outlier",
                      "Ag1000G Populations")
ml.toprint = rbind(c("Chr.", "Region Start", "Region End",
                     "Putative Sweep", "Putative Sweep",
                     "Island Localities", "Mainland Localities",
                     "with Putative Sweep"), ml.toprint)

display = c('s','s','s','s','s','s','s','s','s')
digits  = c( 0,  0,  0,  0,  0,  0,  0,  0,  0 )

caption = paste("Putative sweeps based on H12 statistic present on LVB mainland",
    "but rare or absent on islands.")

xt = xtable(ml.toprint, display=display, digits=digits,
    caption=c(caption, "Putative sweeps present on LVB mainland, but rare or absent on islands"),
    label="table:ml_sweeps")

align(xt) = c('l', 'r', 'r', 'r||', rep('c', 2), '|l', '|p{1.5in}', '||p{1.5in}')

tex.out = "reports/mainland-specific-sweeps.tex"

sink(tex.out)

    cat("\\documentclass{article}",
        "\\usepackage{graphicx}",
        "\\usepackage{longtable}",
        "\\DeclareGraphicsExtensions{.pdf}",
        "\\usepackage[top=0.5in, bottom=0.5in, left=0.5in, right=0.5in]{geometry}",
        "\\usepackage{caption}",
        "\\captionsetup[table]{labelformat=empty}",
        "\\begin{document}", sep="\n")

    print.xtable(xt,
                include.rownames = FALSE,
                size="\\fontsize{9pt}{10pt}\\selectfont",
                tabular.environment = 'longtable', floating = FALSE,
                hline.after=sort(c(1, 1,
                    which(c(ml.toprint[,1], NA) != c(NA, ml.toprint[,1]))-1,
                    1:nrow(ml.toprint),
                    nrow(ml.toprint))),
                caption.placement = "top",
                sanitize.text.function = function(x){x})

    cat("\\end{document}", sep="\n")

sink()

# ----------------------------------------------------------------------------------------
# --- Write tables - site-specific sweeps
# ----------------------------------------------------------------------------------------

site.sm = site.specific.m[,c("chr", "bin.start", "bin.end", "high.pops", "ag1kg")]

site.sm$region = paste(site.sm$bin.start / 1e6, "Mb")

site.agg = aggregate(site.sm$region, by=list(site.sm$high.pops, site.sm$chr),
    FUN=function(x) {paste(x, collapse="; ") })

site.agg.ag1kg = aggregate(site.sm$ag1kg, by=list(site.sm$high.pops, site.sm$chr),
    FUN=function(x) {
        y = x[x!="None"]
        paste(table(y), "also found in", names(table(y)), collapse="; ")
    })

site.agg$ag.info = site.agg.ag1kg$V1

site.agg = merge(site.agg, data.frame(t(site.specific.ct )), by.x="Group.1", by.y="Var2")
site.agg = site.agg[order(site.agg$Group.1, site.agg$Group.2),]

site.toprint = site.agg[,c("Group.1", "Freq", "Group.2", "x", "ag.info")]

later.rows = which(c(site.toprint$Group.1, NA) == c(NA, site.toprint$Group.1))
site.toprint[later.rows,]$Group.1 = ""
site.toprint[later.rows,]$Freq    = ""

site.toprint[site.toprint$ag.info == " also found in ",]$ag.info = ""

names(site.toprint) = c("Site", "Count", "Chr.", "Putative Sweeps",
    "Other sites\\textsuperscript{1}")

site.toprint$Site = as.vector(sapply(site.toprint$Site, simple.cap))

display = c('s','s','s','s','s','s')
digits  = c( 0,  0,  0,  0,  0,  0 )

caption = "Locality-specific (in LVB) putative sweeps based on H12 statistic."

xt = xtable(site.toprint, display=display, digits=digits,
    caption=c(caption, "Locality-specific putative sweeps"),
    label="table:site_sweeps")

align(xt) = c('l', 'l', 'l|', 'l', 'p{2.5in}', 'p{3.5in}')

tex.out = "reports/site-specific-sweeps.tex"

sink(tex.out)

    cat("\\documentclass{article}",
        "\\usepackage{graphicx}",
        "\\usepackage{longtable}",
        "\\DeclareGraphicsExtensions{.pdf}",
        "\\usepackage[top=0.5in, bottom=0.5in, left=0.5in, right=0.5in]{geometry}",
        "\\usepackage{caption}",
        "\\captionsetup[table]{labelformat=empty}",
        "\\begin{document}", sep="\n")

    print.xtable(xt,
                include.rownames = FALSE,
                size="\\fontsize{9pt}{10pt}\\selectfont",
                tabular.environment = 'longtable', floating = FALSE,
                hline.after=sort(c(0,
                    1:nrow(site.toprint),
                    which(site.toprint$Site != "") - 1,
                    nrow(site.toprint))),
                caption.placement = "top",
                sanitize.text.function = function(x){x})

    cat(paste0(
        "\\fontsize{9pt}{10pt}\\selectfont ",
        "\\textsuperscript{1} Ag1000G site codes: ",
        "AOM: Angola [\\emph{coluzzii}]; BFM: Burkina Faso [\\emph{coluzzii}]; ",
        "BFS: Burkina Faso [\\emph{gambiae}]; CMS: Cameroon [\\emph{gambiae}]; ",
        "GAS: Gabon [\\emph{gambiae}]; GNS: Guinea [\\emph{gambiae}]; ",
        "GWA: Guinea-Bissau; UGS: Uganda [\\emph{gambiae}]\n"
        ))
    cat("\\end{document}", sep="\n")

sink()

# ----------------------------------------------------------------------------------------
# --- Check on windows that contain known insecticide genes
# ----------------------------------------------------------------------------------------

genes = read.table("data/insecticide_genes_small.bed")

gene.info = lapply(1:nrow(genes), function (gene.idx) {

    res = isml.test.res[which(isml.test.res$chr == genes[gene.idx,]$V1 &
                              isml.test.res$bin.start < genes[gene.idx,]$V2 &
                              isml.test.res$bin.end > genes[gene.idx,]$V2),]

    if (nrow(res) == 0) {
        res = c(rep("", 3), 0, 5, 0, 4, 1, 0, "None")
    }

    c(genes[gene.idx,], res)
})

gene.info = data.frame(do.call(rbind, gene.info))

gene.info$high.pops = unlist(gene.info$high.pops)

gene.info$high.is.pops = do.call(c, lapply(strsplit(gene.info$high.pops, ";"),
    function (pops) { paste(is.sites[is.sites %in% pops], collapse=";") }))
gene.info$high.ml.pops = do.call(c, lapply(strsplit(gene.info$high.pops, ";"),
    function (pops) { paste(ml.sites[ml.sites %in% pops], collapse=";") }))

gene.info$high.is.pops[gene.info$high.is.pops == ""] = "None"
gene.info$high.ml.pops[gene.info$high.ml.pops == ""] = "None"

# --- Bring Ag1000G info into insecticide sweep exploration

gene.info$ag1kg = do.call(rbind, lapply(1:nrow(gene.info), function (row.idx) {

    this.roi = gene.info[row.idx,]
    ag.overlap = ag.long[ag.long$chrom == this.roi$chr &
                         ag.long$start >= this.roi$bin.start &
                         ag.long$stop <= this.roi$bin.end &
                         is.na(ag.long$value) == FALSE,]
    if (nrow(ag.overlap) > 0) {
        ag.overlap = merge(ag.overlap, ag.cutoffs.df, by.x="variable", "site")
        ag.overlap$is.outlier = ag.overlap$value >= ag.overlap$cutoff

        outlier.counts = aggregate(ag.overlap$is.outlier, by=list(ag.overlap$variable),
            FUN=sum)
        outlier.list = paste0(outlier.counts[outlier.counts$x > 0,]$Group.1,
            collapse=", ")
    } else {
        outlier.list = ""
    }

    if (nchar(outlier.list) == 0) {
        outlier.list = "None"
    }
    outlier.list
}))

# ----------------------------------------------------------------------------------------
# --- Write tables - known insecticide genes
# ----------------------------------------------------------------------------------------

gene.info$is.str = paste(gene.info$is.high.ct, "/", gene.info$is.ct)
gene.info$ml.str = paste(gene.info$ml.high.ct, "/", gene.info$ml.ct)

gene.info = gene.info[gene.info$chr != "",]

genes.toprint = gene.info[,c("V1", "V2", "V3",
    "is.str", "ml.str", "high.is.pops", "high.ml.pops", "ag1kg")]

genes.toprint$V2 = prettyNum(genes.toprint$V2, big.mark=",", scientific=FALSE)
genes.toprint$high.is.pops = gsub(";", "; ", genes.toprint$high.is.pops)
genes.toprint$high.ml.pops = gsub(";", "; ", genes.toprint$high.ml.pops)

# Capitalize
genes.toprint[,6] = do.call(c, lapply(genes.toprint[,6], simple.cap))
genes.toprint[,7] = do.call(c, lapply(genes.toprint[,7], simple.cap))

# Expand Ag1000G abbreviations
for (i in 1:length(ag.dict)) {
    genes.toprint$ag1kg = gsub(names(ag.dict)[i], ag.dict[i], genes.toprint$ag1kg)
}

# Italicize
genes.toprint[,8] = make.italic(genes.toprint[,8])

names(genes.toprint) = c("", "", "Insecticide",
                         "Island Sites with",
                         "Mainland Sites with",
                         "Outlier", "Outlier",
                         "Ag1000G Populations")
genes.toprint = rbind(c("Chr.", "Location", "Gene",
                        "Putative Sweep", "Putative Sweep",
                        "Island Localities", "Mainland Localities",
                        "with Putative Sweep"), genes.toprint)

display = c('s','s','s','s','s','s','s','s','s')
digits  = c( 0,  0,  0,  0,  0,  0,  0,  0,  0 )

caption = paste("Signatures of selective sweeps on known insecticide genes",
    "by site based on H12 statistic.")

xt = xtable(genes.toprint, display=display, digits=digits,
    caption=c(caption, "Insecticide genes' putative sweeps"),
    label="table:insecticide_sweeps")

align(xt) = c('l', 'r', 'r', 'c||', rep('c', 2), 'p{1in}', 'p{1in}', '|p{2in}')

tex.out = "reports/insecticide-sweeps.tex"

sink(tex.out)

    cat("\\documentclass{article}",
        "\\usepackage{graphicx}",
        "\\usepackage{longtable}",
        "\\DeclareGraphicsExtensions{.pdf}",
        "\\usepackage[top=0.5in, bottom=0.5in, left=0.5in, right=0.5in]{geometry}",
        "\\usepackage{caption}",
        "\\captionsetup[table]{labelformat=empty}",
        "\\begin{document}", sep="\n")

    print.xtable(xt,
                include.rownames = FALSE,
                size="\\fontsize{9pt}{10pt}\\selectfont",
                tabular.environment = 'longtable', floating = FALSE,
                hline.after=c(1, 1:nrow(genes.toprint),
                    nrow(genes.toprint)),
                caption.placement = "top",
                sanitize.text.function = function(x){x})

    cat("\\end{document}", sep="\n")

sink()

# ----------------------------------------------------------------------------------------
# --- Output a few stats about sweep sharing
# ----------------------------------------------------------------------------------------

sink("reports/sweep_sharing.txt")

cat("Proportion of island-specific sweeps also found in Uganda Ag1000G pop:\n")
sum(grepl("UGS", is.specific.m$ag1kg)) / length(is.specific.m$ag1kg)
cat("Proportion of island-specific sweeps also found in Gabon Ag1000G pop:\n")
sum(grepl("GAS", is.specific.m$ag1kg)) / length(is.specific.m$ag1kg)

cat("Proportion of mainland-specific sweeps also found in Uganda Ag1000G pop:\n")
sum(grepl("UGS", ml.specific.m$ag1kg)) / length(ml.specific.m$ag1kg)
cat("Proportion of mainland-specific sweeps also found in Gabon Ag1000G pop:\n")
sum(grepl("GAS", ml.specific.m$ag1kg)) / length(ml.specific.m$ag1kg)

sink()
