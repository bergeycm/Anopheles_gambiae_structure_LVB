#!/usr/bin/env Rscript

# ========================================================================================
# --- Parse dadi output
# ========================================================================================

library(ggplot2)
library(scales)
library(xtable)

options(stringsAsFactors=FALSE)
options(scipen=999)

dadi.out = list.files(path="results/dadi_polarized", pattern="dadi.*.out.fix$",
    full.names=TRUE)

# ========================================================================================
# --- Parse and plot 1-pop results
# ========================================================================================

# Remove thinned for now
dadi.out = dadi.out[grep("thin", dadi.out, invert=TRUE)]

# And 2 population model results
dadi.out = dadi.out[grep("2pop", dadi.out, invert=TRUE)]

# And stuff that isn't island or Ag1000g
dadi.out = dadi.out[grepl("island", dadi.out) | grepl("ag1000g", dadi.out)]

# Results must already be fixed:
#    for file in `ls results/dadi_polarized/*out`; do
#        grep "^[id]" $file > $file.fix
#    done
dadi = do.call(rbind, lapply(dadi.out, function(x) {
    cbind(res.file = x, read.table(x, header=TRUE, sep="\t"))
}))

dadi$mode = gsub("dadi\\.([^\\.]*)\\.data", "\\1", gsub("chr3.", "", dadi$input_file))

# ----------------------------------------------------------------------------------------

# --- Remove results wth NaN values for FIM uncertainty estimates in best model

# rows.nan.ct = sapply(1:nrow(dadi), function (row.idx) {
#
#     this.row = dadi[row.idx,]
#     this.best.mod = this.row$max_ll_model
#     fim.col.idx = grep(paste0("FIM_", this.best.mod), names(dadi))
#     nan.ct = sum(is.na(this.row[fim.col.idx]))
# })
#
# dadi.orig = dadi
# dadi = dadi[rows.nan.ct == 0,]

# ----------------------------------------------------------------------------------------

# Find iteration with maximum log likelihood

dadi = do.call(rbind, lapply(unique(dadi$group), function (site) {
    this.site = dadi[dadi$group == site,]
    this.site$is.best = FALSE
    this.site[which(this.site$max_ll == max(this.site$max_ll)),]$is.best = TRUE
    this.site
}))

# Also find three-epoch model with maximum log likelihood

dadi = do.call(rbind, lapply(unique(dadi$group), function (site) {
    this.site = dadi[dadi$group == site,]
    this.site$is.best.three_epoch = FALSE
    ll.three = this.site$ll_three_epoch_uncert
    this.site[which(ll.three == max(ll.three)),]$is.best.three_epoch = TRUE
    this.site
}))

# ----------------------------------------------------------------------------------------

# --- Write file of best iterations, so we can grab the best optimization figures

best.iters = dadi[dadi$is.best,]$res.file

write.table(best.iters, file="results/dadi.best_iterations.txt",
    sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

# ----------------------------------------------------------------------------------------

dadi.is = dadi[dadi$mode=="island",]

# --- Fix site names
dadi.is$group[dadi.is$group == "BUGALAIS"] = "BUGALA (I)"
dadi.is$group[dadi.is$group == "BUGALAML"] = "BUGALA (M)"
dadi.is$group[dadi.is$group == "KAZZI"]    = "KAAZI"
dadi.is$group[dadi.is$group == "MITYANA"]  = "WAMALA"

islands  = c("BANDA", "BUGALA (I)", "BUKASA", "NSADZI", "SSERINYA")
mainland = c("BUGALA (M)", "BUWAMA", "KAAZI", "KIYINDI", "WAMALA")

dadi.is$group = factor(dadi.is$group, levels=c(mainland, islands))

dadi.is$island_mainland = "Mainland"
dadi.is[dadi.is$group %in% islands,]$island_mainland = "Island"

# save.image("parsing.dadi.Rdata")

AIC_cols = which(grepl("AIC_", names(dadi)))

min_AIC_model = gsub("AIC_", "", names(dadi.is[AIC_cols])[
    apply(dadi.is[AIC_cols], MARGIN=1, which.min)
])

# Check all results: do AIC and ll models pick best?
min_AIC_model = gsub("AIC_", "", names(dadi[AIC_cols])[
    apply(dadi[AIC_cols], MARGIN=1, which.min)
])

if (sum(min_AIC_model != dadi$max_ll_model) > 0) {
    quit("Log likelihood and AIC select different best models.", status=1)
}

# ----------------------------------------------------------------------------------------
# --- Scale to real units
# ----------------------------------------------------------------------------------------

# --- Set some constants:

# Per-chromosome arm size
L = list()
L["2L"]    = 49364325
L["2R"]    = 61545105
L["3L"]    = 41963435
L["3R"]    = 53200684

# And we need to subtract bits that were exlcuded for being heterozygous or inversions
inv = read.table("data/inversion_simple.bed")
het = read.table("data/heterochromatin.bed")
het$V4 = "het"
inv.het = rbind(inv, het)

names(inv.het) = c("chr", "start", "end", "name")
inv.het$length = inv.het$end - inv.het$start

mask.total = aggregate(inv.het$length, by=list(chr=inv.het$chr), FUN=sum)
mask.total.genome = sum(mask.total[mask.total$chr != "X",]$x)

L.fix = lapply(names(L), function(x) {
    L[[x]] = L[[x]] - mask.total[mask.total$chr == x, 2]
})
names(L.fix) = names(L)

L.genome = sum(do.call(c, L.fix))
L.chr3 = L.fix[["3L"]] + L.fix[["3R"]]

# Ag100G had: 30,001,805

# Mutation rate per year from Miles et al 2017
mu = 3.5e-9

# Generation time, 11 per year
g = 1/11

models = gsub("AIC_", "", names(dadi.is)[grepl("AIC_", names(dadi.is))])

# Do conversion to real units for each model:
#    snm         : theta
#    two_epoch   : theta nu T
#    growth      : theta nu T
#    bottle      : theta nuB nuF T
#    three_epoch : theta nuB nuF TB TF

for (m in models) {

    # --- Compute effective population size - ancestral
    # theta = 4 * mu * L * Nref
    dadi.is[[paste0(m, "_Na")]]     =
        dadi.is[[paste0("theta_", m)]]         / (4 * mu * L.chr3)
    dadi.is[[paste0(m, "_GIM_Na")]] =
        dadi.is[[paste0("GIM_", m, "_theta")]] / (4 * mu * L.chr3)
    dadi.is[[paste0(m, "_FIM_Na")]] =
        dadi.is[[paste0("FIM_", m, "_theta")]] / (4 * mu * L.chr3)

    # --- Convert T to real time in years
    # time in years = 2 * Nref * T * g
    for (T.type in c("T", "TB", "TF")) {

        if (paste0(m, "_", T.type) %in% names(dadi.is)) {

            dadi.is[[paste0(m, "_", T.type, ".real")]]     =
                2 * dadi.is[[paste0(m, "_Na")]] * dadi.is[[paste0(m, "_", T.type)]] * g
            dadi.is[[paste0(m, "_GIM_", T.type, ".real")]] =
                2 * dadi.is[[paste0(m, "_Na")]] *
                dadi.is[[paste0("GIM_", m, "_", T.type)]] * g
            dadi.is[[paste0(m, "_FIM_", T.type, ".real")]] =
                2 * dadi.is[[paste0(m, "_Na")]] *
                dadi.is[[paste0("FIM_", m, "_", T.type)]] * g
        }
    }
}

# ----------------------------------------------------------------------------------------
# --- Plot population size history - scaled to real units
# ----------------------------------------------------------------------------------------

min.x = 1000
max.x = 10e5
cut.on.left =
    min(dadi.is$three_epoch_uncert_TF.real) < min.x
cut.on.right =
    max(dadi.is$three_epoch_uncert_TF.real + dadi.is$three_epoch_uncert_TB.real) > max.x
if (cut.on.left | cut.on.right) {
    warning("Plotting region cut.")
}

p = ggplot(subset(dadi.is, dadi.is$is.best.three_epoch & !is.na(dadi.is$group)),
        aes(color=island_mainland)) +
    geom_segment(aes(
        x = min.x,
        xend = three_epoch_uncert_TF.real,
        y = three_epoch_uncert_nuF * three_epoch_uncert_Na,
        yend = three_epoch_uncert_nuF * three_epoch_uncert_Na)) +
    geom_segment(aes(
        x = three_epoch_uncert_TF.real,
        xend = three_epoch_uncert_TF.real,
        y = three_epoch_uncert_nuF * three_epoch_uncert_Na,
        yend = three_epoch_uncert_nuB * three_epoch_uncert_Na)) +
    geom_segment(aes(
        x = three_epoch_uncert_TF.real,
        xend = three_epoch_uncert_TF.real + three_epoch_uncert_TB.real,
        y = three_epoch_uncert_nuB * three_epoch_uncert_Na,
        yend = three_epoch_uncert_nuB * three_epoch_uncert_Na)) +
    geom_segment(aes(
        x = three_epoch_uncert_TF.real + three_epoch_uncert_TB.real,
        xend = three_epoch_uncert_TF.real + three_epoch_uncert_TB.real,
        y = three_epoch_uncert_nuB * three_epoch_uncert_Na,
        yend = three_epoch_uncert_Na)) +
    geom_segment(aes(
        x = three_epoch_uncert_TF.real + three_epoch_uncert_TB.real,
        xend = max.x,
        y = three_epoch_uncert_Na,
        yend = three_epoch_uncert_Na)) +
    scale_x_continuous(name   = "Time (years before present)",
                       breaks = c(10e2, 10e3, 10e4, 10e5, 10e6),
                       limits = c(min.x, max.x),
                       labels = trans_format('log10',math_format(10^.x))) +
    scale_y_continuous(name   = expression("N"["e"]),
                       breaks = c(10e4, 10e5, 10e6, 10e7, 10e8),
                       limits = c(10e4,10e8),
                       expand = c(0,0),
                       labels = trans_format('log10',math_format(10^.x))) +
    annotation_logticks(short = unit(.5,"mm"),
                        mid   = unit(2,"mm"),
                        long  = unit(3,"mm"),
                        scaled = FALSE, sides='bl', color='grey') +
    coord_trans(y = "log10", x = "log10") +
    facet_wrap( ~ group, nrow=2) +
    theme_bw() +
    theme(panel.grid.major = element_line(linetype=3),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank()) +
    scale_color_discrete(guide=FALSE)

ggsave(p, file="reports/dadi_three_epoch_pop_history.pdf", height=6, width=12)

# ----------------------------------------------------------------------------------------
# --- Write table of key results
# ----------------------------------------------------------------------------------------

simpleCap = function (x) {
    s = strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
        sep="", collapse=" ")
}

signif.3 = function(num) {
    format(c(signif(num, 3)), big.mark=",", trim=TRUE)
}

dadi.is.out = do.call(rbind, lapply(c(islands, mainland), function (site) {

    site = as.character(site)
    this.site.res = list()
    this.site.res[["site"]] = gsub("\\(m", "(M", gsub("\\(i", "(I", simpleCap(site)))

    # Blank for the table hack
    this.site.res[["site_GIM"]] = ""
    this.site.res[["site_FIM"]] = ""

    this.site = dadi.is[which(dadi.is$group == site & dadi.is$is.best),]
    mod = this.site$max_ll_model
    this.site.res[["model"]] = simpleCap(gsub("_", " ", gsub("_uncert", "", mod)))
    # Blank for the table hack
    this.site.res[["model_GIM"]] = ""
    this.site.res[["model_FIM"]] = ""

    this.site.res[["ll"]] = signif.3(this.site$max_ll)
    this.site.res[["ll_GIM"]] = ""
    this.site.res[["ll_FIM"]] = ""

    this.site.mod = this.site[grepl(mod, names(this.site))]

    # --- Na
    Na.vals = this.site.mod[grepl("Na$", names(this.site.mod))]
    this.site.res[["Na"]] = signif.3(Na.vals[1])
    Na.gim.tmp = as.numeric(Na.vals[1]) + c(-1,1) * 1.96 * as.numeric(Na.vals[2])
    this.site.res[["Na_GIM"]] = paste0("(", signif.3(Na.gim.tmp[1]), ", ",
                                            signif.3(Na.gim.tmp[2]), ")")
    Na.fim.tmp = as.numeric(Na.vals[1]) + c(-1,1) * 1.96 * as.numeric(Na.vals[3])
    this.site.res[["Na_FIM"]] = paste0("(", signif.3(Na.fim.tmp[1]), ", ",
                                            signif.3(Na.fim.tmp[2]), ")")

    # --- NB
    NB.vals = NA
    this.site.res[["NB"]] = "-"
    this.site.res[["NB_GIM"]] = ""
    this.site.res[["NB_FIM"]] = ""
    if (sum(grepl("nuB", names(this.site.mod)))) {
        NB.vals = this.site.mod[grepl("nuB$", names(this.site.mod))]

        this.site.res[["NB"]] = signif.3(NB.vals[1])
        NB.gim.tmp = as.numeric(NB.vals[1]) + c(-1,1) * 1.96 * as.numeric(NB.vals[2])
        this.site.res[["NB_GIM"]] = paste0("(", signif.3(NB.gim.tmp[1]), ", ",
                                                signif.3(NB.gim.tmp[2]), ")")
        NB.fim.tmp = as.numeric(NB.vals[1]) + c(-1,1) * 1.96 * as.numeric(NB.vals[3])
        this.site.res[["NB_FIM"]] = paste0("(", signif.3(NB.fim.tmp[1]), ", ",
                                                signif.3(NB.fim.tmp[2]), ")")
    }

    # --- NF
    NF.vals = NA
    this.site.res[["NF"]] = "-"
    this.site.res[["NF_GIM"]] = ""
    this.site.res[["NF_FIM"]] = ""
    if (sum(grepl("nu", names(this.site.mod)))) {
        if (sum(grepl("nuF", names(this.site.mod)))) {
            NF.vals = this.site.mod[grepl("nuF$", names(this.site.mod))]
        } else {
            NF.vals = this.site.mod[grepl("nu$", names(this.site.mod))]
        }

        this.site.res[["NF"]] = signif.3(NF.vals[1])
        NF.gim.tmp = as.numeric(NF.vals[1]) + c(-1,1) * 1.96 * as.numeric(NF.vals[2])
        this.site.res[["NF_GIM"]] = paste0("(", signif.3(NF.gim.tmp[1]), ", ",
                                                signif.3(NF.gim.tmp[2]), ")")
        NF.fim.tmp = as.numeric(NF.vals[1]) + c(-1,1) * 1.96 * as.numeric(NF.vals[3])
        this.site.res[["NF_FIM"]] = paste0("(", signif.3(NF.fim.tmp[1]), ", ",
                                                signif.3(NF.fim.tmp[2]), ")")
    }

    # --- TB
    TB.vals = NA
    this.site.res[["TB"]] = "-"
    this.site.res[["TB_GIM"]] = ""
    this.site.res[["TB_FIM"]] = ""
    if (sum(grepl("_T", names(this.site.mod)))) {
        if (sum(grepl("TB.real", names(this.site.mod)))) {
            TB.vals = this.site.mod[grepl("TB.real$", names(this.site.mod))]

            this.site.res[["TB"]] = signif.3(TB.vals[1])
            TB.gim.tmp = as.numeric(TB.vals[1]) + c(-1,1) * 1.96 * as.numeric(TB.vals[2])
            this.site.res[["TB_GIM"]] = paste0("(", signif.3(TB.gim.tmp[1]), ", ",
                                                    signif.3(TB.gim.tmp[2]), ")")
            TB.fim.tmp = as.numeric(TB.vals[1]) + c(-1,1) * 1.96 * as.numeric(TB.vals[3])
            this.site.res[["TB_FIM"]] = paste0("(", signif.3(TB.fim.tmp[1]), ", ",
                                                    signif.3(TB.fim.tmp[2]), ")")
        }
    }

    # --- TF
    TF.vals = NA
    this.site.res[["TF"]] = "-"
    this.site.res[["TF_GIM"]] = ""
    this.site.res[["TF_FIM"]] = ""
    if (sum(grepl("_T", names(this.site.mod)))) {
        if (sum(grepl("TF.real", names(this.site.mod)))) {
            TF.vals = this.site.mod[grepl("TF.real$", names(this.site.mod))]
        } else {
            TF.vals = this.site.mod[grepl("T.real$", names(this.site.mod))]
        }

        this.site.res[["TF"]] = signif.3(TF.vals[1])
        TF.gim.tmp = as.numeric(TF.vals[1]) + c(-1,1) * 1.96 * as.numeric(TF.vals[2])
        this.site.res[["TF_GIM"]] = paste0("(", signif.3(TF.gim.tmp[1]), ", ",
                                                signif.3(TF.gim.tmp[2]), ")")
        TF.fim.tmp = as.numeric(TF.vals[1]) + c(-1,1) * 1.96 * as.numeric(TF.vals[3])
        this.site.res[["TF_FIM"]] = paste0("(", signif.3(TF.fim.tmp[1]), ", ",
                                                signif.3(TF.fim.tmp[2]), ")")
    }

    return(this.site.res)
}))

dadi.is.out.main = dadi.is.out[,seq(from=1, to=ncol(dadi.is.out), by=3)]
dadi.is.out.GIM  = dadi.is.out[,seq(from=2, to=ncol(dadi.is.out), by=3)]
dadi.is.out.FIM  = dadi.is.out[,seq(from=3, to=ncol(dadi.is.out), by=3)]

# Now using just FIM
dadi.is.out.sort = do.call(rbind, lapply(1:nrow(dadi.is.out.main), function (line.num) {
    return(rbind(dadi.is.out.main[line.num,],
                 dadi.is.out.FIM [line.num,]))
}))

dadi.is.out.sort = data.frame(dadi.is.out.sort)

dadi.is.out.sort$model = gsub("Snm", "Standard Neutral", dadi.is.out.sort$model)

# ----------------------------------------------------------------------------------------
# --- Make xtable
# ----------------------------------------------------------------------------------------

names(dadi.is.out.sort) = c("Site", "Best Model", "LnL",
                           "Ancestral Size",
                           "Intermediate Size", "Final Size",
                           "Intermediate Duration",
                           "Time Since Change")

display = c('s','s','s','s','s','s','s','s','s')
digits  = c( 0,  0,  0,  0,  0,  0,  0,  0,  0 )

caption.1pop = paste0("Results of single population demographic inference with ",
    "$\\delta$a$\\delta$i for all sampling sites considered individually. ",
    "Numbers in parentheses are bounds of 95\\% confidence interval computed using ",
    "Fisher information matrix and 100 bootstrap replicates of 1 Mb from the dataset.")

xt = xtable(dadi.is.out.sort, display=display, digits=digits,
    caption=c(caption.1pop, "$\\delta$a$\\delta$i single population results."),
    label="table:dadi_1pop")

align(xt) = c('l', 'p{1in}||', 'l', rep('r', 6))

tex.out = "reports/dadi.island.out.tex"

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
                hline.after=c(0, 0,
                    seq(from=2, to=nrow(dadi.is.out.sort), by=2),
                    length(islands) * 2,
                    length(c(islands, mainland)) * 2),
                caption.placement = "top",
                sanitize.text.function = function(x){x})

    cat("\\end{document}", sep="\n")

sink()

# ----------------------------------------------------------------------------------------
# --- Test for significant difference in...
# ----------------------------------------------------------------------------------------

sink("reports/dadi_island_is-ml_differences.txt")

is.three_epoch.bool = dadi.is$is.best.three_epoch & dadi.is$island_mainland == "Island"
ml.three_epoch.bool = dadi.is$is.best.three_epoch & dadi.is$island_mainland == "Mainland"

# --- ...ancestral population size
cat("Ancestral population size (three epoch):\n")
paste("Island median:",   median(dadi.is[is.three_epoch.bool,]$three_epoch_uncert_Na))
paste("Mainland median:", median(dadi.is[ml.three_epoch.bool,]$three_epoch_uncert_Na))

wilcox.test(dadi.is[dadi.is$is.best.three_epoch,]$three_epoch_uncert_Na ~
            dadi.is[dadi.is$is.best.three_epoch,]$island_mainland)

# --- ...final population size
cat("Final population size (three epoch):\n")
tmp.final.is = dadi.is[is.three_epoch.bool,]$three_epoch_uncert_Na *
               dadi.is[is.three_epoch.bool,]$three_epoch_uncert_nuF
tmp.final.ml = dadi.is[ml.three_epoch.bool,]$three_epoch_uncert_Na *
               dadi.is[ml.three_epoch.bool,]$three_epoch_uncert_nuF
paste("Island median:",   median(tmp.final.is))
paste("Mainland median:", median(tmp.final.ml))

wilcox.test(dadi.is[dadi.is$is.best.three_epoch,]$three_epoch_uncert_Na *
            dadi.is[dadi.is$is.best.three_epoch,]$three_epoch_uncert_nuF ~
            dadi.is[dadi.is$is.best.three_epoch,]$island_mainland)

# --- ...time since size change
cat("Time since size change (three epoch):\n")
paste("Island median:",
                        median(dadi.is[is.three_epoch.bool,]$three_epoch_uncert_TF.real))
paste("Mainland median:",
                        median(dadi.is[ml.three_epoch.bool,]$three_epoch_uncert_TF.real))

wilcox.test(dadi.is[dadi.is$is.best.three_epoch,]$three_epoch_uncert_TF.real ~
            dadi.is[dadi.is$is.best.three_epoch,]$island_mainland)

sink()

# save.image("end_of_1pop.Rdata")

# ========================================================================================
# --- Parse and plot 2-pop results
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Plot 2-population model results
# ----------------------------------------------------------------------------------------

dadi.2pop.out = list.files(path="results/dadi_polarized",
                           pattern="dadi.island.2pop.*.out.fix",
                           full.names=TRUE)

# Remove empty files
dadi.2pop.out = dadi.2pop.out[
    sapply(dadi.2pop.out, function (file) { file.info(file)$size != 0 })
]

col.names.2pop = c("input_file", "pop1", "pop2", "sample_size1", "sample_size2",
                   "ll", "AIC", "theta",
                   "s", "nu1", "nu2",
                   "T", "m12", "m21", "pmisid",
                   "GIM_s", "GIM_nu1", "GIM_nu2",
                   "GIM_T", "GIM_m12", "GIM_m21", "GIM_pmisid",
                   "GIM_theta",
                   "FIM_s", "FIM_nu1", "FIM_nu2",
                   "FIM_T", "FIM_m12", "FIM_m21", "FIM_pmisid",
                   "FIM__theta", "max_iter")

dadi.2pop = do.call(rbind, lapply(dadi.2pop.out, function(x) {

    orig = read.table(x, header=TRUE, sep="\t")

    # Add NA migration columns if they don't exist
    if (! "m12" %in% names(orig)) {
        for (m.col in col.names.2pop[grepl("m[12]", col.names.2pop)]) {
            orig[[m.col]] = NA
        }
    }
    cbind(res.file = x, orig[,col.names.2pop])
}))

dadi.2pop$mode = gsub("chr3.", "",
    gsub("dadi\\.(.*)\\.data", "\\1", dadi.2pop$input_file))

dadi.2pop$both.pops = paste(dadi.2pop$pop1, dadi.2pop$pop2, sep="-")

# ----------------------------------------------------------------------------------------

# --- Remove results wth NaN values for FIM uncertainty estimates

dadi.2pop = do.call(rbind,
    lapply(unique(dadi.2pop$both.pops), function (pair) {
        this.pair = dadi.2pop[dadi.2pop$both.pops == pair,]
        pops = strsplit(pair, split="-")[[1]]
        nan.ct = rowSums(dadi.2pop[dadi.2pop$pop1 == pops[1] &
                                   dadi.2pop$pop2 == pops[2],
                                   grepl("FIM", names(dadi.2pop))] == "NaN", na.rm=TRUE)
        this.pair[nan.ct == 0,]
    })
)

# Migration rejected before filtering for NaNs
table(is.na(dadi.2pop$m12))

# ----------------------------------------------------------------------------------------

# --- Find iteration with maximum log likelihood

dadi.2pop = do.call(rbind,
    lapply(unique(dadi.2pop$both.pops), function (pair) {
        this.pair = dadi.2pop[dadi.2pop$both.pops == pair,]
        this.pair$is.best = FALSE
        this.pair[which(this.pair$ll == max(this.pair$ll)),]$is.best = TRUE
        this.pair
    })
)

# Migration rejected after filtering for NaNs
table(is.na(dadi.2pop$m12))

# ----------------------------------------------------------------------------------------

# --- Write file of best iterations, so we can grab the best optimization figures

best.iters.2pop = dadi.2pop[dadi.2pop$is.best,]$res.file

write.table(best.iters.2pop, file="results/dadi.best_iterations.2pop.txt",
    sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

# ----------------------------------------------------------------------------------------

# --- Fix site names
fix.names = function (names) {

    sapply(names, function (x) {
        x = gsub("BUGALAIS", "BUGALA (I)", x)
        x = gsub("BUGALAML", "BUGALA (M)", x)
        x = gsub("KAZZI",    "KAAZI",      x)
        x = gsub("MITYANA",  "WAMALA",     x)
    })
}

dadi.2pop$pop1      = fix.names(dadi.2pop$pop1)
dadi.2pop$pop2      = fix.names(dadi.2pop$pop2)
dadi.2pop$both.pops = fix.names(dadi.2pop$both.pops)

islands  = c("BANDA", "BUGALA (I)", "BUKASA", "NSADZI", "SSERINYA")
mainland = c("BUGALA (M)", "BUWAMA", "KAAZI", "KIYINDI", "WAMALA")

# ----------------------------------------------------------------------------------------

# IM:
# https://github.com/paulirish/dadi/blob/master/dadi/Demographics2D.py
# Isolation-with-migration model with exponential pop growth.
#     s:    *Fraction* of pop 1 after split. (Pop 2 has size 1-s.)
#     nu1:  Final size of pop 1.
#     nu2:  Final size of pop 2.
#     T:    Time in the past of split (in units of 2*Na generations)
#     m12:  Migration from pop 2 to pop 1 (2*Na*m12)
#     m21:  Migration from pop 1 to pop 2

dadi.2pop$Na     = dadi.2pop$theta      / (4 * mu * L.chr3)
dadi.2pop$GIM_Na = dadi.2pop$GIM_theta  / (4 * mu * L.chr3)
dadi.2pop$FIM_Na = dadi.2pop$FIM__theta / (4 * mu * L.chr3)   # Typo in dadi script: "__"

# Final size of pop1
dadi.2pop$nu1.real = dadi.2pop$Na * dadi.2pop$nu1
# Final size of pop2
dadi.2pop$nu2.real = dadi.2pop$Na * dadi.2pop$nu2

# Convert T to real time in years
dadi.2pop$T.real     = 2 * dadi.2pop$Na * dadi.2pop$T     * g
dadi.2pop$GIM_T.real = 2 * dadi.2pop$Na * dadi.2pop$GIM_T * g
dadi.2pop$FIM_T.real = 2 * dadi.2pop$Na * dadi.2pop$FIM_T * g

# And get proportion of migrants
# "fraction of individuals in each generation in population i who are new migrants from
#  population j." Multiply by population size to get number of individuals each generation
#  that are migrating from population j to population i."
dadi.2pop$m12.real = dadi.2pop$m12 / (2 * dadi.2pop$Na)
dadi.2pop$m21.real = dadi.2pop$m21 / (2 * dadi.2pop$Na)

dadi.2pop$GIM_m12.real = dadi.2pop$GIM_m12 / (2 * dadi.2pop$Na)
dadi.2pop$GIM_m21.real = dadi.2pop$GIM_m21 / (2 * dadi.2pop$Na)

dadi.2pop$FIM_m12.real = dadi.2pop$FIM_m12 / (2 * dadi.2pop$Na)
dadi.2pop$FIM_m21.real = dadi.2pop$FIM_m21 / (2 * dadi.2pop$Na)

# And get effective number of individuals migrating from,
# e.g., population i to j each generation
dadi.2pop$m12.real.inds = dadi.2pop$m12.real * dadi.2pop$nu1.real
dadi.2pop$m21.real.inds = dadi.2pop$m21.real * dadi.2pop$nu2.real

dadi.2pop$GIM_m12.real.inds = dadi.2pop$GIM_m12.real * dadi.2pop$nu1.real
dadi.2pop$GIM_m21.real.inds = dadi.2pop$GIM_m21.real * dadi.2pop$nu2.real

dadi.2pop$FIM_m12.real.inds = dadi.2pop$FIM_m12.real * dadi.2pop$nu1.real
dadi.2pop$FIM_m21.real.inds = dadi.2pop$FIM_m21.real * dadi.2pop$nu2.real

# ----------------------------------------------------------------------------------------

# Figure out which comparison class we're looking at
island.island.boolean =     dadi.2pop$pop1 %in% islands &
                            dadi.2pop$pop2 %in% islands
mainland.mainland.boolean = dadi.2pop$pop1 %in% mainland &
                            dadi.2pop$pop2 %in% mainland
dadi.2pop$comparison.type = "island.mainland"
dadi.2pop[island.island.boolean,]$comparison.type     = "island.island"
dadi.2pop[mainland.mainland.boolean,]$comparison.type = "mainland.mainland"

# ----------------------------------------------------------------------------------------

# --- See how migration differs between comparison classes

sink("reports/dadi_2pop_is-ml_differences.txt")

# --- Median migration (setting no-migration models to zero)
cat("Median migration:\n")

m.real.vals = c(dadi.2pop[dadi.2pop$is.best,]$m12.real.inds,
                dadi.2pop[dadi.2pop$is.best,]$m21.real.inds)
m.real.vals.withzero = m.real.vals
m.real.vals.withzero[is.na(m.real.vals.withzero)] = 0

m.real.inds.median =  aggregate(m.real.vals.withzero,
    by=list(rep(dadi.2pop[dadi.2pop$is.best,]$comparison.type, 2)),
    function (x) median(x, na.rm=TRUE))
m.real.inds.median

# --- Median migration excluding zeros
cat("Median migration excluding zeros:\n")

m.real.inds.median.gt0 = aggregate(m.real.vals,
    by=list(rep(dadi.2pop[dadi.2pop$is.best,]$comparison.type,2)),
    function (x) median(x, na.rm=TRUE))
m.real.inds.median.gt0

# --- Median migration excluding zeros
cat("Proportion of no migration comparisons:\n")

m.real.inds.none.prop =  aggregate(m.real.vals,
    by=list(rep(dadi.2pop[dadi.2pop$is.best,]$comparison.type,2)),
    function (x) sum(is.na(x)) / length(x))
m.real.inds.none.prop

# --- Median time since split excluding zero-migration models
cat("Median time since split excluding zeros:\n")

best.and.with.mig = dadi.2pop$is.best & !is.na(dadi.2pop$m12)

t.real.inds.median.gt0 = aggregate(dadi.2pop[best.and.with.mig,]$T.real,
    by=list(dadi.2pop[best.and.with.mig,]$comparison.type),
    function (x) median(x, na.rm=TRUE))
t.real.inds.median.gt0

# --- Ancestral population size excluding zero-migration models
cat("Ancestral population size excluding zeros:\n")

na.median.gt0 = aggregate(dadi.2pop[best.and.with.mig,]$Na,
    by=list(dadi.2pop[best.and.with.mig,]$comparison.type),
    function (x) median(x, na.rm=TRUE))
na.median.gt0

# --- Final population size excluding zero-migration models
cat("Final population size excluding zeros:\n")

island.nu =   c(dadi.2pop[best.and.with.mig & dadi.2pop$pop1 %in% islands ,]$nu1.real,
                dadi.2pop[best.and.with.mig & dadi.2pop$pop2 %in% islands ,]$nu2.real)
mainland.nu = c(dadi.2pop[best.and.with.mig & dadi.2pop$pop1 %in% mainland,]$nu1.real,
                dadi.2pop[best.and.with.mig & dadi.2pop$pop2 %in% mainland,]$nu2.real)

cat(paste("- Island final pop:",   median(island.nu)), "\n")
cat(paste("- Mainland final pop:", median(mainland.nu)), "\n")

sink()

# ----------------------------------------------------------------------------------------

# --- Test for asymmetric migration

dadi.2pop$asym.bool = dadi.2pop$is.best &
                      ((dadi.2pop$m12 - 1.96 * dadi.2pop$FIM_m12 > dadi.2pop$m21) |
                       (dadi.2pop$m21 - 1.96 * dadi.2pop$FIM_m21 > dadi.2pop$m12))

asym = dadi.2pop[which(dadi.2pop$asym.bool),]

asym.cols = c("pop1", "pop2", "m12", "m21", "FIM_m12", "FIM_m21",
              "m12.real.inds", "m21.real.inds")

asym.sort = asym[order(-abs(asym$m12.real.inds - asym$m21.real.inds)), asym.cols]

write.table(asym.sort, file="results/asymmetrical.migration.txt",
    sep="\t", quote=FALSE, row.names=FALSE)

asym.sort$less.m.donor.pop = asym.sort$pop1
asym.sort[asym.sort$m21 > asym.sort$m12,]$less.m.donor.pop =
    asym.sort$pop2[asym.sort$m21 > asym.sort$m12]
asym.sort$higher.m.donor.pop = asym.sort$pop1
asym.sort[asym.sort$m21 < asym.sort$m12,]$higher.m.donor.pop =
    asym.sort$pop2[asym.sort$m21 < asym.sort$m12]

# ----------------------------------------------------------------------------------------

dadi.2pop.out = do.call(rbind, lapply(unique(dadi.2pop$both.pops), function (pair) {

    this.pair = dadi.2pop[which(dadi.2pop$both.pops == pair & dadi.2pop$is.best),]

    this.pair.res = list()
    this.pair.res[["sites"]] = paste(
        gsub("\\(m", "(M", gsub("\\(i", "(I", simpleCap(this.pair$pop1))),
        gsub("\\(m", "(M", gsub("\\(i", "(I", simpleCap(this.pair$pop2))),
        sep=" - ")

    # Blank for the table hack
    this.pair.res[["sites_GIM"]] = ""
    this.pair.res[["sites_FIM"]] = ""

    #   this.pair.res[["ll"]] = signif.3(this.pair$ll)
    #   this.pair.res[["ll_GIM"]] = ""
    #   this.pair.res[["ll_FIM"]] = ""

    # --- Na
    Na.vals = this.pair[grepl("Na$", names(this.pair))]
    this.pair.res[["Na"]] = signif.3(Na.vals[1])
    Na.gim.tmp = as.numeric(Na.vals[1]) + c(-1,1) * 1.96 * as.numeric(Na.vals[2])
    this.pair.res[["Na_GIM"]] = paste0("(", signif.3(Na.gim.tmp[1]), ", ",
                                            signif.3(Na.gim.tmp[2]), ")")
    Na.fim.tmp = as.numeric(Na.vals[1]) + c(-1,1) * 1.96 * as.numeric(Na.vals[3])
    this.pair.res[["Na_FIM"]] = paste0("(", signif.3(Na.fim.tmp[1]), ", ",
                                            signif.3(Na.fim.tmp[2]), ")")

    # --- s
    s.vals = this.pair[grepl("^s$", names(this.pair)) |
                       grepl("[GF]IM_s$", names(this.pair))]
    this.pair.res[["s"]] = signif.3(s.vals[1])
    s.gim.tmp = as.numeric(s.vals[1]) + c(-1,1) * 1.96 * as.numeric(s.vals[2])
    this.pair.res[["s_GIM"]] = paste0("(", signif.3(s.gim.tmp[1]), ", ",
                                            signif.3(s.gim.tmp[2]), ")")
    s.fim.tmp = as.numeric(s.vals[1]) + c(-1,1) * 1.96 * as.numeric(s.vals[3])
    this.pair.res[["s_FIM"]] = paste0("(", signif.3(s.fim.tmp[1]), ", ",
                                            signif.3(s.fim.tmp[2]), ")")

    # --- NF1
    NF1.vals = NA
    this.pair.res[["NF1"]] = "-"
    this.pair.res[["NF1_GIM"]] = ""
    this.pair.res[["NF1_FIM"]] = ""

    NF1.vals = this.pair[grepl("nu1$", names(this.pair))]

    this.pair.res[["NF1"]] = signif.3(NF1.vals[1])
    NF1.gim.tmp = as.numeric(NF1.vals[1]) + c(-1,1) * 1.96 * as.numeric(NF1.vals[2])
    this.pair.res[["NF1_GIM"]] = paste0("(", signif.3(NF1.gim.tmp[1]), ", ",
                                             signif.3(NF1.gim.tmp[2]), ")")
    NF1.fim.tmp = as.numeric(NF1.vals[1]) + c(-1,1) * 1.96 * as.numeric(NF1.vals[3])
    this.pair.res[["NF1_FIM"]] = paste0("(", signif.3(NF1.fim.tmp[1]), ", ",
                                             signif.3(NF1.fim.tmp[2]), ")")

    # --- NF2
    NF2.vals = NA
    this.pair.res[["NF2"]] = "-"
    this.pair.res[["NF2_GIM"]] = ""
    this.pair.res[["NF2_FIM"]] = ""

    NF2.vals = this.pair[grepl("nu2$", names(this.pair))]

    this.pair.res[["NF2"]] = signif.3(NF2.vals[1])
    NF2.gim.tmp = as.numeric(NF2.vals[1]) + c(-1,1) * 1.96 * as.numeric(NF2.vals[2])
    this.pair.res[["NF2_GIM"]] = paste0("(", signif.3(NF2.gim.tmp[1]), ", ",
                                             signif.3(NF2.gim.tmp[2]), ")")
    NF2.fim.tmp = as.numeric(NF2.vals[1]) + c(-1,1) * 1.96 * as.numeric(NF2.vals[3])
    this.pair.res[["NF2_FIM"]] = paste0("(", signif.3(NF2.fim.tmp[1]), ", ",
                                             signif.3(NF2.fim.tmp[2]), ")")

    # --- T
    T.vals = NA
    this.pair.res[["T"]] = "-"
    this.pair.res[["T_GIM"]] = ""
    this.pair.res[["T_FIM"]] = ""

    T.vals = this.pair[grepl("T.real$", names(this.pair))]

    this.pair.res[["T"]] = signif.3(T.vals[1])
    T.gim.tmp = as.numeric(T.vals[1]) + c(-1,1) * 1.96 * as.numeric(T.vals[2])
    this.pair.res[["T_GIM"]] = paste0("(", signif.3(T.gim.tmp[1]), ", ",
                                           signif.3(T.gim.tmp[2]), ")")
    T.fim.tmp = as.numeric(T.vals[1]) + c(-1,1) * 1.96 * as.numeric(T.vals[3])
    this.pair.res[["T_FIM"]] = paste0("(", signif.3(T.fim.tmp[1]), ", ",
                                           signif.3(T.fim.tmp[2]), ")")

    # --- m12
    m12.RI.vals = NA
    this.pair.res[["m12.RI"]] = "None"
    this.pair.res[["m12.RI_GIM"]] = ""
    this.pair.res[["m12.RI_FIM"]] = ""

    m12.RI.vals = this.pair[grepl("m12.real.inds$", names(this.pair))]

    if (!is.na(m12.RI.vals[1])) {
        this.pair.res[["m12.RI"]] = signif.3(m12.RI.vals[1])
        m12.RI.gim.tmp = as.numeric(m12.RI.vals[1]) + c(-1,1) * 1.96 * as.numeric(m12.RI.vals[2])
        this.pair.res[["m12.RI_GIM"]] = paste0("(", signif.3(m12.RI.gim.tmp[1]), ", ",
                                               signif.3(m12.RI.gim.tmp[2]), ")")
        m12.RI.fim.tmp = as.numeric(m12.RI.vals[1]) + c(-1,1) * 1.96 * as.numeric(m12.RI.vals[3])
        this.pair.res[["m12.RI_FIM"]] = paste0("(", signif.3(m12.RI.fim.tmp[1]), ", ",
                                               signif.3(m12.RI.fim.tmp[2]), ")")
    }

    # --- m21
    m21.RI.vals = NA
    this.pair.res[["m21.RI"]] = "None"
    this.pair.res[["m21.RI_GIM"]] = ""
    this.pair.res[["m21.RI_FIM"]] = ""

    m21.RI.vals = this.pair[grepl("m21.real.inds$", names(this.pair))]

    if (!is.na(m21.RI.vals[1])) {
        this.pair.res[["m21.RI"]] = signif.3(m21.RI.vals[1])
        m21.RI.gim.tmp = as.numeric(m21.RI.vals[1]) + c(-1,1) * 1.96 * as.numeric(m21.RI.vals[2])
        this.pair.res[["m21.RI_GIM"]] = paste0("(", signif.3(m21.RI.gim.tmp[1]), ", ",
                                               signif.3(m21.RI.gim.tmp[2]), ")")
        m21.RI.fim.tmp = as.numeric(m21.RI.vals[1]) + c(-1,1) * 1.96 * as.numeric(m21.RI.vals[3])
        this.pair.res[["m21.RI_FIM"]] = paste0("(", signif.3(m21.RI.fim.tmp[1]), ", ",
                                               signif.3(m21.RI.fim.tmp[2]), ")")
    }

    this.pair.res[["comparison.type"]]     = this.pair$comparison.type
    this.pair.res[["comparison.type_GIM"]] = this.pair$comparison.type
    this.pair.res[["comparison.type_FIM"]] = this.pair$comparison.type

    return(this.pair.res)
}))

dadi.2pop.out.main = dadi.2pop.out[,seq(from=1, to=ncol(dadi.2pop.out), by=3)]
dadi.2pop.out.GIM  = dadi.2pop.out[,seq(from=2, to=ncol(dadi.2pop.out), by=3)]
dadi.2pop.out.FIM  = dadi.2pop.out[,seq(from=3, to=ncol(dadi.2pop.out), by=3)]

# Now using just FIM
dadi.2pop.out.sort = do.call(rbind, lapply(1:nrow(dadi.2pop.out.main), function (line.num) {
    return(rbind(dadi.2pop.out.main[line.num,],
                 dadi.2pop.out.FIM [line.num,]))
}))

dadi.2pop.out.sort = data.frame(dadi.2pop.out.sort)

names(dadi.2pop.out.sort) = c(
    "Localities",
    #"lnL",
    "$N_a$", "\\% Pop. 1 in Split",
    "Pop. 1 $\\nu_F$", "Pop. 2 $\\nu_F$",
    "Time since split",
    "$m_{12}$", "$m_{21}$",
    "comparison.type")

# Split by comparison type
dadi.2pop.out.splt = split(dadi.2pop.out.sort, unlist(dadi.2pop.out.sort$comparison.type))

# ----------------------------------------------------------------------------------------

display = c(rep('s', 9))
digits  = c(rep(0, 9))

tmp = lapply(names(dadi.2pop.out.splt), function (cmp.type) {

    caption.2pop = paste0("Results of two population demographic inference with IM model in ",
        "$\\delta$a$\\delta$i when comparing ", gsub("\\.", " to ", cmp.type), " localities.",
        "Numbers in parentheses are bounds of 95\\% confidence interval computed using ",
        "Fisher information matrix and 100 bootstrap replicates of 1 Mb from the dataset.")

    xt = xtable(dadi.2pop.out.splt[[cmp.type]][,-ncol(dadi.2pop.out.sort)],
        display=display, digits=digits,
        caption=c(caption.2pop, "$\\delta$a$\\delta$i two population IM results"),
        label=paste0("table:dadi_2pop_", gsub("\\.", "_", cmp.type)))

    align(xt) = c('r', 'r|', rep('r', 7))

    tex.out = paste0("reports/dadi.2pop.", cmp.type, ".out.tex")

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
                    tabular.environment = 'longtable', floating = FALSE,
                    hline.after=c(0, 0,
                        seq(from=2, to=nrow(dadi.2pop.out.splt[[cmp.type]]), by=2),
                        nrow(dadi.2pop.out.splt[[cmp.type]])),
                    caption.placement = "top",
                    sanitize.text.function = function(x){x},
                    size="\\fontsize{9pt}{10pt}\\selectfont")

        cat("\\end{document}", sep="\n")

    sink()
})

# ----------------------------------------------------------------------------------------

# Compile *.tex files outside of R

# cd reports
# for tex in $( ls *.tex ); do
#     pdflatex $tex
# done
# cd ..

# ----------------------------------------------------------------------------------------

dadi.2pop.m = cbind(
    do.call(rbind, lapply(dadi.2pop.out.main[,1], function (x) {
        strsplit(x, split=" - ")[[1]]
    })),
    data.frame(dadi.2pop.out.main)[,c("m12.RI", "m21.RI")],
    data.frame(dadi.2pop.out.FIM)[,c("m12.RI_FIM", "m21.RI_FIM")]
)

names(dadi.2pop.m)[1:2] = c("pop1", "pop2")

dadi.2pop.m = data.frame(sapply(dadi.2pop.m, unlist))

write.table(dadi.2pop.m, "reports/dadi_migration_matrix.txt", sep="\t", row.names=FALSE)

dadi.2pop.m$m12.RI[dadi.2pop.m$m12.RI %in% c("None", 0)] = 1e-3
dadi.2pop.m$m21.RI[dadi.2pop.m$m21.RI %in% c("None", 0)] = 1e-3

# --- Write sorted file with highest populations (for ease of access)

dadi.dists.towrite = dadi.2pop.m
dadi.dists.towrite$m12.RI = as.numeric(dadi.dists.towrite$m12.RI)
dadi.dists.towrite$m21.RI = as.numeric(dadi.dists.towrite$m21.RI)
dadi.dists.towrite$highest.m = apply(dadi.dists.towrite[,3:4], 1, max)
dadi.dists.towrite = dadi.dists.towrite[order(-dadi.dists.towrite$highest.m),]

write.table(dadi.dists.towrite, "reports/dadi_top_migration.txt",
    sep="\t", row.names=FALSE)

# ----------------------------------------------------------------------------------------

# Make matrix of migration values, m

islands.lc  = gsub("\\(m", "(M", gsub("\\(i", "(I", unlist(lapply(islands,  simpleCap))))
mainland.lc = gsub("\\(m", "(M", gsub("\\(i", "(I", unlist(lapply(mainland, simpleCap))))

sites = unique(c(dadi.2pop.m$pop1, dadi.2pop.m$pop2))

dadi.dists = rbind(data.frame(pop1=dadi.2pop.m$pop1, pop2=dadi.2pop.m$pop2, dist=dadi.2pop.m$m21.RI),
                   data.frame(pop1=dadi.2pop.m$pop2, pop2=dadi.2pop.m$pop1, dist=dadi.2pop.m$m12.RI),
                   data.frame(pop1=sites,            pop2=sites,            dist=NA))
dadi.dists = dadi.dists[order(dadi.dists$pop1, dadi.dists$pop2),]

dadi.dists.m = matrix(dadi.dists$dist, nrow=length(sites))

colnames(dadi.dists.m) = rownames(dadi.dists.m) = sites

new.order = sites[order(sites %in% mainland.lc)]
dadi.dists.m = dadi.dists.m[new.order, new.order]

# ----------------------------------------------------------------------------------------

dadi.dists$pop1 = factor(dadi.dists$pop1, levels = c(islands.lc, mainland.lc))
dadi.dists$pop2 = factor(dadi.dists$pop2, levels = c(islands.lc, mainland.lc))

dadi.dists$dist = as.numeric(gsub(",", "", dadi.dists$dist))

p = ggplot(dadi.dists, aes(pop1, pop2)) +
    geom_tile(aes(fill = log(dist)), color = "black") +
    geom_text(data=dadi.dists,
        aes(pop1, pop2,
            label = round(dist, digits=2),
            color=dadi.dists$dist),
        size=rel(5)) +
    coord_fixed() +
    labs(x='Migration Source', y='Migration Destination') +
    theme_grey(base_size=9) +
    theme(panel.background = element_blank(),
        panel.grid.major=element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14, face="bold")) +
    scale_y_discrete(limits = rev(levels(factor(dadi.dists$pop2)))) +
    scale_fill_gradientn(colors=c("white", "darkred"),
        na.value = "black")

# Add lines to separate islands and mainland
p = p +
    geom_segment(x=length(islands.lc) + 0.5,
                 xend=length(islands.lc) + 0.5,
                 y=0.5,
                 yend=length(c(islands.lc, mainland.lc)) + 0.5,
                 col='black', lwd=2) +
    geom_segment(x=0.5,
                 xend=length(c(islands.lc, mainland.lc)) + 0.5,
                 y=length(mainland.lc) + 0.5,
                 yend=length(mainland.lc) + 0.5,
                 col='black', lwd=2)

ggsave(p, file="reports/dadi_migration_matrix.pdf", height=7, width=7)

# ----------------------------------------------------------------------------------------

# --- Test for differences in split time and population size

sink("reports/dadi_2pop_is-ml_differences.time-size.txt")

# --- Median time since split
cat("Median time since split:\n")

best = dadi.2pop$is.best

t.real.inds.median.gt0 = aggregate(dadi.2pop[best,]$T.real,
    by=list(dadi.2pop[best,]$comparison.type),
    function (x) median(x, na.rm=TRUE))
t.real.inds.median.gt0

# --- Mean time since split
cat("Mean time since split:\n")

best = dadi.2pop$is.best

t.real.inds.mean.gt0 = aggregate(dadi.2pop[best,]$T.real,
    by=list(dadi.2pop[best,]$comparison.type),
    function (x) mean(x, na.rm=TRUE))
t.real.inds.mean.gt0

# --- Ancestral population size
cat("Ancestral population size:\n")

na.median.gt0 = aggregate(dadi.2pop[best,]$Na,
    by=list(dadi.2pop[best,]$comparison.type),
    function (x) median(x, na.rm=TRUE))
na.median.gt0

# --- Final population size
cat("Final population size:\n")

island.nu =   c(dadi.2pop[best & dadi.2pop$pop1 %in% islands ,]$nu1.real,
                dadi.2pop[best & dadi.2pop$pop2 %in% islands ,]$nu2.real)
mainland.nu = c(dadi.2pop[best & dadi.2pop$pop1 %in% mainland,]$nu1.real,
                dadi.2pop[best & dadi.2pop$pop2 %in% mainland,]$nu2.real)

cat(paste("- Island final pop:",   median(island.nu)), "\n")
cat(paste("- Mainland final pop:", median(mainland.nu)), "\n")

sink()
