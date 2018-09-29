#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Find CSS outliers
# ----------------------------------------------------------------------------------------

css = read.table("results/css.all.txt")

names(css) = c("chr", "bin.start", "pop1", "pop2", "css")

cutoff = quantile(css$css, 0.99)

# ----------------------------------------------------------------------------------------

is.sites = c("BANDA", "BUGALAIS", "BUKASA", "NSADZI", "SSERINYA")
ml.sites = c("BUGALAML", "BUWAMA", "KAZZI", "KIYINDI", "MITYANA")

# ----------------------------------------------------------------------------------------

# ignore if het   

# ID singleton sweeps (or doubleton) = only comparisons with one site outliers   

# ----------------------------------------------------------------------------------------

# --- ID island-mainland comparisons

is.ct = rowSums(cbind(css$pop1 %in% is.sites,
                      css$pop2 %in% is.sites))
css$is.ml = is.ct == 1
css$ml.ml = is.ct == 0

bins = unique(css[,c(1:2)])

assess.bin.for.is.ml = function (chr, bin.start) {

    this.bin.df = css[css$chr == chr & css$bin.start == bin.start,]

    conting.tbl = table(factor(this.bin.df$is.ml,        levels = c(TRUE, FALSE)),
                        factor(this.bin.df$css > cutoff, levels = c(TRUE, FALSE)))

    chi = chisq.test(conting.tbl)

    # Proportion of IS-ML comparisons that are high
    prop.isml.high    = conting.tbl[1] / (conting.tbl[1] + conting.tbl[3])
    # Proportion of non-IS-ML comparisons that are high
    prop.notisml.high = conting.tbl[2] / (conting.tbl[2] + conting.tbl[4])

    # Proportion of ML-ML comparisons that are high
    conting.tbl.mlml = table(factor(this.bin.df$ml.ml,        levels = c(TRUE, FALSE)),
                             factor(this.bin.df$css > cutoff, levels = c(TRUE, FALSE)))
    prop.mlml.high = conting.tbl.mlml[1] / (conting.tbl.mlml[1] + conting.tbl.mlml[3])

    res = c(chr, bin.start, as.vector(conting.tbl),
        prop.isml.high, prop.notisml.high, prop.mlml.high, chi$p.value)

    names(res) = c("chr", "bin.start",
        "notisml.nothigh", "isml.nothigh", "notisml.high", "isml.high",
        "prop.isml.high", "prop.notisml.high", "prop.mlml.high",
        "chi.p")

    return(res)
}

isml.test.res = data.frame(do.call(rbind, lapply(1:nrow(bins), function (bin.idx) {
    assess.bin.for.is.ml(bins[bin.idx,1], bins[bin.idx,2])
})))

save.image('css.outlier.finding.Rdata')

isml.test.res$chi.p = as.numeric(isml.test.res$chi.p)

high = isml.test.res[which(isml.test.res$chi.p < 0.01),]

write.table(high, file="tmp.high.isml.chisq.txt",
    sep="\t", quote=FALSE, row.names=FALSE)

very.high = high[high$chi.p < 1e-4,]
very.high[order(very.high$chr, very.high$bin.start),c(1:2,7:8,10)]

