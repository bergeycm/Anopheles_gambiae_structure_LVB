#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Summarize Fst bootstrapping results to get p-value
# ----------------------------------------------------------------------------------------

in.reals = list.files("results/bootstrapped_Fst/", pattern="real.log")
pairs = unique(sapply(in.reals, function (x) { strsplit(x, split="\\.")[[1]][2] }))

compute.p = function (pair) {

    in.real = paste0("results/bootstrapped_Fst/3L.", pair, ".real.log")

    real = as.numeric(system(paste0("grep 'weighted' ", in.real, " | cut -d':' -f 2"),
        intern=TRUE))

    in.boots = paste0("results/bootstrapped_Fst/3L.", pair, ".weighted.fst")
    boots = read.table(in.boots)

    p = sum(boots$V1 > real) / length(boots$V1)

    if (p == 0) {
        p = paste0("p < ", 1/length(boots$V1))
    }
}

pvals = data.frame(comparison = pairs, p = sapply(pairs, compute.p))

write.table(pvals, file="reports/Fst_pvals.txt",
    quote=FALSE, row.names=FALSE, sep="\t")
