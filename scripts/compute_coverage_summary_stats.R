#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Compute coverage average, SD, range
# ----------------------------------------------------------------------------------------

cov = read.table("reports/all.pass.snp.flt.idepth", header=TRUE)

cat(paste("Average coverage", mean(cov$MEAN_DEPTH)), "\n")

cat(paste("S.D. coverage", sd(cov$MEAN_DEPTH)), "\n")

cat(paste("Range coverage", range(cov$MEAN_DEPTH)[1], range(cov$MEAN_DEPTH)[2]), "\n")