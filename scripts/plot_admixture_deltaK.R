#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Plot ADMIXTURE delta K from Best K
# ----------------------------------------------------------------------------------------

options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly = TRUE)
deltak.in = args[1]
out.pdf   = args[2]

# Set K to NA if none specified
if (length(args) == 2) {
    args[3] = NA
}
ideal.k  = args[3]

admix_deltak = read.table(deltak.in, header=FALSE)

names(admix_deltak) = c("K", "DeltaK")

admix_deltak = admix_deltak[order(admix_deltak$K),]

pdf(file=out.pdf)

    plot(admix_deltak$DeltaK ~ admix_deltak$K, type="b", pch=19,
        main="ADMIXTURE Delta K Plot",
        xlab="K", ylab="Delta K")

    if (is.na(ideal.k)) {
        ideal.k = admix_deltak[which(admix_deltak$DeltaK == max(admix_deltak$DeltaK)),]$K
    }
    points(admix_deltak$DeltaK[which(admix_deltak$K == ideal.k)] ~ ideal.k,
        col="red", cex=2)
    text(ideal.k, admix_deltak$DeltaK[which(admix_deltak$K == ideal.k)],
        labels=round(admix_deltak$DeltaK[ideal.k], digits=3),
        cex=0.8, pos=2, col="red")

dev.off()
