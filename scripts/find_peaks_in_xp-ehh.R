#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Find peaks in XP-EHH
# ----------------------------------------------------------------------------------------

# module load r/3.4
# Rscript scripts/find_peaks_in_xp-ehh.R > results/xpehh_peaks.txt

cutoff = 0.999

# ----------------------------------------------------------------------------------------

# --- Bump around 34Mb on 2L

xp = read.table(pipe(paste0("awk '{ if ($1 == \"id\" || ",
                            "$2 > 31000000 && $2 < 36000000) ",
                            "print $2,$9 }' ",
                            "results/selscan/xp-ehh.BUKASA-MITYANA.2L.xpehh.out.norm")),
                sep=" ", header=TRUE)

extreme.cutoff = quantile(xp$normxpehh, probs=c(cutoff))

xp$is.outlier = abs(xp$normxpehh) > extreme.cutoff
mid.extreme = median(xp[xp$is.outlier,]$pos)

xp = rbind(xp[xp$is.outlier,],
           xp[sample(nrow(xp[!xp$is.outlier,]), size=1000),])

pdf("reports/xpehh_peak_2L_34Mb.pdf")
    plot(xp$normxpehh ~ xp$pos,
        col=c("grey", "red")[1+xp$is.outlier],
        main="Bump around 34Mb on 2L",
        sub=paste0("Peak: ", mid.extreme),
        xlab="Position", ylab="Normalized XP-EHH")
    abline(h = c(1,-1) * extreme.cutoff, col = 'grey', lty=3)
    abline(v = mid.extreme, col = 'red', lwd=2)
dev.off()

cat(paste0("2L_near34Mb\t2L:", round(mid.extreme)), "\n")

# ----------------------------------------------------------------------------------------

# --- Bump around 44Mb on 2L

xp = read.table(pipe(paste0("awk '{ if ($1 == \"id\" || ",
                            "$2 > 41000000 && $2 < 46000000) ",
                            "print $2,$9 }' ",
                            "results/selscan/xp-ehh.BANDA-BUKASA.2L.xpehh.out.norm")),
                sep=" ", header=TRUE)

extreme.cutoff = quantile(xp$normxpehh, probs=c(cutoff))

xp$is.outlier = abs(xp$normxpehh) > extreme.cutoff
mid.extreme = median(xp[xp$is.outlier,]$pos)

xp = rbind(xp[xp$is.outlier,],
           xp[sample(nrow(xp[!xp$is.outlier,]), size=1000),])

pdf("reports/xpehh_peak_2L_44Mb.pdf")
    plot(xp$normxpehh ~ xp$pos,
        col=c("grey", "red")[1+xp$is.outlier],
        main="Bump around 44Mb on 2L",
        sub=paste0("Peak: ", mid.extreme),
        xlab="Position", ylab="Normalized XP-EHH")
    abline(h = c(1,-1) * extreme.cutoff, col = 'grey', lty=3)
    abline(v = mid.extreme, col = 'red', lwd=2)
dev.off()

cat(paste0("2L_near44Mb\t2L:", round(mid.extreme)), "\n")

# ----------------------------------------------------------------------------------------

# --- CYP cluster on 2R

xp = read.table(pipe(paste0("awk '{ if ($1 == \"id\" || ",
                            "$2 > 26000000 && $2 < 31000000) ",
                            "print $2,$9 }' ",
                            "results/selscan/xp-ehh.BANDA-BUKASA.2R.xpehh.out.norm")),
                sep=" ", header=TRUE)

extreme.cutoff = quantile(xp$normxpehh, probs=c(cutoff))

xp$is.outlier = abs(xp$normxpehh) > extreme.cutoff
mid.extreme = median(xp[xp$is.outlier,]$pos)

xp = rbind(xp[xp$is.outlier,],
           xp[sample(nrow(xp[!xp$is.outlier,]), size=1000),])

pdf("reports/xpehh_peak_2R_CYP.pdf")
    plot(xp$normxpehh ~ xp$pos,
        col=c("grey", "red")[1+xp$is.outlier],
        main="CYP cluster on 2R",
        sub=paste0("Peak: ", mid.extreme),
        xlab="Position", ylab="Normalized XP-EHH")
    abline(h = c(1,-1) * extreme.cutoff, col = 'grey', lty=3)
    abline(v = mid.extreme, col = 'red', lwd=2)
dev.off()

cat(paste0("2R_CYP\t2R:", round(mid.extreme)), "\n")

# ----------------------------------------------------------------------------------------

# --- Bump near 57Mb on 2R

xp = read.table(pipe(paste0("awk '{ if ($1 == \"id\" || ",
                            "$2 > 56500000 && $2 < 57500000) ",
                            "print $2,$9 }' ",
                            "results/selscan/xp-ehh.BANDA-BUWAMA.2R.xpehh.out.norm")),
                sep=" ", header=TRUE)

extreme.cutoff = quantile(xp$normxpehh, probs=c(cutoff))

xp$is.outlier = abs(xp$normxpehh) > extreme.cutoff
mid.extreme = median(xp[xp$is.outlier,]$pos)

xp = rbind(xp[xp$is.outlier,],
           xp[sample(nrow(xp[!xp$is.outlier,]), size=1000),])

pdf("reports/xpehh_peak_2R_57Mb.pdf")
    plot(xp$normxpehh ~ xp$pos,
        col=c("grey", "red")[1+xp$is.outlier],
        main="Bump around 53Mb on 2R",
        sub=paste0("Peak: ", mid.extreme),
        xlab="Position", ylab="Normalized XP-EHH")
    abline(h = c(1,-1) * extreme.cutoff, col = 'grey', lty=3)
    abline(v = mid.extreme, col = 'red', lwd=2)
dev.off()

cat(paste0("2R_near57Mb\t2R:", round(mid.extreme)), "\n")

# ----------------------------------------------------------------------------------------

# --- Bump near 7Mb on 3L

xp = read.table(pipe(paste0("awk '{ if ($1 == \"id\" || ",
                            "$2 > 500000 && $2 < 20000000) ",
                            "print $2,$9 }' ",
                            "results/selscan/xp-ehh.BUKASA-SSERINYA.3L.xpehh.out.norm")),
                sep=" ", header=TRUE)

extreme.cutoff = quantile(xp$normxpehh, probs=c(cutoff))

xp$is.outlier = abs(xp$normxpehh) > extreme.cutoff
mid.extreme = median(xp[xp$is.outlier,]$pos)

xp = rbind(xp[xp$is.outlier,],
           xp[sample(nrow(xp[!xp$is.outlier,]), size=1000),])

pdf("reports/xpehh_peak_3L_7Mb.pdf")
    plot(xp$normxpehh ~ xp$pos,
        col=c("grey", "red")[1+xp$is.outlier],
        main="Bump around 7Mb on 3L",
        sub=paste0("Peak: ", mid.extreme),
        xlab="Position", ylab="Normalized XP-EHH")
    abline(h = c(1,-1) * extreme.cutoff, col = 'grey', lty=3)
    abline(v = mid.extreme, col = 'red', lwd=2)
dev.off()

cat(paste0("3L_7Mb\t3L:", round(mid.extreme)), "\n")

# ----------------------------------------------------------------------------------------

# --- CYP cluster on 3L

xp = read.table(pipe(paste0("awk '{ if ($1 == \"id\" || ",
                            "$2 > 14200000 && $2 < 14700000) ",
                            "print $2,$9 }' ",
                            "results/selscan/xp-ehh.BUGALAIS-KIYINDI.3L.xpehh.out.norm")),
                sep=" ", header=TRUE)

extreme.cutoff = quantile(xp$normxpehh, probs=c(cutoff))

xp$is.outlier = abs(xp$normxpehh) > extreme.cutoff
mid.extreme = median(xp[xp$is.outlier,]$pos)

xp = rbind(xp[xp$is.outlier,],
           xp[sample(nrow(xp[!xp$is.outlier,]), size=1000),])

pdf("reports/xpehh_peak_3L_CYP.pdf")
    plot(xp$normxpehh ~ xp$pos,
        col=c("grey", "red")[1+xp$is.outlier],
        main="CYP cluster on 3L",
        sub=paste0("Peak: ", mid.extreme),
        xlab="Position", ylab="Normalized XP-EHH")
    abline(h = c(1,-1) * extreme.cutoff, col = 'grey', lty=3)
    abline(v = mid.extreme, col = 'red', lwd=2)
dev.off()

cat(paste0("3L_CYP\t3L:", round(mid.extreme)), "\n")

# ----------------------------------------------------------------------------------------

# --- GSTE on 3R

xp = read.table(pipe(paste0("awk '{ if ($1 == \"id\" || ",
                            "$2 > 28000000 && $2 < 29000000) ",
                            "print $2,$9 }' ",
                            "results/selscan/xp-ehh.KAZZI-MITYANA.3R.xpehh.out.norm")),
                sep=" ", header=TRUE)

extreme.cutoff = quantile(xp$normxpehh, probs=c(cutoff))

xp$is.outlier = abs(xp$normxpehh) > extreme.cutoff
mid.extreme = median(xp[xp$is.outlier,]$pos)

xp = rbind(xp[xp$is.outlier,],
           xp[sample(nrow(xp[!xp$is.outlier,]), size=1000),])

pdf("reports/xpehh_peak_3R_GSTE.pdf")
    plot(xp$normxpehh ~ xp$pos,
        col=c("grey", "red")[1+xp$is.outlier],
        main="GSTE cluster on 3R",
        sub=paste0("Peak: ", mid.extreme),
        xlab="Position", ylab="Normalized XP-EHH")
    abline(h = c(1,-1) * extreme.cutoff, col = 'grey', lty=3)
    abline(v = mid.extreme, col = 'red', lwd=2)
dev.off()

cat(paste0("3R_GSTE\t3R:", round(mid.extreme)), "\n")

# ----------------------------------------------------------------------------------------

# --- Bump around 46Mb on 3R

xp = read.table(pipe(paste0("awk '{ if ($1 == \"id\" || ",
                            "$2 > 43000000 && $2 < 47000000) ",
                            "print $2,$9 }' ",
                            "results/selscan/xp-ehh.BANDA-KAZZI.3R.xpehh.out.norm")),
                sep=" ", header=TRUE)

extreme.cutoff = quantile(xp$normxpehh, probs=c(cutoff))

xp$is.outlier = abs(xp$normxpehh) > extreme.cutoff
mid.extreme = median(xp[xp$is.outlier,]$pos)

xp = rbind(xp[xp$is.outlier,],
           xp[sample(nrow(xp[!xp$is.outlier,]), size=1000),])

pdf("reports/xpehh_peak_3R_46Mb.pdf")
    plot(xp$normxpehh ~ xp$pos,
        col=c("grey", "red")[1+xp$is.outlier],
        main="Bump around 46Mb on 3R",
        sub=paste0("Peak: ", mid.extreme),
        xlab="Position", ylab="Normalized XP-EHH")
    abline(h = c(1,-1) * extreme.cutoff, col = 'grey', lty=3)
    abline(v = mid.extreme, col = 'red', lwd=2)
dev.off()

cat(paste0("3R_near46Mb\t3R:", round(mid.extreme)), "\n")

# ----------------------------------------------------------------------------------------

# --- Bump around 9Mb on X

xp = read.table(pipe(paste0("awk '{ if ($1 == \"id\" || ",
                            "$2 > 6000000 && $2 < 12000000) ",
                            "print $2,$9 }' ",
                            "results/selscan/xp-ehh.BANDA-BUWAMA.X.xpehh.out.norm")),
                sep=" ", header=TRUE)

extreme.cutoff = quantile(xp$normxpehh, probs=c(cutoff))

xp$is.outlier = abs(xp$normxpehh) > extreme.cutoff
mid.extreme = median(xp[xp$is.outlier,]$pos)

xp = rbind(xp[xp$is.outlier,],
           xp[sample(nrow(xp[!xp$is.outlier,]), size=1000),])

pdf("reports/xpehh_peak_X_9Mb.pdf")
    plot(xp$normxpehh ~ xp$pos,
        col=c("grey", "red")[1+xp$is.outlier],
        main="Bump around 9Mb on X",
        sub=paste0("Peak: ", mid.extreme),
        xlab="Position", ylab="Normalized XP-EHH")
    abline(h = c(1,-1) * extreme.cutoff, col = 'grey', lty=3)
    abline(v = mid.extreme, col = 'red', lwd=2)
dev.off()

cat(paste0("X_near9Mb\tX:", round(mid.extreme)), "\n")

# ----------------------------------------------------------------------------------------

# --- CYP9K1 on X

xp = read.table(pipe(paste0("awk '{ if ($1 == \"id\" || ",
                            "$2 > 13000000 && $2 < 17000000) ",
                            "print $2,$9 }' ",
                            "results/selscan/xp-ehh.BUGALAIS-KAZZI.X.xpehh.out.norm")),
                sep=" ", header=TRUE)

extreme.cutoff = quantile(xp$normxpehh, probs=c(cutoff))

xp$is.outlier = abs(xp$normxpehh) > extreme.cutoff
mid.extreme = median(xp[xp$is.outlier,]$pos)

xp = rbind(xp[xp$is.outlier,],
           xp[sample(nrow(xp[!xp$is.outlier,]), size=1000),])

pdf("reports/xpehh_peak_X_CYP.pdf")
    plot(xp$normxpehh ~ xp$pos,
        col=c("grey", "red")[1+xp$is.outlier],
        main="CYP9K1 on X",
        sub=paste0("Peak: ", mid.extreme),
        xlab="Position", ylab="Normalized XP-EHH")
    abline(h = c(1,-1) * extreme.cutoff, col = 'grey', lty=3)
    abline(v = mid.extreme, col = 'red', lwd=2)
dev.off()

cat(paste0("X_CYP9K1\tX:", round(mid.extreme)), "\n")

# ----------------------------------------------------------------------------------------

# --- Bump around 16Mb on X

xp = read.table(pipe(paste0("awk '{ if ($1 == \"id\" || ",
                            "$2 > 15000000 && $2 < 18000000) ",
                            "print $2,$9 }' ",
                            "results/selscan/xp-ehh.BANDA-NSADZI.X.xpehh.out.norm")),
                sep=" ", header=TRUE)

extreme.cutoff = quantile(xp$normxpehh, probs=c(cutoff))

xp$is.outlier = abs(xp$normxpehh) > extreme.cutoff
mid.extreme = median(xp[xp$is.outlier,]$pos)

xp = rbind(xp[xp$is.outlier,],
           xp[sample(nrow(xp[!xp$is.outlier,]), size=1000),])

pdf("reports/xpehh_peak_X_17Mb.pdf")
    plot(xp$normxpehh ~ xp$pos,
        col=c("grey", "red")[1+xp$is.outlier],
        main="Bump around 16Mb on X",
        sub=paste0("Peak: ", mid.extreme),
        xlab="Position", ylab="Normalized XP-EHH")
    abline(h = c(1,-1) * extreme.cutoff, col = 'grey', lty=3)
    abline(v = mid.extreme, col = 'red', lwd=2)
dev.off()

cat(paste0("X_near16Mb\tX:", round(mid.extreme)), "\n")
