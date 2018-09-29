#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Compute moving average for XP-EHH
# ----------------------------------------------------------------------------------------

args = commandArgs(trailingOnly = TRUE)
in.file = args[1]  # e.g. results/selscan/xp-ehh.BUWAMA-KAZZI.X.xpehh.out.norm

write(paste0("Averaging file ", in.file), stderr())

window.size = 10000

xpehh = read.table(in.file, header=TRUE, row.names=NULL)

# Infer chr from input file name
chr = strsplit(in.file, "\\.")[[1]][3]

xpehh$row.names = paste(chr, xpehh$pos, sep=":")
names(xpehh)[1] = "locus"

# Figure out end of last window
upper.bound = window.size * ceiling(max(xpehh$pos) / window.size)

win.mins = seq(from=1, to=upper.bound, by=window.size)
averages = vector(length = length(win.mins))

for (i in 1:length(averages)) {
    averages[i] = mean(xpehh[xpehh$pos >= win.mins[i] &
        xpehh$pos < win.mins[i+1],]$normxpehh)
    if (i %% 10 == 0) {
        write(paste0("On averaging step ", i, " of ", length(averages), "..."),
            stderr())
    }
}

xpehh = data.frame(cbind(win.mins, averages))
names(xpehh) = c("pos", "normxpehh")

write.table(xpehh, paste0(in.file, ".avg"), quote=FALSE, sep="\t", row.names=FALSE)
