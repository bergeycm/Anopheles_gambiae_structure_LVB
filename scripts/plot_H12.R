#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Plot H12
# ----------------------------------------------------------------------------------------

library(ggplot2)

chrs = c("2L", "2R", "3L", "3R", "X")

h12 = do.call(rbind, lapply(chrs, function (chr) {

    h12.this = read.table(paste0("results/chr", chr, ".pass.snp.phased.H12.out.txt"),
        sep="\t")
    h12.this$chr = chr
    h12.this
}))

names(h12) = c("win.center", "win.left", "win.right", "k.haplotypes", "hap.freq.spec",
               "ind", "H1.het", "H2", "H12", "H2_H1", "chr")

h12.sm = unique(rbind(h12[seq(from=1, to=nrow(h12), by=1000),],
                      h12[h12$H12 >= quantile(h12$H12, 0.90),]))

p = ggplot(h12.sm, aes(win.center, H12, color=k.haplotypes)) +
    geom_line() +
    facet_grid(chr ~ .) +
    xlab("Window center") +
    guides(color=FALSE) +
    scale_x_continuous(labels=function(x) paste(x/1000000, "Mb"), expand=c(0,0)) +
    theme_bw()

ggsave("reports/H12.pdf", p, height=8, width=12)
