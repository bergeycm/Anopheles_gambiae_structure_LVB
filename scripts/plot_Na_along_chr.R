#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Plot Na across chromsome 3
# ----------------------------------------------------------------------------------------

library(ggplot2)

options(stringsAsFactors=FALSE)

na.out = list.files(path="results/", pattern="dadi.chr3.island.*.theta_along_chr3.txt",
    full.names=TRUE)

na = do.call(rbind, lapply(na.out, function (x) {
    tmp = read.table(x)
    tmp$pop = strsplit(x, split="\\.")[[1]][4]
    tmp
}))

na$chr = do.call(c, lapply(na$V1, function (x) { strsplit(x, split="_")[[1]][1] }))

na$V1 = gsub(".*_", "", na$V1)
na$V2 = gsub(".*_", "", na$V2)

names(na) = c("start", "end", "theta", "pop", "chr")

na$start = as.numeric(na$start)

# --- Fix site names
na$pop[na$pop == "BUGALAIS"] = "BUGALA (I)"
na$pop[na$pop == "BUGALAML"] = "BUGALA (M)"
na$pop[na$pop == "KAZZI"]    = "KAAZI"
na$pop[na$pop == "MITYANA"]  = "WAMALA"

islands  = c("BANDA", "BUGALA (I)", "BUKASA", "NSADZI", "SSERINYA")
mainland = c("BUGALA (M)", "BUWAMA", "KAAZI", "KIYINDI", "WAMALA")

na$pop = factor(na$pop, levels=c(mainland, islands))

na$island_mainland = "Mainland"
na[na$pop %in% islands,]$island_mainland = "Island"

na$chr = factor(na$chr, levels = c("3R", "3L"))

p = ggplot(na, aes(start, theta, group=pop,
        col=factor(island_mainland, levels=c("Mainland", "Island")))) +
    geom_line() +
    facet_grid(.~chr, scales="free_x", space="free_x") +
    #geom_vline(data=data.frame(start=28597652, chr="3R"),
    #    aes(xintercept=start), lty=2, col="black") +
    xlab("Genome position (Mb)") +
    ylab(expression(paste(theta~"inferred with"~delta,"a",delta,
        "i in 1Mb-long random windows"))) +
    scale_x_continuous(labels = function(x) { x/1000000 },
                       breaks = seq(from=0e6, to=60e6, by=10e6),
                       expand = c(0,0)) +
    scale_color_manual(values=c("#F8766D", "#00BFC4"),
                       name="Site type",
                       breaks=c("Mainland", "Island")) +
    theme_bw() +
    theme(strip.background = element_rect(fill="white"))

ggsave(p, file="reports/theta_along_chromosome.pdf", width=12, height=6)
