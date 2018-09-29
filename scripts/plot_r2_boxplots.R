#!/usr/bin/env Rscript

# ========================================================================================
# --- Plot LD measure, r^2
# ========================================================================================

library(ggplot2)
library(xtable)

options(stringsAsFactors=FALSE)

# Just 3 for now
r2.out = list.files(path="results/", pattern="chr3.*\\.*.list.hap.ld$", full.names=TRUE)

r2.all = do.call(rbind, lapply(r2.out, function(x) {
    this.r2 = read.table(x, header=TRUE, sep="\t")
    this.r2$pop = x
    this.r2
}))

r2.all$pop = gsub(".*chr3[LR]\\.(.*)\\.list.hap.ld", "\\1", r2.all$pop)

# Remove unsplit Bugala
r2.all = r2.all[r2.all$pop != "BUGALA",]

r2.all$island_mainland = "Island"
is.sites = c("BANDA", "BUGALAIS", "BUKASA", "NSADZI", "SSERINYA")
ml.sites = c("BUGALAML", "BUWAMA", "KAZZI", "KIYINDI", "MITYANA")
r2.all[r2.all$pop %in% ml.sites,]$island_mainland = "Mainland"

r2.all[r2.all$pop == "BUGALAIS",]$pop = "BUGALA (I)"
r2.all[r2.all$pop == "BUGALAML",]$pop = "BUGALA (M)"
is.sites[is.sites == "BUGALAIS"] = "BUGALA (I)"
ml.sites[ml.sites == "BUGALAML"] = "BUGALA (M)"

# --- Correct for differing sample sizes by subtracting 1/n
# --- where n is the number of sampled chromosomes
r2.all$R.2.adj = r2.all$R.2 - (1 / r2.all$N_CHR)

p = ggplot(r2.all, aes(factor(pop, levels=rev(c(is.sites, ml.sites))), R.2.adj,
        color=factor(island_mainland, levels=c("Mainland", "Island")))) +
    geom_boxplot(outlier.shape=NA) +
    theme_bw() +
    xlab("") +
    ylab(expression("R"^2)) +
    scale_color_discrete(guide=FALSE) +
    scale_y_continuous(limits = (c(0,0.25))) +
    coord_flip()

ggsave(p, filename="reports/r2_boxplot.pdf", height=3, width=4)

# Also plot island vs. mainland
p = ggplot(r2.all, aes(island_mainland, R.2.adj,
        color=factor(island_mainland, levels=rev(c("Mainland", "Island"))))) +
    geom_boxplot(outlier.shape=NA) +
    theme_bw() +
    xlab("") +
    ylab(expression("R"^2)) +
    scale_color_discrete(guide=FALSE) +
    scale_y_continuous(limits = (c(0,0.25))) +
    coord_flip()

ggsave(p, filename="reports/r2_boxplot.mainland-island.pdf", height=3, width=4)

# Plot LD decay

simple.cap = function(x) {
    s = strsplit(x, " ")[[1]]
    lstr = paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
        sep="", collapse=" ")

    gsub("\\(i\\)", "(I)", gsub("\\(m\\)", "(M)", lstr))
}

r2.all$pop = as.vector(sapply(r2.all$pop, simple.cap))
is.sites = as.vector(sapply(is.sites, simple.cap))
ml.sites = as.vector(sapply(ml.sites, simple.cap))

r2.all$dist = abs(r2.all$POS2 - r2.all$POS1)

p = ggplot(r2.all, aes(dist, R.2.adj,
        color=factor(island_mainland, levels=c("Mainland", "Island")))) +
    stat_smooth(se = FALSE,  aes(fill = pop), lwd=0.5) +
    guides(fill=FALSE) +
    scale_color_discrete(name="") +
    ylab(expression("R"^2)) +
    xlab("Physical distance (bp)") +
    scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x)),
        limits = c(1, max(r2.all$dist)),
        expand=c(0,0)
    ) +
    annotation_logticks(sides = "b") +
    theme_bw(base_size=15) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black"),
          legend.justification=c(1,1),
          legend.position=c(0.9,0.9),
          legend.title=element_blank(),
          legend.background = element_rect(color    = 'black',
                                           fill     = 'white',
                                           linetype = 'solid',
                                           size     = 0))

ggsave(p, filename="reports/r2_decay.pdf", height=3, width=4)
