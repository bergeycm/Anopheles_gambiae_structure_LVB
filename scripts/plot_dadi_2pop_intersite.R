#!/usr/bin/env Rscript

# ========================================================================================
# --- Parse dadi output
# ========================================================================================

library(ggplot2)
library(xtable)

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Plot 2-population model results - Between sites
# ----------------------------------------------------------------------------------------

dadi.2pop.out = list.files(path="results/dadi",
                           pattern="dadi.island.2pop.*.out",
                           full.names=TRUE)

# Skip two lines with citation info in them when bootstrapping
dadi = do.call(rbind, lapply(dadi.2pop.out, function(x) {
    read.table(x, header=TRUE, sep="\t", skip=2)
}))

# ----------------------------------------------------------------------------------------

# Make matrix of migration values, m

sites = unique(c(dadi$pop1, dadi$pop2))

dadi.dists = rbind(data.frame(pop1=dadi$pop1, pop2=dadi$pop2, dist=dadi$m12),
                   data.frame(pop1=dadi$pop2, pop2=dadi$pop1, dist=dadi$m21),
                   data.frame(pop1=sites,     pop2=sites,     dist=NA))
dadi.dists = dadi.dists[order(dadi.dists$pop1, dadi.dists$pop2),]

dadi.dists.m = matrix(dadi.dists$dist, nrow=length(sites))

colnames(dadi.dists.m) = rownames(dadi.dists.m) = sites

island.sites   = c("BANDA", "BUGALAIS", "BUKASA", "NSADZI", "SSERINYA")
mainland.sites = c("BUGALAML", "BUWAMA", "KAZZI", "KIYINDI", "MITYANA")

new.order = sites[order(sites %in% mainland.sites)]
dadi.dists.m = dadi.dists.m[new.order, new.order]

# ----------------------------------------------------------------------------------------

# Reorder
dadi.dists$pop1 = gsub("BUGALAIS", "BUGALA (I)", dadi.dists$pop1)
dadi.dists$pop1 = gsub("BUGALAML", "BUGALA (M)", dadi.dists$pop1)
dadi.dists$pop2 = gsub("BUGALAIS", "BUGALA (I)", dadi.dists$pop2)
dadi.dists$pop2 = gsub("BUGALAML", "BUGALA (M)", dadi.dists$pop2)

island.sites[island.sites     == "BUGALAIS"] = "BUGALA (I)"
mainland.sites[mainland.sites == "BUGALAML"] = "BUGALA (M)"

dadi.dists$pop1 = factor(dadi.dists$pop1, levels = c(island.sites, mainland.sites))
dadi.dists$pop2 = factor(dadi.dists$pop2, levels = c(island.sites, mainland.sites))

color="darkred"

p = ggplot(dadi.dists, aes(pop1, pop2)) +
    geom_tile(aes(fill = dist), color = "black") +
    geom_text(data=dadi.dists,
        aes(pop1, pop2,
           label = round(dist, digits=2)),
        color=c(color, "white")[1 + as.numeric(dadi.dists$dist > 1)],
        size=rel(5)) +
    scale_fill_gradient(low = "white", high = color) +
    coord_fixed() +
    labs(x='', y='') +
    theme_grey(base_size=9) +
    theme(panel.background = element_blank(),
        panel.grid.major=element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text=element_text(size=12, face="bold")) +
    scale_y_discrete(limits = rev(levels(factor(dadi.dists$pop2))))

# Add lines to separate islands and mainland
p = p +
    geom_segment(x=length(island.sites) + 0.5,
                 xend=length(island.sites) + 0.5,
                 y=0.5,
                 yend=length(c(island.sites, mainland.sites)) + 0.5,
                 col='black', lwd=2) +
    geom_segment(x=0.5,
                 xend=length(c(island.sites, mainland.sites)) + 0.5,
                 y=length(mainland.sites) + 0.5,
                 yend=length(mainland.sites) + 0.5,
                 col='black', lwd=2)

ggsave(p, file="reports/dadi_migration_matrix.pdf", height=7, width=7)

# ========================================================================================

# Make matrix of migration values, scaled by ancestral pop size theta   

dadi.dists = rbind(data.frame(pop1=dadi$pop1, pop2=dadi$pop2, dist=dadi$m12 / (dadi$theta * 2)),
                   data.frame(pop1=dadi$pop2, pop2=dadi$pop1, dist=dadi$m21 / (dadi$theta * 2)),
                   data.frame(pop1=sites,     pop2=sites,     dist=NA))
dadi.dists = dadi.dists[order(dadi.dists$pop1, dadi.dists$pop2),]

dadi.dists.m = matrix(dadi.dists$dist, nrow=length(sites))

colnames(dadi.dists.m) = rownames(dadi.dists.m) = sites

island.sites   = c("BANDA", "BUGALAIS", "BUKASA", "NSADZI", "SSERINYA")
mainland.sites = c("BUGALAML", "BUWAMA", "KAZZI", "KIYINDI", "MITYANA")

new.order = sites[order(sites %in% mainland.sites)]
dadi.dists.m = dadi.dists.m[new.order, new.order]

# ----------------------------------------------------------------------------------------

# Reorder
dadi.dists$pop1 = gsub("BUGALAIS", "BUGALA (I)", dadi.dists$pop1)
dadi.dists$pop1 = gsub("BUGALAML", "BUGALA (M)", dadi.dists$pop1)
dadi.dists$pop2 = gsub("BUGALAIS", "BUGALA (I)", dadi.dists$pop2)
dadi.dists$pop2 = gsub("BUGALAML", "BUGALA (M)", dadi.dists$pop2)

island.sites[island.sites     == "BUGALAIS"] = "BUGALA (I)"
mainland.sites[mainland.sites == "BUGALAML"] = "BUGALA (M)"

dadi.dists$pop1 = factor(dadi.dists$pop1, levels = c(island.sites, mainland.sites))
dadi.dists$pop2 = factor(dadi.dists$pop2, levels = c(island.sites, mainland.sites))

color="darkred"

p = ggplot(dadi.dists, aes(pop1, pop2)) +
    geom_tile(aes(fill = dist), color = "black") +
    geom_text(data=dadi.dists,
        aes(pop1, pop2,
           label = round(dist, digits=2)),
        color=c(color, "white")[1 + as.numeric(dadi.dists$dist > 1)],
        size=rel(5)) +
    scale_fill_gradient(low = "white", high = color) +
    coord_fixed() +
    labs(x='', y='') +
    theme_grey(base_size=9) +
    theme(panel.background = element_blank(),
        panel.grid.major=element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text=element_text(size=12, face="bold")) +
    scale_y_discrete(limits = rev(levels(factor(dadi.dists$pop2))))

# Add lines to separate islands and mainland
p = p +
    geom_segment(x=length(island.sites) + 0.5,
                 xend=length(island.sites) + 0.5,
                 y=0.5,
                 yend=length(c(island.sites, mainland.sites)) + 0.5,
                 col='black', lwd=2) +
    geom_segment(x=0.5,
                 xend=length(c(island.sites, mainland.sites)) + 0.5,
                 y=length(mainland.sites) + 0.5,
                 yend=length(mainland.sites) + 0.5,
                 col='black', lwd=2)

ggsave(p, file="reports/dadi_migration_matrix_scaled.pdf", height=7, width=7)

# ----------------------------------------------------------------------------------------

# prior_onegrow_mig:
#     https://github.com/paulirish/dadi/blob/master/dadi/Demographics2D.py
# Model with growth, split, bottleneck in pop2 , exp recovery, migration
#     nu1F: The ancestral population size after growth. (Its initial size is
#           defined to be 1.)
#     nu2B: The bottleneck size for pop2
#     nu2F: The final size for pop2
#     m: The scaled migration rate
#     Tp: The scaled time between ancestral population growth and the split.
#     T: The time between the split and present

dadi.2pop$Na = dadi.2pop$theta / (4 * mu * L.chr3)
# Size of ancestral pop after split
dadi.2pop$nu1F.real = dadi.2pop$Na * dadi.2pop$nu1F
# Size of pop2 - bottleneck
dadi.2pop$nu2B.real = dadi.2pop$Na * dadi.2pop$nu2B
# Size of pop2 - final
dadi.2pop$nu2F.real = dadi.2pop$Na * dadi.2pop$nu2F
# Convert T to real time in years
dadi.2pop$T.real  = 2 * dadi.2pop$Na * dadi.2pop$T  * g
dadi.2pop$Tp.real = 2 * dadi.2pop$Na * dadi.2pop$Tp * g

# And get proportion of migrants
# "fraction of individuals in each generation in population i who are new migrants from
#  population j." Multiply by population size to get number of individuals each generation
#  that are migrating from population j to population i."
dadi.2pop$m.real = dadi.2pop$m / (2 * dadi.2pop$Na)
                                             
dadi.2pop$m.real.inds = dadi.2pop$m.real * dadi.2pop$Na
# Isn't this the equivalent of dividing m by 2?
                                             

dadi.2pop.out = dadi.2pop[,c(5,6,7,11,17,23,8,14,9,15,10,16,12,18,13,19,29)]
names(dadi.2pop.out) = c("Log Likelihood", "AIC",
    "$\theta$", "$m$", "$m$ GIM uncertainty", "$N_a",
    "Post-growth anc. size", "$\nu 1F$ GIM uncertainty",
    "Pop2 bottleneck size", "$\nu 2B$ GIM uncertainty",
    "Pop. 2 Final Size", "$\nu 2F$ GIM uncertainty",
    "Time from growth to split", "$Tp$ GIM uncertainty",
    "Time since split", "$T$ GIM uncertainty",
    "Migrants")
#dadi.2pop.out$note = c("Including Bugala", "Excluding Bugala")

#display = c('s', rep('f', 11), 's')
#digits  = c(0, 3, 0, 0, -3, 0, 0, 0, 0, 0, 0, -3, 0)
#xt = xtable(dadi.2pop.out, display=display, digits=digits)
#align(xt) = rep('r', 13)

dadi.2pop.out.t = data.frame(row.names(t(dadi.2pop.out)), t(dadi.2pop.out))
#names(dadi.2pop.out.t) = dadi.2pop.out$note

xt = xtable(dadi.2pop.out.t,
    display=c('s','s','s'), digits=c(0,0,0))
align(xt) = c('l','r','r')

tex.out = "reports/dadi.2pop.out.tex"

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
                include.rownames = TRUE,
                tabular.environment = 'longtable', floating = FALSE,
                caption.placement = "top")

    cat("\\end{document}", sep="\n")

sink()

# ----------------------------------------------------------------------------------------

# Compile *.tex files outside of R

# cd reports
# for tex in $( ls *.tex ); do
#     pdflatex $tex
# done
# cd ..
