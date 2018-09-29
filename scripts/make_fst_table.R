#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Make Fst table
# ----------------------------------------------------------------------------------------

library(xtable)

library(ggplot2)
library(reshape2)

fst    = read.table("data/all.pass.snp.flt.eigen.fst.se.out.fst.txt",
                    row.names=NULL)[,c(-1)]
fst.sd = read.table("data/all.pass.snp.flt.eigen.fst.se.out.sd.txt",
                    row.names=NULL)[,c(-1)]
fst.z  = read.table("data/all.pass.snp.flt.eigen.fst.se.out.fstZ.txt",
                    row.names=NULL)[,c(-1)]

# Make Z-scores symmetric
fst.z[upper.tri(fst.z)] = t(fst.z)[upper.tri(fst.z)]

# Get all into correct order

row.names(fst) = row.names(fst.sd) = row.names(fst.z) = names(fst)

islands  = c("BANDA", "BUGALAIS", "BUKASA", "NSADZI", "SSERINYA")
mainland = c("BUGALAML", "BUWAMA", "KAZZI", "KIYINDI", "MITYANA")

fst    = fst   [c(islands, mainland), c(islands, mainland)]
fst.sd = fst.sd[c(islands, mainland), c(islands, mainland)]
fst.z  = fst.z [c(islands, mainland), c(islands, mainland)]

fst    = fst    / 1000
fst.sd = fst.sd / 1000000
fst.z  = fst.z  / 1000

#   # Combine and put SD in parentheses. Also lose leading zero
#   combo = `dim<-`(sprintf('%0.3f (%0.1e)', as.matrix(fst), as.matrix(fst.sd)),
#       dim(as.matrix(fst)))
#   combo = gsub("e-0", "e-", combo)

# Now just Fst (removed SD in parentheses)
combo = `dim<-`(sprintf('%0.3f', as.matrix(fst)),
    dim(as.matrix(fst)))
combo = gsub("e-0", "e-", combo)

# Put z-score in lower triangle
combo[lower.tri(combo)] = fst.z[lower.tri(combo)]
combo[combo == 9.999] = "> 10"

colnames(combo) = row.names(combo) = c(islands, mainland)

row.names(combo)[row.names(combo) == "BUGALAIS"] = "BUGALA (I)"
row.names(combo)[row.names(combo) == "BUGALAML"] = "BUGALA (M)"

colnames(combo)[colnames(combo) == "BUGALAIS"] = "BUGALA (I)"
colnames(combo)[colnames(combo) == "BUGALAML"] = "BUGALA (M)"

islands[islands   == "BUGALAIS"] = "BUGALA (I)"
mainland[mainland == "BUGALAIS"] = "BUGALA (I)"
islands[islands   == "BUGALAML"] = "BUGALA (M)"
mainland[mainland == "BUGALAML"] = "BUGALA (M)"

row.names(combo)[row.names(combo) == "KAZZI"] = "KAAZI"
row.names(combo)[row.names(combo) == "MITYANA"] = "WAMALA"

colnames(combo)[colnames(combo) == "KAZZI"] = "KAAZI"
colnames(combo)[colnames(combo) == "MITYANA"] = "WAMALA"

mainland[mainland == "KAZZI"]   = "KAAZI"
mainland[mainland == "MITYANA"] = "WAMALA"

diag(combo) = NA

# Add names to right hand side
combo = cbind(combo, row.names(combo))
colnames(combo)[ncol(combo)] = "-"

align.str = c('l',
    rep('l', length(islands) ), "|",
    rep('l', length(mainland)), "|",
    'l')

xt = xtable(combo,
    align=align.str, digits=rep(0, 1+ncol(combo)),
    label = c("table:fst_hudson"),
    caption = c(paste("$F_{ST}$ between sites computed using Hudson's estimator.",
                      "$Z$-scores shown in lower diagonal.")))
print(xt, include.rownames=FALSE,
    hline.after=c(-1,0,0,length(islands), nrow(combo), nrow(combo)),
    size="\\fontsize{9pt}{10pt}\\selectfont",
    file = "reports/hudson_fst.tex")

# ----------------------------------------------------------------------------------------

# --- Do heatmap while we're here

combo.long = melt(combo[,-c(ncol(combo))], id.vars = c("-"))

# Change order to separate islands and mainland
combo.long$Var1 = factor(combo.long$Var1, levels = c(islands, mainland))
combo.long$Var2 = factor(combo.long$Var2, levels = c(islands, mainland))

# Set Fst for Z-score cells to 0 to fill them white
combo.long$fst = as.numeric(combo.long$value)
combo.long[as.vector(lower.tri(combo)),]$fst = 0

# Fix capitalization
simple.cap = function(x) {
    s = strsplit(x, " ")[[1]]
    lstr = paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
        sep="", collapse=" ")

    gsub("\\(i\\)", "(I)", gsub("\\(m\\)", "(M)", lstr))
}

levels(combo.long$Var1) = sapply(levels(combo.long$Var1), simple.cap)
levels(combo.long$Var2) = sapply(levels(combo.long$Var2), simple.cap)

p = ggplot(combo.long, aes(Var1, Var2)) +
    geom_tile(aes(fill = fst), color = "light grey") +
    geom_text(data=combo.long,
        aes(Var1, Var2,
           label = value),
        color=c("black", "white")[1 + as.numeric(combo.long$fst > 0.01)],
        size=rel(5)) +
    scale_fill_gradient(low = "white", high = "black") +
    coord_fixed() +
    labs(x='', y='') +
    theme_grey(base_size=9) +
    theme(panel.background = element_blank(),
        panel.grid.major=element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text=element_text(size=12, face="bold")) +
    scale_y_discrete(limits = rev(levels(factor(combo.long$Var2))))

# Add lines to separate islands and mainland
p = p +
    geom_segment(x=length(islands) + 0.5,
                 xend=length(islands) + 0.5,
                 y=0.5,
                 yend=length(mainland) + 0.5,
                 col='black', lwd=1) +
    geom_segment(x=0.5,
                 xend=length(islands) + 0.5,
                 y=length(mainland) + 0.5,
                 yend=length(mainland) + 0.5,
                 col='black', lwd=1)

ggsave(p, filename=paste0("reports/fst_heatmap_hudson.pdf"))
