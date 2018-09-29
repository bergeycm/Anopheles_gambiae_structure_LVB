#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ========================================================================================
# --- Do DAPC
# ========================================================================================

library(adegenet)
library(ape)
library(RColorBrewer)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
in.prefix = args[1]   # e.g. "results/dapc_subset_snps.chr3L"

# ----------------------------------------------------------------------------------------
# --- Read in data
# ----------------------------------------------------------------------------------------

# Already created high Fst SNPs in PLINK RAW format so that it can be read by read.PLINK()
# as a genlight object

gen = read.PLINK(paste0(in.prefix, ".raw"),
    map.file = paste0(in.prefix, ".map"),
    parallel = TRUE, n.cores = 12)

# Remove ssese99.PE from gen object
if ("ssese99.PE" %in% gen$ind.names) {
    gen = gen[-which(gen$ind.names == "ssese99.PE"),]
}

# Fix SNP names
map.file = read.table(paste0(in.prefix, ".map"))
snp.names = map.file$V2

# Test
# gen = gen[,1:10000]
# map.file = map.file[1:10000,]
# snp.names = map.file$V2

# ----------------------------------------------------------------------------------------
# --- Do PCA
# ----------------------------------------------------------------------------------------

pca = glPca(gen, nf=10,
    parallel = TRUE, n.cores = 12)

# --- Set pop variable equal to site

ind.info = read.table("data/ssese_individual_info_simple_bugala_split.txt")

ind.info$V3[ind.info$V3 == "BUGALAIS"] = "BUGALA (I)"
ind.info$V3[ind.info$V3 == "BUGALAML"] = "BUGALA (M)"

is.sites = c("BANDA", "BUGALA (I)", "BUKASA", "NSADZI", "SSERINYA")
ml.sites = c("BUGALA (M)", "BUWAMA", "KAZZI", "KIYINDI", "MITYANA")

pops.pca = do.call(c, lapply(row.names(pca$scores), function(x) {
    as.character(ind.info[ind.info$V2 == x, 3])
}))

gen$pop = factor(pops.pca, levels = c(is.sites, ml.sites))

# --- Plot PCA and tree with colors to represent position

pdf(paste0(in.prefix, ".pca.colorplot.pdf"), height=20, width=20)

    pca.cols = colorplot(pca$scores,pca$scores, transp=TRUE, cex=4)
    abline(h=0, v=0, col="grey")
    # add.scatter.eig(pca$eig[1:40],2,1,2, posi="topright", inset=.05, ratio=.3)

    tre = nj(dist(as.matrix(gen)))

    plot.new()

    par(mar=c(0, 0, 0, 0))
    plot(tre, no.margin=TRUE, direction='upwards', cex=0.5)
    tiplabels(pch=20, col=pca.cols, cex=1)

dev.off()

# ----------------------------------------------------------------------------------------
# --- Cluster
# ----------------------------------------------------------------------------------------

grp = find.clusters(gen, n.pca = 100,
    choose.n.clust = FALSE, criterion = "diffNgroup", max.n.clust = 20, glPca=pca)

ideal.k = length(table(grp$grp))

# --- Plot BIC for each cluster

bic.df = data.frame(K = 1:length(grp$Kstat),
                    BIC = grp$Kstat,
                    is.best = 1:length(grp$Kstat) == ideal.k)
ideal.k.label = paste0("K = ", ideal.k)
p = ggplot(bic.df, aes(K, BIC)) +
    geom_point() +
    geom_line() +
    geom_point(data=subset(bic.df, bic.df$is.best), aes(K, BIC), col='red', cex=2) +
    geom_text(data=subset(bic.df, bic.df$is.best), aes(K, BIC, label=ideal.k.label),
        col='red', cex=3, hjust=-0.75) +
    theme_bw()

ggsave(p, file = paste0(in.prefix, ".BIC.pdf"))

# --- Compare actual populations to clusters

clust.tbl = table(pop(gen), grp$grp)
clust.tbl = clust.tbl[,order(-colSums(clust.tbl))]
clust.tbl

clust.tbl.is = clust.tbl[is.sites,]
clust.tbl.ml = clust.tbl[ml.sites,]

# Make rows with just zeros
clust.tbl.is.0 = clust.tbl.is
clust.tbl.ml.0 = clust.tbl.ml
clust.tbl.is.0[clust.tbl.is.0 != 0] = 0
clust.tbl.ml.0[clust.tbl.ml.0 != 0] = 0

clust.tbl.is = rbind(clust.tbl.is, clust.tbl.ml.0)
clust.tbl.ml = rbind(clust.tbl.is.0, clust.tbl.ml)

library(adegraphics)

pdf(paste0(in.prefix, ".clusters.by.site.pdf"), height=6, width=6)

    table.value(clust.tbl.is,
        labelsx=paste("Cluster", 1:ideal.k),
        method="size",
        breaks=c(0,4,8,12,16,20),
        col=c("#00B6EB", "#00B6EB"))

    table.value(clust.tbl.ml,
        labelsx=paste("Cluster", 1:ideal.k),
        method="size",
        breaks=c(0,4,8,12,16,20),
        col=c("#F8766D", "#F8766D"),
        add=TRUE)

dev.off()

# ----------------------------------------------------------------------------------------
# --- Do DAPC
# ----------------------------------------------------------------------------------------

dapc1 = dapc(gen, grp$grp, n.pca=100, n.da=100, glPca=pca)

pdf(paste0(in.prefix, ".dapc.scatter.pdf"), height=6, width=6)

    scatter(dapc1, scree.da=FALSE, bg="white",
        pch=20, cell=0, cstar=0, col = brewer.pal(ideal.k,"Set3"), solid=0.4,
        cex=3, clab=0, leg=TRUE, txt.leg=paste("Cluster",1:ideal.k))

dev.off()

save.image(file=paste0(in.prefix, ".DAPCcheckpoint.Rdata"))

# ----------------------------------------------------------------------------------------
# --- Plot PCA loadings
# ----------------------------------------------------------------------------------------

#   pdf(paste0(in.prefix, ".pca.loadings.pdf"), height=12, width=20)
#
#       #loadingplot(abs(pca$loadings), threshold=quantile(abs(pca$loadings), 0.9995), axis=1)
#
#       contrib.1 = loadingplot(dapc1$var.contr, axis = 1,
#           thres = 0.07, lab.jitter = 1)
#
#       #contrib.2 = loadingplot(dapc1$var.contr, axis = 2,
#       #    thres = 0.07, lab.jitter = 1)
#
#   dev.off()
