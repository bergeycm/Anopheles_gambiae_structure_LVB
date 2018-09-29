#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Get geographic distances between sites and correlate with Fst
# ----------------------------------------------------------------------------------------

library(ggplot2)
library(ade4)
library(plyr)

is.sites = c("BANDA", "BUGALA", "BUGALAIS", "BUKASA", "NSADZI", "SSERINYA")
ml.sites = c("BUGALAML", "BUWAMA", "KAZZI", "KIYINDI", "MITYANA")

# --- Function to calculate distance in kilometers between two points. From:
# https://conservationecology.wordpress.com/2013/06/30/distance-between-two-points-in-r/

earth.dist = function (long1, lat1, long2, lat2) {
    rad  = pi / 180
    a1   = lat1  * rad
    a2   = long1 * rad
    b1   = lat2  * rad
    b2   = long2 * rad
    dlon = b2 - a2
    dlat = b1 - a1
    a    = (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
    c    = 2 * atan2(sqrt(a), sqrt(1 - a))
    R    = 6378.145
    d    = R * c
    return(d)
}

# --- Read in GPS points

gps = read.table("data/site_gps_coordinates.txt")
gps = gps[gps$V1 != "Entebbe",]

pop.pairs = combn(gps$V1, 2)

# --- Get geographic distance matrix

dists = apply(pop.pairs, 2, function (pair) {

    lat.A = gps[gps$V1 == pair[1],]$V2
    lat.B = gps[gps$V1 == pair[2],]$V2
    lon.A = gps[gps$V1 == pair[1],]$V3
    lon.B = gps[gps$V1 == pair[2],]$V3

    this.dist = earth.dist (lon.A, lat.A, lon.B, lat.B)

    return(list(site.A = pair[1], site.B = pair[2], dist=this.dist))
})

dists.df = data.frame(do.call(rbind, dists))
dists.df$dist   = unlist(dists.df$dist)
dists.df$site.A = unlist(dists.df$site.A)
dists.df$site.B = unlist(dists.df$site.B)

# ----------------------------------------------------------------------------------------

# --- Read in FST for all SNPs

fst = read.table("data/all.pass.snp.flt.eigen.fst.se.out.fst.txt", row.names=NULL)[,-c(1)]
fst = fst / 1000

fst.df = cbind(expand.grid(names(fst), names(fst)), do.call(c, c(matrix(fst))))

names(fst.df) = c("site.A", "site.B", "fst")

# Let Mweena (site on Bugala) represent island-like Bugala
# and Bugoma (site on Bugala) represent mainland-like Bugala
dists.df$site.A[dists.df$site.A == "Bugala_Mweena"] = "BUGALAIS"
dists.df$site.B[dists.df$site.B == "Bugala_Mweena"] = "BUGALAIS"
dists.df$site.A[dists.df$site.A == "Bugala_Bugoma"] = "BUGALAML"
dists.df$site.B[dists.df$site.B == "Bugala_Bugoma"] = "BUGALAML"

dists.df$site.A = toupper(dists.df$site.A)
dists.df$site.B = toupper(dists.df$site.B)

# ----------------------------------------------------------------------------------------

# --- Merge distance and Fst data.frames

both = merge(dists.df, fst.df)

both$island.ct = both$site.A %in% is.sites + both$site.B %in% is.sites

both$island.ct = factor(both$island.ct)
levels(both$island.ct) = list("Mainland-Mainland" = "0",
                              "Mainland-Island"   = "1",
                              "Island-Island"     = "2")

# --- Plot FST against geographic distance

# Neat solution for putting lm equations on facets stolen from:
# https://stackoverflow.com/a/19701636

lm_eqn = function(df){
    m = lm(fst ~ dist, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2*', p ='~p,
         list(a = as.numeric(format(coef(m)[1], digits = 2)),
              b = as.numeric(format(coef(m)[2], digits = 2)),
             r2 = format(summary(m)$r.squared,  digits = 3),
              p = format(summary(m)$coefficients[2,4], digits=4)
))
    as.character(as.expression(eq));
}

eq = ddply(both,.(island.ct),lm_eqn)

p = ggplot(both, aes(dist, fst, col=island.ct)) +
    geom_smooth(aes(col=island.ct), method = "lm", se=FALSE,
        lwd=1, linetype=2, formula = y ~ x) +
    geom_point() +
    geom_text(data=eq, aes(x = 125, y = 0.0275, label=V1),
        parse = TRUE, inherit.aes=FALSE) +
    xlab("Distance between sites (km)") +
    facet_grid(factor(island.ct) ~ .) +
    ylab(expression(F[ST])) +
    scale_color_discrete(name="Comparison Type:",
                         breaks=c("Mainland-Mainland",
                                  "Mainland-Island",
                                  "Island-Island"),
                         labels=c("Mainland-Mainland",
                                  "Mainland-Island",
                                  "Island-Island")) +
    theme_bw()

ggsave(p, file="reports/fst_vs_distance.pdf", width=12, height=8)

# --- Do Mantel test and linear regression

both.sp = split(both, both$island.ct)

do.mantel = function (df) {

    matrix.dim = 1 + length(unique(df$site.A))

    dist.mat = matrix(NA, ncol=matrix.dim, nrow=matrix.dim)
    dist.mat[lower.tri(dist.mat)] = df$dist
    row.names(dist.mat) = colnames(dist.mat) = unique(c(df$site.A, df$site.B))

    fst.mat = matrix(NA, ncol=matrix.dim, nrow=matrix.dim)
    fst.mat[lower.tri(fst.mat)] = df$fst
    row.names(fst.mat) = colnames(fst.mat) = unique(c(df$site.A, df$site.B))

    mantel.rtest(as.dist(dist.mat), as.dist(fst.mat), nrepet = 9999)
}

sink("reports/mantel_results.txt")
    do.mantel(both)
    fit = lm(fst ~ dist, data=both)
    summary(fit)
sink()

sink("reports/mantel_results.Mainland-Mainland.txt")
    do.mantel(both.sp[["Mainland-Mainland"]])
    fit = lm(fst ~ dist, data=both.sp[["Mainland-Mainland"]])
    summary(fit)
sink()

sink("reports/mantel_results.Island-Island.txt")
    do.mantel(both.sp[["Island-Island"]])
    fit = lm(fst ~ dist, data=both.sp[["Island-Island"]])
    summary(fit)
sink()

sink("reports/mantel_results.Mainland-Island.txt")
    # Mantel impossible
    fit = lm(fst ~ dist, data=both.sp[["Mainland-Island"]])
    summary(fit)
sink()
