#!/usr/bin/env Rscript

# ========================================================================================
# --- Make sampling map and plot of haplotypes
# ========================================================================================

library(RgoogleMaps)
library(mapplots)
library(seqinr)
library(plotrix)
library(RColorBrewer)

options(stringsAsFactors = FALSE)

site.coords = read.table("data/site_gps_coordinates.txt", header=FALSE, sep="\t")
names(site.coords) = c("site", "lat", "lon")

lat.width  = max(site.coords$lat) - min(site.coords$lat)
lon.height = max(site.coords$lon) - min(site.coords$lon)

map.lat.ylim = range(site.coords$lat)
map.lon.xlim = range(site.coords$lon)

center = c(mean(map.lat.ylim), mean(map.lon.xlim))
zoom = 10

# Fix names
site.coords$site[site.coords$site == "Kazzi"] = "Kaazi"
site.coords$site[site.coords$site == "Mityana"] = "Wamala"

# ----------------------------------------------------------------------------------------
# --- Add info to points
# ----------------------------------------------------------------------------------------

site.coords$size = "small"
site.coords$col = "red"
site.coords$char = ""

markers = cbind.data.frame(site.coords$lat, site.coords$lon,
                           site.coords$size, site.coords$col, site.coords$char)

names(markers) = c("lat", "lon", "size", "col", "char")

# ----------------------------------------------------------------------------------------
# --- Plot sattelite map showing sites
# ----------------------------------------------------------------------------------------

ss.map = GetMap.bbox(lonR = range(lon), latR = range(lat),
                     center = center,
                     zoom = zoom,
                     markers = markers,
                     maptype = "terrain",     # or hybrid
                     destfile = "results/ssese_sites.png",
                     SCALE=10)

# ----------------------------------------------------------------------------------------
# --- Read in individual information, add geographic coordinates
# ----------------------------------------------------------------------------------------

ind.info = read.table("data/ssese_individual_info_simple.txt")
names(ind.info) = c("samp.id", "seq.id", "island", "site", "sex", "species")
ind.info = ind.info[ind.info$species == "AG",]

# Expand site for Bugala to include village
ind.info[ind.info$island == "BUGALA",]$island =
    paste(ind.info$island, ind.info$site, sep="_")[ind.info$island == "BUGALA"]

# Fix names
ind.info[ind.info$island == "KAZZI",  ]$island = "KAAZI"
ind.info[ind.info$island == "MITYANA",]$island = "WAMALA"

ind.group.index = match(ind.info$island, toupper(site.coords$site))
ind.info.coords = data.frame(ind.info,
                             lat = site.coords$lat[ind.group.index],
                             lon = site.coords$lon[ind.group.index])

# Bring in haplotypes of each individual for sweeps on 2L and X
haps.2L = read.table(paste0("results/for_dendrograms/",
    "chr2L.2L_34Mb.1e+05.haploid.with-ag1000g.sans-snpeff.haps.haplotypes.txt"))
haps.X = read.table(paste0("results/for_dendrograms/",
    "chrX.X_9Mb.1e+05.haploid.with-ag1000g.sans-snpeff.haps.haplotypes.txt"))

names(haps.2L) = names(haps.X) = c("hap", "ind", "hap.index", "pop")

# Fix names
haps.2L[haps.2L$pop == "KAZZI",  ]$pop = "KAAZI"
haps.2L[haps.2L$pop == "MITYANA",]$pop = "WAMALA"
haps.X[haps.X$pop == "KAZZI",  ]$pop = "KAAZI"
haps.X[haps.X$pop == "MITYANA",]$pop = "WAMALA"

haps.both = merge(haps.2L, haps.X, by=c("hap", "ind", "pop"),
    suffixes=c(".2L_34Mb", ".X_9Mb"))

ind.info.coords.haps = merge(ind.info.coords, haps.both, by.x="seq.id", by.y="ind")

# Set boring haps to 0
# 2L
hap.names.2L = names(table(ind.info.coords.haps$hap.index.2L_34Mb))
boring.haps.2L = hap.names.2L[table(ind.info.coords.haps$hap.index.2L_34Mb) == 1]
boring.haps.2L.bool = ind.info.coords.haps$hap.index.2L_34Mb %in% boring.haps.2L
ind.info.coords.haps[boring.haps.2L.bool,]$hap.index.2L_34Mb = 0
# X
hap.names.X = names(table(ind.info.coords.haps$hap.index.X_9Mb))
boring.haps.X = hap.names.X[table(ind.info.coords.haps$hap.index.X_9Mb) == 1]
boring.haps.X.bool = ind.info.coords.haps$hap.index.X_9Mb %in% boring.haps.X
ind.info.coords.haps[boring.haps.X.bool,]$hap.index.X_9Mb = 0

# Number of individuals a row represents, which is 1.
ind.info.coords.haps$fake.count = 1

# Factor
ind.info.coords.haps$hap.factor.2L_34Mb =
    as.factor(as.character((ind.info.coords.haps$hap.index.2L_34Mb)))
ind.info.coords.haps$hap.factor.X_9Mb   =
    as.factor(as.character((ind.info.coords.haps$hap.index.X_9Mb)))

# Make things that should be factors, factors

# ----------------------------------------------------------------------------------------
# --- Get all data needed to plot pie charts (coords and taxon counts) together
# ----------------------------------------------------------------------------------------

xyz.2L = make.xyz(x=ind.info.coords.haps$lon,
                  y=ind.info.coords.haps$lat,
                  z=ind.info.coords.haps$fake.count,
                  group=ind.info.coords.haps$hap.factor.2L_34Mb)

xyz.X  = make.xyz(x=ind.info.coords.haps$lon,
                  y=ind.info.coords.haps$lat,
                  z=ind.info.coords.haps$fake.count,
                  group=ind.info.coords.haps$hap.factor.X_9Mb)

# ----------------------------------------------------------------------------------------
# --- Write table of haplotypes by population
# ----------------------------------------------------------------------------------------

# Get count of haplotypes - LVB only
hap.count.2L = table(ind.info.coords.haps[,c("pop", "hap.factor.2L_34Mb")])
hap.count.X  = table(ind.info.coords.haps[,c("pop", "hap.factor.X_9Mb")])

# Get count of haplotypes - with Ag1000G
# Set boring haps to 0
hap.names.2L = names(table(haps.both$hap.index.2L_34Mb))
boring.haps.2L = hap.names.2L[table(haps.both$hap.index.2L_34Mb) <= 3]
boring.haps.2L.bool = haps.both$hap.index.2L_34Mb %in% boring.haps.2L
haps.both[boring.haps.2L.bool,]$hap.index.2L_34Mb = 0

hap.names.X = names(table(haps.both$hap.index.X_9Mb))
boring.haps.X = hap.names.X[table(haps.both$hap.index.X_9Mb) <= 3]
boring.haps.X.bool = haps.both$hap.index.X_9Mb %in% boring.haps.X
haps.both[boring.haps.X.bool,]$hap.index.X_9Mb = 0

hap.count.2L.full = table(haps.both[,c("pop", "hap.index.2L_34Mb")])
hap.count.X.full  = table(haps.both[,c("pop", "hap.index.X_9Mb")])

write.table(hap.count.2L.full, file="reports/haplotype_count_by_pop.2L_34Mb.txt")
write.table(hap.count.X.full,  file="reports/haplotype_count_by_pop.X_9Mb.txt")

# ----------------------------------------------------------------------------------------
# --- Make map of sampling localities
# ----------------------------------------------------------------------------------------

# Correct xyz$z object containing coordinates and counts
# Use 0 rather than NA

xyz.2L$z[is.na(xyz.2L$z)] = 0
xyz.X$z[is.na(xyz.X$z)] = 0

zoom = 9

png (filename = "results/ssese_sampling_sites.png", width=1280, height=1280)

# Make Google map without markers
ssese.map.nomrk = GetMap.bbox(
    lonR = range(lon), latR = range(lat),
    center = center,
    zoom = zoom,
    maptype = "terrain",
    SCALE=2,
    destfile="results/ssese_sites_no_markers.png",
    path=paste0("&style=feature:road|visibility:off",
                "&style=feature:poi|visibility:off",
                "&style=feature:administrative|visibility:off",
                "&style=feature:landscape|element:labels|visibility:off",
                "&style=feature:transit|visibility:off",
                "&style=feature:water|element:labels|visibility:off",
                "&style=feature:water|element:geometry.fill|",
                "visibility:on|hue:0x006a6f|saturation:-55|lightness:-65")
)

# Get count of individuals per site
site.coords$site.upper = toupper(site.coords$site)
ind.counts = merge(data.frame(table(ind.info.coords$island)), site.coords,
    by.x="Var1", by.y="site.upper")

# Remove extraneous Bugala labels
ind.counts$to.print = ind.counts$site
ind.counts$to.print[ind.counts$to.print == "Bugala_Bugoma"] = ""
ind.counts$to.print[ind.counts$to.print == "Bugala_Mweena"] = ""
ind.counts$to.print[ind.counts$to.print == "Bugala_Lutoboka"] = "Bugala"

ind.counts$type = "Mainland"
islands = c("Banda", "Bugala", "Bukasa", "Nsadzi", "Sserinya", "")
ind.counts[ind.counts$to.print %in% islands,]$type = "Island"

# Add points scaled by site
scale.factor = 4
PlotOnStaticMap(ssese.map.nomrk,
                lat = ind.counts$lat, lon = ind.counts$lon,
                cex = (ind.counts$Freq / scale.factor) + 1,
                pch=16,
                col="black",
                verbose=0)

PlotOnStaticMap(ssese.map.nomrk,
                lat = ind.counts$lat, lon = ind.counts$lon,
                cex = ind.counts$Freq / scale.factor,
                pch=16,
                col=c("#00BFC4", "#F8766D")[factor(ind.counts$type)],
                verbose=0,
                add=TRUE)

# Up and right for all except:
# Sserinya - down right
# Banda    -      right
# Bukasa   - down right
# Bugala   - down left
# Kiyindi  - up   left
ind.counts$lat.adj = 0.04
ind.counts[ind.counts$to.print == "Sserinya",]$lat.adj = -0.05
ind.counts[ind.counts$to.print == "Bukasa",]$lat.adj   = -0.05
ind.counts[ind.counts$to.print == "Bugala",]$lat.adj   = -0.05
ind.counts$lon.adj = 0.05
ind.counts[ind.counts$to.print == "Sserinya",]$lon.adj = 0.08
ind.counts[ind.counts$to.print == "Banda",]$lon.adj    = 0.08
ind.counts[ind.counts$to.print == "Bugala",]$lon.adj   = -0.15
ind.counts[ind.counts$to.print == "Kiyindi",]$lon.adj  = -0.08

TextOnStaticMap(ssese.map.nomrk,
                lat = -0.0025 + ind.counts$lat + ind.counts$lat.adj,
                lon = -0.0025 + ind.counts$lon + ind.counts$lon.adj,
                labels = ind.counts$to.print,
                col = c("white", "black")[1 + (ind.counts$to.print %in% islands)],
                cex = 3,
                add = TRUE)
TextOnStaticMap(ssese.map.nomrk,
                lat = ind.counts$lat + ind.counts$lat.adj,
                lon = ind.counts$lon + ind.counts$lon.adj,
                labels = ind.counts$to.print,
                col = c("black", "white")[1 + (ind.counts$to.print %in% islands)],
                cex = 3,
                add = TRUE)

# Add scale bar
# http://www.edwilliams.org/gccalc.htm
lat0 = -0.7855117210902682
lon0 = 32.9180375726562
lat1 = -0.785505
lon1 = 32.69344

PlotArrowsOnStaticMap(ssese.map.nomrk,
                      lat0 = lat0,
                      lon0 = lon0,
                      lat1 = lat1,
                      lon1 = lon1,
                      FUN = segments,
                      add = TRUE,
                      col = '#DDE3D5',
                      lwd = 3)
TextOnStaticMap(ssese.map.nomrk,
                lat = lat0 - 0.025,
                lon = (lon0 + lon1) / 2,
                labels = "25 km",
                col = '#DDE3D5',
                cex = 2,
                add = TRUE)

# Add labels for Lake Victoria and Uganda
TextOnStaticMap(ssese.map.nomrk,
                lat = -0.45,
                lon = 33,
                labels = expression(italic('Lake Victoria')),
                col = '#DDE3D5',
                cex = 5,
                add = TRUE)

TextOnStaticMap(ssese.map.nomrk,
                lat = 0.6,
                lon = center[2],
                labels = "UGANDA",
                col = 'black',
                cex = 6,
                add = TRUE)

legend.pts = seq(from=5, to=20, by=5)
legend("topleft",
    inset = c(0.05, 0.025),
    legend = legend.pts,
    pch = 21,
    cex = 2,
    pt.cex = (legend.pts / scale.factor) + 1,
    bg = "#DDE3D5",
    box.col = "#006a6f",
    box.lwd = 2,
    text.col = "black",
    col = "#006a6f",
    pt.bg = "black",
    horiz = TRUE,
    title = "Sample Count:")

dev.off()

# ----------------------------------------------------------------------------------------
# --- Get ratio of haplotype counts between island and mainland (for X sweep only)
# ----------------------------------------------------------------------------------------

islands = c(islands, "BUGALAIS")

ind.info.coords.haps.isl  =
    ind.info.coords.haps[ind.info.coords.haps$pop %in% toupper(islands),]
ind.info.coords.haps.main =
    ind.info.coords.haps[!(ind.info.coords.haps$pop %in% toupper(islands)),]

x.cols = c("island","hap.index.X_9Mb")

hap.count.isl.X  = colSums(table(ind.info.coords.haps.isl[,x.cols]))
hap.count.main.X = colSums(table(ind.info.coords.haps.main[,x.cols]))

hap.freq.isl.X  = hap.count.isl.X  / sum(hap.count.isl.X)
hap.freq.main.X = hap.count.main.X / sum(hap.count.main.X)

hap.freq.isl.X  = data.frame(hap.num   = names(hap.freq.isl.X),
                             hap.count = hap.count.isl.X,
                             hap.freq  = hap.freq.isl.X)
hap.freq.main.X = data.frame(hap.num   = names(hap.freq.main.X),
                             hap.count = hap.count.main.X,
                             hap.freq  = hap.freq.main.X)

hap.freq.both.X = merge(hap.freq.isl.X, hap.freq.main.X,
    by="hap.num",
    all=TRUE, suffixes=c(".isl", ".main"))

# Set NA counts to zero
hap.freq.both.X$hap.count.isl [is.na(hap.freq.both.X$hap.count.isl) ] = 0
hap.freq.both.X$hap.freq.isl  [is.na(hap.freq.both.X$hap.freq.isl)  ] = 0
hap.freq.both.X$hap.count.main[is.na(hap.freq.both.X$hap.count.main)] = 0
hap.freq.both.X$hap.freq.main [is.na(hap.freq.both.X$hap.freq.main) ] = 0

# Exclude misc. haplotypes
hap.freq.both.X[hap.freq.both.X$hap.num == "0",]$hap.freq.isl  = NA
hap.freq.both.X[hap.freq.both.X$hap.num == "0",]$hap.freq.main = NA

hap.freq.both.X$isl.to.main.ratio =
    hap.freq.both.X$hap.freq.isl / (hap.freq.both.X$hap.freq.main + 1e-256)

island.specific.hap =
    hap.freq.both.X$hap.num[which.max(hap.freq.both.X$isl.to.main.ratio)]

# Write table with island-specific hap counts
island.specific.hap.freq = hap.count.X.full[,island.specific.hap]

write.table(t(t(island.specific.hap.freq)),
    file="reports/haplotype_count_by_pop.X_9Mb.island-specific-hap.txt",
    col.names=FALSE)

write.table(hap.freq.both.X, "reports/haplotype_count_by_group.X_9Mb.txt",
    sep="\t", quote=FALSE, row.names=FALSE)

# Also write table of individuals with island-specific hap
write.table(haps.X[haps.X$hap.index == island.specific.hap,],
    file="reports/individuals_with_island_specific_X_haplotype.txt",
    quote=FALSE, row.names=FALSE)

# ----------------------------------------------------------------------------------------
# --- Get ratio of haplotype counts between Uganda and rest of Africa (for 2L sweep only)
# ----------------------------------------------------------------------------------------

pop.names = row.names(hap.count.2L.full)
ugandans = c("Uganda", pop.names[pop.names == toupper(pop.names)])

hap.count.2L.full.ug  =
    hap.count.2L.full[row.names(hap.count.2L.full) %in% ugandans,]
hap.count.2L.full.rest =
    hap.count.2L.full[!(row.names(hap.count.2L.full) %in% ugandans),]

hap.count.ug.2L   = colSums(hap.count.2L.full.ug)
hap.count.rest.2L = colSums(hap.count.2L.full.rest)

hap.freq.ug.2L   = hap.count.ug.2L / sum(hap.count.ug.2L)
hap.freq.rest.2L = hap.count.rest.2L / sum(hap.count.rest.2L)

hap.freq.ug.2L   = data.frame(hap.num   = names(hap.freq.ug.2L),
                              hap.count = hap.count.ug.2L,
                              hap.freq  = hap.freq.ug.2L)
hap.freq.rest.2L = data.frame(hap.num   = names(hap.freq.rest.2L),
                              hap.count = hap.count.rest.2L,
                              hap.freq  = hap.freq.rest.2L)

hap.freq.both.2L = merge(hap.freq.ug.2L, hap.freq.rest.2L,
    by="hap.num",
    all=TRUE, suffixes=c(".ug", ".rest"))

# Set NA counts to zero
hap.freq.both.2L$hap.count.ug  [is.na(hap.freq.both.2L$hap.count.ug)  ] = 0
hap.freq.both.2L$hap.freq.ug   [is.na(hap.freq.both.2L$hap.freq.ug)   ] = 0
hap.freq.both.2L$hap.count.rest[is.na(hap.freq.both.2L$hap.count.rest)] = 0
hap.freq.both.2L$hap.freq.rest [is.na(hap.freq.both.2L$hap.freq.rest) ] = 0

# Exclude misc. haplotypes
hap.freq.both.2L[hap.freq.both.2L$hap.num == "0",]$hap.freq.ug   = NA
hap.freq.both.2L[hap.freq.both.2L$hap.num == "0",]$hap.freq.rest = NA

# Add tiny amount to avoid dividing by zero
hap.freq.both.2L$ug.to.rest.ratio =
    hap.freq.both.2L$hap.freq.ug / (hap.freq.both.2L$hap.freq.rest + 1e-256)

# There are several, but this is the highest in Uganda
hap.freq.both.2L.notzero = hap.freq.both.2L[hap.freq.both.2L$hap.num != 0,]
uganda.specific.hap =
    hap.freq.both.2L.notzero$hap.num[which.max(hap.freq.both.2L.notzero$hap.count.ug)]

# Write table with uganda-specific hap counts
uganda.specific.hap.freq = hap.count.2L.full[,uganda.specific.hap]

write.table(t(t(uganda.specific.hap.freq)),
    file="reports/haplotype_count_by_pop.2L_34Mb.uganda-specific-hap.txt",
    col.names=FALSE)

write.table(hap.freq.both.2L, "reports/haplotype_count_by_group.2L_34Mb.txt",
    sep="\t", quote=FALSE, row.names=FALSE)

# Also write table of individuals with island-specific hap
write.table(haps.2L[haps.2L$hap.index == uganda.specific.hap,],
    file="reports/individuals_with_uganda_specific_2L_haplotype.txt",
    sep="\t", quote=FALSE, row.names=FALSE)

# ----------------------------------------------------------------------------------------
# --- Plot pi charts of haplotypes
# ----------------------------------------------------------------------------------------

# Make Google map without markers (similar to above, but colors tweaked)
ssese.map.nomrk.blk = GetMap.bbox(
    lonR = range(lon), latR = range(lat),
    center = center,
    zoom = zoom,
    maptype = "terrain",
    SCALE=2,
    destfile="results/ssese_sites_no_markers_black.png",
    path=paste0("&style=feature:road|visibility:off",
                "&style=feature:poi|visibility:off",
                "&style=feature:administrative|visibility:off",
                "&style=feature:landscape|element:labels|visibility:off",
                "&style=feature:transit|visibility:off",
                "&style=feature:water|element:labels|visibility:off",
                "&style=feature:water|element:geometry.fill|",
                "visibility:on|color:0x5d6061|saturation:0|lightness:0")
)

# Put input data into lists to parallelize below code
xyz.l = list()
xyz.l[["2L"]] = xyz.2L
xyz.l[["X"]]  = xyz.X

for (chr in c("2L", "X")) {

    if (chr == "2L") {
        hap.out.code = "2L_34Mb"
    } else {
        hap.out.code = "X_9Mb"
    }

    png (filename = paste0("results/ssese_haplotypes_sweep_", hap.out.code, ".png"),
        width=1280, height=1280)

    # Display Google map from before on screen
    stat.map = PlotOnStaticMap(ssese.map.nomrk.blk)

    # Set colors and make them semi-transparent
    #col = c(brewer.pal(9, "Set1"), "#42f4cb")
    num.cols = ncol(xyz.l[[chr]]$z) - 1
    col = c("#ffffff", brewer.pal(num.cols, "Set1"))

    for (i in 1:length(col)) {
        col[i] = col2alpha(color=col[i], alpha = 0.75)
    }

    # Get individual counts for sizing

    total.ind = as.vector(rowSums(xyz.l[[chr]]$z, na.rm=TRUE))

    max.ind = max(total.ind)
    min.ind = min(total.ind)

    max.size = 1500
    min.size = (min.ind / max.ind) * max.size

    pie.sizes = ((max.size - min.size) *
        (total.ind - min.ind)) / (max.ind - min.ind)
    pie.sizes = pie.sizes + min.size

    # Convert areas to radius so as not to misrepresent differences.
    pie.sizes = sqrt(pie.sizes / pi)

    # Get xpos and ypos coordintates relative to map, not in degrees
    # with RgoogleMaps' LatLon2XY.centered

    # Adjust Banda longitude
    banda.lon = site.coords[site.coords$site == "Banda",]$lon
    chr.adj = paste0(chr, ".adj")
    xyz.l[[chr.adj]] = xyz.l[[chr]]
    xyz.l[[chr.adj]]$x[xyz.l[[chr.adj]]$x == banda.lon] = banda.lon + 0.05

    # Adjust Bugala (Mweena) latitude
    mweena.lat = site.coords[site.coords$site == "Bugala_Mweena",]$lat
    chr.adj = paste0(chr, ".adj")
    xyz.l[[chr.adj]]$x[xyz.l[[chr.adj]]$x == mweena.lat] = mweena.lat - 0.05

    coord.convert = LatLon2XY.centered(ssese.map.nomrk.blk,
        lat = xyz.l[[chr.adj]]$y, lon = xyz.l[[chr.adj]]$x, zoom=zoom)

    if (chr == "X") {
        # Set island-speicific hap to black
        col[which(island.specific.hap == colnames(xyz.l[[chr]]$z))] = "black"
    } else {
        # Set Uganda-speicific hap to black
        col[which(uganda.specific.hap == colnames(xyz.l[[chr]]$z))] = "black"
    }

    for (i in 1:length(xyz.l[[chr.adj]]$x)) {

        # Only pass correct number of colors if not all group represented in pie chart
        # That is remove colors for missing pie slices for this location
        this.col = col[xyz.l[[chr.adj]]$z[i,] != 0]

        floating.pie(coord.convert$newX[i], coord.convert$newY[i], xyz.l[[chr.adj]]$z[i,],
            radius=pie.sizes[i], col=this.col, border="black")
    }

    # Add labels (Should place these more carefully, but fine for now.
    #text(coord.convert$newX, coord.convert$newY,
    #    site.coords[order(site.coords$lon),]$site, pos=c(3,3,2,3,3,3,2,3,3), offset=1.5)

    dev.off()
}

# ----------------------------------------------------------------------------------------
# --- Plot pi charts of ADMIXTURE ancestry
# ----------------------------------------------------------------------------------------

fam.file = paste0("data/ssese_with_ag1000g/",
     "chr3.pass.snp.phased.ag1000g.strict.thinned_for_ADMIXTURE.withM.fam")

adm.files = list.files(path="data/ssese_with_ag1000g/adm_subsets/",
    pattern="CLUMPAK_withM.*Q$",
    full.names=TRUE)

adm.inds = read.table(fam.file)$V1

ind.info = read.table("data/ssese_individual_info_simple_bugala_split.txt")

# Bring in plate info to infer sequencing ID (e.g. ssese42.PE)
plates = read.table("data/sample_to_seq_id_mapping.txt")
names(plates) = c("mosquito.id", "seq.id")

# Reduce individual info to just sequenced animals
ind.info = ind.info[ind.info$V1 %in% plates$mosquito.id,]

ind.info.ag1000g = read.table("data/ag1000g.phase1.ar3/samples.all.txt",
    sep="\t", quote="", header=TRUE)

for (adm.file in adm.files) {

    adm = read.table(adm.file)

    this.k = ncol(adm)

    pal = read.table(paste0("data/ssese_with_ag1000g/adm_subsets/CLUMPAK_withM.", this.k,
        ".colors.txt"), comment="'", header=TRUE)

    names(adm) = paste0("ADM_", 1:this.k)
    adm = cbind(adm.inds, adm)

    adm$pops = do.call(c, lapply(adm$adm.inds, function(x) {
        pop.ss = as.character(ind.info[ind.info$V2 == x, 3])
        pop.ag = paste(
            as.character(ind.info.ag1000g[ind.info.ag1000g$ox_code == x,]$country),
            as.character(ind.info.ag1000g[ind.info.ag1000g$ox_code == x,]$m_s),
            sep="_")
        if (length(pop.ss)) {
            return(pop.ss)
        } else if (length(pop.ag)) {
            return(pop.ag)
        } else if (grepl("AD", x)) {
            return("AD")
        } else {
            return(x)
        }
    }))

    # Average ancestry within each population
    grp.avg.anc = do.call(cbind, lapply(1:this.k, function (adm.comp) {
        aggregate(adm[[paste0("ADM_", adm.comp)]], by=list(adm$pops), FUN=mean)
    }))

    grp.avg.anc = grp.avg.anc[,c(1,2 * (1:this.k))]

    site.coords.simp = site.coords

    # Let Mweena represent island-like Bugala
    # and Bugoma represent mainland-like Bugala
    site.coords.simp$site.upper[site.coords.simp$site.upper == "BUGALA_MWEENA"] = "BUGALAIS"
    site.coords.simp$site.upper[site.coords.simp$site.upper == "BUGALA_BUGOMA"] = "BUGALAML"

    # Merge
    grp.avg.anc.coords = merge(site.coords.simp, grp.avg.anc,
        by.x="site.upper", by.y="Group.1")

    anc.components = grp.avg.anc.coords[,grepl("x", names(grp.avg.anc.coords))]

    row.names(anc.components) = paste(grp.avg.anc.coords$lon, grp.avg.anc.coords$lat,
        sep=", ")

    # Bring in number of individuals (divide number of chromosomes by 2)
    anc.comp.ind = merge(anc.components,
        cbind(names(rowSums(xyz.l[["X"]]$z)), rowSums(xyz.l[["X"]]$z) / 2),
        by='row.names', sort=FALSE)

    # Abort if order of pops gets messed up (which should never happen)
    if (sum(row.names(anc.components) != anc.comp.ind$Row.names) > 0) {
        cat("Order of sites disrupted during merge. Aborting.\n")
        quit()
    }

    ind.count = as.numeric(anc.comp.ind$V2)

    anc.components = as.matrix(
        grp.avg.anc.coords[,grepl("^x", names(grp.avg.anc.coords))]
    )

    row.names(anc.components) = paste(grp.avg.anc.coords$lon, grp.avg.anc.coords$lat,
        sep=", ")

    xyz.adm = list(x = grp.avg.anc.coords$lon,
                   y = grp.avg.anc.coords$lat,
                   z = anc.components)

    png (filename = paste0(adm.file, ".map.png"),
        width=1280, height=1280)

    # Display Google map from before on screen
    stat.map = PlotOnStaticMap(ssese.map.nomrk.blk)

    # Set colors and make them semi-transparent
    #col = c(brewer.pal(9, "Set1"), "#42f4cb")
    num.cols = ncol(xyz.adm$z) - 1
    col = pal$colors

    # Leaving alpha alone now
    # for (i in 1:length(col)) {
    #     col[i] = col2alpha(color=col[i], alpha = 0.75)
    # }

    # Get individual counts for sizing

    max.ind = max(ind.count)
    min.ind = min(ind.count)

    max.size = 1500
    min.size = (min.ind / max.ind) * max.size

    pie.sizes = ((max.size - min.size) *
        (ind.count - min.ind)) / (max.ind - min.ind)
    pie.sizes = pie.sizes + min.size

    # Convert areas to radius so as not to misrepresent differences.
    pie.sizes = sqrt(pie.sizes / pi)

    # Get xpos and ypos coordintates relative to map, not in degrees
    # with RgoogleMaps' LatLon2XY.centered

    # Adjust Banda longitude
    banda.lon = site.coords[site.coords$site == "Banda",]$lon
    chr.adj = paste0(chr, ".adj")
    xyz.adm.adj = xyz.adm
    xyz.adm.adj$x[xyz.adm.adj$x == banda.lon] = banda.lon + 0.05

    # Adjust Bugala (Mweena) latitude
    mweena.lat = site.coords[site.coords$site == "Bugala_Mweena",]$lat
    chr.adj = paste0(chr, ".adj")
    xyz.adm.adj$x[xyz.adm.adj$x == mweena.lat] = mweena.lat - 0.05

    coord.convert = LatLon2XY.centered(ssese.map.nomrk.blk,
        lat = xyz.adm.adj$y, lon = xyz.adm.adj$x, zoom=zoom)

    for (i in 1:length(xyz.adm.adj$x)) {

        # Only pass correct number of colors if not all group represented in pie chart
        # That is remove colors for missing pie slices for this location
        this.col = col[xyz.adm.adj$z[i,] != 0]

        floating.pie(coord.convert$newX[i], coord.convert$newY[i], xyz.adm.adj$z[i,],
            radius=pie.sizes[i], col=this.col, border="black")
    }

    # Add labels (Should place these more carefully, but fine for now.
    #text(coord.convert$newX, coord.convert$newY,
    #    site.coords[order(site.coords$lon),]$site, pos=c(3,3,2,3,3,3,2,3,3), offset=1.5)

    dev.off()
}
