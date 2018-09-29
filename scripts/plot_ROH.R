#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Plot ROH
# ----------------------------------------------------------------------------------------

library(ggplot2)

# Bring in file with number of ROH segments, NSEG, by individual
roh.indiv = read.table("results/ROH/all.pass.snp.flt.noinv.hom.indiv", header=TRUE)

# Bring in info on all ROH segments to check length
hom = read.table("results/ROH/all.pass.snp.flt.noinv.hom", header=TRUE)

# Assert that no small ROH regions are included.
# (This is all hom is used for.)
if (min(hom$KB) < 100) {
    stop("ERROR: Summary file includes small (<100 kb) ROH fragments. Check PLINK call.")
}

# FROH = length of ROHs (only those > 100kb) divided by genome size

# AgamP4 Golden Path Length. Contains X-chromosome
genome.len = 273109044

roh.indiv$FROH = (roh.indiv$KB * 1000) / genome.len

# ----------------------------------------------------------------------------------------

# Bring in individual info

ind.info = read.csv("data/ssese_individual_info_bugala_split.csv")

# Bring in plate info to infer sequencing ID (e.g. ssese42.PE)
plates = read.table("data/sample_to_seq_id_mapping.txt")
names(plates) = c("mosquito.id", "seq.id")

# Reduce individual info to just sequenced animals
ind.info = ind.info[ind.info$mosquito_id %in% plates$mosquito.id,]

roh = merge(roh.indiv, plates, by.x="FID", by.y="seq.id")

island   = as.character(roh$FID)
site     = as.character(roh$FID)
sex      = as.character(roh$FID)

for (i in 1:dim(ind.info)[1]) {

    this.ind = ind.info$mosquito_id[i]
    this.ind.info = ind.info[ind.info$mosquito_id == this.ind,]

    this.seq.id = plates$seq.id[plates$mosquito.id == this.ind]

    island = replace(island, island==this.seq.id, this.ind.info$island)
    site   = replace(site,   site==this.seq.id,   this.ind.info$site_name)
    sex    = replace(sex,    sex==this.seq.id,    this.ind.info$sex)
}

roh$island = island
roh$site   = site
roh$sex    = sex

# Remove unsplit Bugala
roh = roh[roh$island != "BUGALA",]

islands = c("BANDA", "BUGALA", "BUGALAIS", "BUKASA", "NSADZI", "SSERINYA")
mainland = sort(unique(roh$island[!roh$island %in% islands]))

roh$is.island = FALSE
roh[roh$island %in% islands,]$is.island = TRUE

roh[roh$island == "BUGALAIS",]$island = "BUGALA (I)"
roh[roh$island == "BUGALAML",]$island = "BUGALA (M)"
islands[islands == "BUGALAIS"] = "BUGALA (I)"
mainland[mainland == "BUGALAML"] = "BUGALA (M)"

roh[roh$island == "KAZZI",]$island = "KAAZI"
roh[roh$island == "MITYANA",]$island = "WAMALA"
mainland[mainland == "KAZZI"] = "KAAZI"
mainland[mainland == "MITYANA"] = "WAMALA"

# ----------------------------------------------------------------------------------------

p = ggplot(roh, aes(FROH, NSEG, color=is.island, pch=(island=="BANDA"))) +
    geom_point() +
    theme_bw(base_size=15) +
    scale_color_discrete(guide=FALSE) +
    scale_shape_discrete(guide=FALSE) +
    xlab(expression("F"["ROH (100kb)"])) +
    ylab("Number of ROH")

ggsave(p, filename="reports/roh_plot.pdf", height=7, width=7)

# ----------------------------------------------------------------------------------------

simple.cap = function(x) {
    s = strsplit(x, " ")[[1]]
    lstr = paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
        sep="", collapse=" ")

    gsub("\\(i\\)", "(I)", gsub("\\(m\\)", "(M)", lstr))
}

roh$island = as.vector(sapply(roh$island, simple.cap))
islands = as.vector(sapply(islands, simple.cap))
mainland = as.vector(sapply(mainland, simple.cap))

p = ggplot(roh, aes(factor(island, levels=rev(c(islands, mainland))), FROH,
        fill=is.island)) +
    geom_boxplot(outlier.size=0.5, outlier.shape = NA) +
    xlab("") +
    ylab(expression(atop("F"["ROH"],""))) +
    scale_fill_discrete(guide=FALSE) +
    coord_flip() +
    theme_bw(base_size=15) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black"))

ggsave(p, filename="reports/roh_boxplot.pdf", height=3, width=4)
