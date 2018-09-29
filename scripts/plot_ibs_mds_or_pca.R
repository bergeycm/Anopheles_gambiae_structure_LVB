#!/usr/bin/Rscript

options(stringsAsFactors = FALSE)

library(ggplot2)
library(grid)
library(gridExtra)
library(lattice)
library(scales)

args = commandArgs(trailingOnly = TRUE)
# MDS or PCA input
mdspca.in = args[1] # e.g. "data/all.pass.snp.flt.mds" or
                    #      "data/chr3.pass.snp.flt.noinv.eigenvec"
if (length(args) == 1) {
    split.bugala = FALSE
} else {
    split.bugala = as.logical(args[2]) # Should Bugala be split? Default is FALSE
}

ind.info = read.csv("data/ssese_individual_info.csv")

# Bring in plate info to infer sequencing ID (e.g. ssese42.PE)
plates = read.table("data/sample_to_seq_id_mapping.txt")
names(plates) = c("mosquito.id", "seq.id")

# Reduce individual info to just sequenced animals
ind.info = ind.info[ind.info$mosquito_id %in% plates$mosquito.id,]

if (split.bugala) {
    # Add site to BUGALA and SSERINYA
    ind.info[ind.info$island == "BUGALA",]$island = paste(
        ind.info[ind.info$island == "BUGALA",]$island,
        ind.info[ind.info$island == "BUGALA",]$site_name, sep=" - ")
    ind.info[ind.info$island == "SSERINYA",]$island = paste(
        ind.info[ind.info$island == "SSERINYA",]$island,
        ind.info[ind.info$island == "SSERINYA",]$site_name, sep=" - ")
}

# Figure out if we're plotting MDS or PCA
is.mds = FALSE
if (grepl("mds", mdspca.in)) {
    is.mds = TRUE
}

if (is.mds) {
    mdspca = read.table(mdspca.in, header = T)
    mdspca = merge(mdspca, plates, by.x="FID", by.y="seq.id")
} else {
    mdspca = read.table(mdspca.in)
    # Name as C to be parallel with MDS column names
    names(mdspca) = c("IID", "FID", paste0("C", 1:(ncol(mdspca) - 2)))
    mdspca = merge(mdspca, plates, by.x="FID", by.y="seq.id")
}

island   = as.character(mdspca$FID)
site     = as.character(mdspca$FID)
sex      = as.character(mdspca$FID)

for (i in 1:dim(ind.info)[1]) {

    this.ind = ind.info$mosquito_id[i]
    this.ind.info = ind.info[ind.info$mosquito_id == this.ind,]

    this.seq.id = plates$seq.id[plates$mosquito.id == this.ind]

    island = replace(island, island==this.seq.id, this.ind.info$island)
    site   = replace(site,   site==this.seq.id,   this.ind.info$site_name)
    sex    = replace(sex,    sex==this.seq.id,    this.ind.info$sex)
}

mdspca$island = island
mdspca$site   = site
mdspca$sex    = sex

if (split.bugala == FALSE) {
    islands = c("BANDA", "BUGALA", "BUKASA", "NSADZI", "SSERINYA")
} else {
    islands = c("BANDA", "BUKASA",
        "BUGALA - BUGOMA", "BUGALA - LUTOBOKA", "BUGALA - MWEENA",
        "BUKASA", "NSADZI", "SSERINYA - BBOSA", "SSERINYA - KASISA")
}

mdspca$is.island = FALSE
mdspca[mdspca$island %in% islands,]$is.island = TRUE

if (is.mds) {
    out.file = gsub("results/(.*).mds$", "reports/\\1.ibs_mds_plot.pdf", mdspca.in)
    out.file = gsub("data/(.*).mds$", "reports/\\1.ibs_mds_plot.pdf", out.file)
} else {
    out.file = gsub("data/(.*).eigenvec$", "reports/\\1.pca_plot.pdf", mdspca.in)
    # Not actually necessary, since this is only for plotting Ssese-only PCA
    #out.file = gsub("/ssese_with_ag1000g", "", out.file)
}

if (split.bugala) {
    out.file = gsub(".pdf$", ".bugala.pdf", out.file)
}

# Fix names (KAAZI and WAMALA)
mdspca$island[mdspca$island == "KAZZI"]   = "KAAZI"
mdspca$island[mdspca$island == "MITYANA"] = "WAMALA"

# Reorder
if (split.bugala == FALSE) {
    mdspca$island = factor(mdspca$island,
                        levels=c("BANDA", "BUKASA",
                                 "BUGALA",
                                 "NSADZI", "SSERINYA",
                                 "BUWAMA", "KAAZI", "KIYINDI", "WAMALA"))
} else {
    mdspca$island = factor(mdspca$island,
                        levels=c("BANDA", "BUKASA",
                                 "BUGALA - BUGOMA", "BUGALA - LUTOBOKA", "BUGALA - MWEENA",
                                 "NSADZI", "SSERINYA - BBOSA", "SSERINYA - KASISA",
                                 "BUWAMA", "KAAZI", "KIYINDI", "WAMALA"))
}

if (split.bugala == FALSE) {
    legend.colors = c(rep("#00BFC4", 5), rep("#F8766D", 4))
} else {
    legend.colors = c("grey", "red", "green", "blue",
                      rep("grey", 6), "purple", "orange")
}

# Add proportion of variance explained, if PCA
prop.explained = rep("", (ncol(mdspca) - 7))
if (grepl("eigenvec$", mdspca.in)) {
    val.in = gsub("vec$", "val", mdspca.in)
    val = read.table(val.in)
    prop.explained = paste0(" (", percent(val$V1 / sum(val$V1)), ")")
}

if (is.mds) {
    component.name = "Component "
} else {
    component.name = "PC"
}

# Lowercase letters
simple.cap = function(x) {
    s = strsplit(x, " ")[[1]]
    lstr = paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
        sep="", collapse=" ")

    gsub("\\(i\\)", "(I)", gsub("\\(m\\)", "(M)", lstr))
}

levels(mdspca$island) = do.call(c, lapply(levels(mdspca$island), simple.cap))

# For Ssese-only MDS or PCA, only plot first two components

p = ggplot(mdspca, aes(C1, C2, col=island, pch=island))

if (split.bugala == FALSE) {
    p = p + geom_point(cex=3)
} else {
    p = p +
        geom_point(data=subset(mdspca, !grepl("Bugala", mdspca$island), cex=3)) +
        geom_point(data=subset(mdspca,  grepl("Bugala", mdspca$island), cex=3))
}

p = p +
    xlab(paste0(component.name, "1", prop.explained[1])) +
    ylab(paste0(component.name, "2", prop.explained[2])) +
    coord_fixed() +
    scale_shape_manual(values=seq(0, length(table(mdspca$island)))) +
    scale_color_manual(name=NULL,
                         breaks=levels(mdspca$island),
                         values=legend.colors) +
    guides(shape = guide_legend(override.aes = list(color=legend.colors), title="Site"),
           color=FALSE) +
    theme_bw(base_size=15) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=2))

ggsave(p, file=out.file, height=7, width=7)

# --- Write file of BUGALA individuals split between mainland- and island-like groups

if (split.bugala) {

    bugala.ml = mdspca[grepl("Bugala", mdspca$island) & mdspca$C2 < 0,]
    bugala.is = mdspca[grepl("Bugala", mdspca$island) & mdspca$C2 >= 0,]

    write.table(bugala.ml$mosquito.id, file="data/ssese.samples.is-BUGALAML.txt",
        quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
    write.table(bugala.is$mosquito.id, file="data/ssese.samples.is-BUGALAIS.txt",
        quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")

    write.table(bugala.ml$IID, file="data/ssese.seqids.is-BUGALAML.txt",
        quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
    write.table(bugala.is$IID, file="data/ssese.seqids.is-BUGALAIS.txt",
        quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")

    ind.info.bugala = read.csv("data/ssese_individual_info.csv")

    ind.info.bugala[ind.info.bugala$mosquito_id %in% bugala.ml$mosquito.id,]$island =
        "BUGALAML"
    ind.info.bugala[ind.info.bugala$mosquito_id %in% bugala.is$mosquito.id,]$island =
        "BUGALAIS"

    write.csv(ind.info.bugala, file="data/ssese_individual_info_bugala_split.csv",
        quote=FALSE, row.names=FALSE)

    ind.info.bugala.simple = read.table("data/ssese_individual_info_simple.txt")

    ind.info.bugala.simple[ind.info.bugala.simple$V1 %in% bugala.ml$mosquito.id,]$V3 =
        "BUGALAML"
    ind.info.bugala.simple[ind.info.bugala.simple$V1 %in% bugala.is$mosquito.id,]$V3 =
        "BUGALAIS"

    write.table(ind.info.bugala.simple,
        file="data/ssese_individual_info_simple_bugala_split.txt",
        quote=FALSE, col.names=FALSE, row.names=FALSE, sep=" ")
}

if (split.bugala == FALSE) {

    # --- Graph with inset for full genome plots

    if (grepl("noinv", mdspca.in)) {

        inset.x.range = range(mdspca[!mdspca$is.island,]$C1)
        inset.y.range = range(mdspca[!mdspca$is.island,]$C2)

        if (is.mds) {
            inset.x.range = inset.x.range + c(-0.001, 0.001)
            inset.y.range = inset.y.range + c(-0.001, 0.001)
        } else {
            inset.x.range = inset.x.range + c(-0.025, 0.025)
            inset.y.range = inset.y.range + c(-0.025, 0.025)
        }

        p.noleg = ggplot(mdspca, aes(C1, C2, col=is.island, pch=island)) +
            geom_rect(xmin=inset.x.range[1], xmax=inset.x.range[2],
                      ymin=inset.y.range[1], ymax=inset.y.range[2],
                      fill=NA, col="lightgrey", lty=2) +
            geom_point(cex=3) +
            xlab(paste0(component.name, "1", prop.explained[1])) +
            ylab(paste0(component.name, "2", prop.explained[2])) +
            coord_fixed() +
            scale_shape_manual(values=seq(0, length(table(mdspca$island)))) +
            scale_color_discrete(name=NULL,
                                 breaks=c(TRUE, FALSE),
                                 labels=c("Island", "Mainland")) +
            guides(shape = FALSE, color=FALSE) +
            theme_bw()

        p.inset = ggplot(mdspca, aes(C1, C2, col=is.island, pch=factor(island))) +
            geom_point(cex=3) +
            xlab(paste0(component.name, "1", prop.explained[1])) +
            ylab(paste0(component.name, "2", prop.explained[2])) +
            coord_fixed() +
            scale_shape_manual(values=seq(0, length(table(mdspca$island)))) +
            scale_color_discrete(name=NULL,
                                 breaks=c(TRUE, FALSE),
                                 labels=c("Island", "Mainland")) +
            guides(shape = FALSE, color=FALSE) +
            theme_bw() +
            scale_x_continuous(limits=inset.x.range, expand = c(0, 0)) +
            scale_y_continuous(limits=inset.y.range, expand = c(0, 0)) +
            theme(axis.title = element_blank(),
                  axis.text  = element_blank(),
                  axis.ticks = element_blank(),
                  panel.border = element_rect(color="lightgrey", fill=NA,
                    linetype=3, size=2))

        # Function to extract legend
        # From: http://stackoverflow.com/questions/11883844/
        # inserting-a-table-under-the-legend-in-a-ggplot2-histogram
        get.legend = function(mygplot) {
            mygplot = mygplot + theme(legend.box.margin=margin(1,1,5,1, unit="cm"))
            tmp = ggplot_gtable(ggplot_build(mygplot))
            leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
            legend = tmp$grobs[[leg]]
            return(legend)
        }

        legend = get.legend(p)

        lay = rbind(c(1,1,1,2,2),
                    c(1,1,1,2,2),
                    c(1,1,1,3,4),
                    c(1,1,1,3,4))

        p.full = grid.arrange(p.noleg, p.inset, legend, grob(), layout_matrix = lay)
        ggsave(filename=gsub("pdf", "inset.pdf", out.file), p.full)

    }

    # Export mds now that info on individuals has been added.

    if (is.mds) {
        out.txt = gsub("results/(.*).mds$", "reports/\\1.mds.txt", mdspca.in)
        out.txt = gsub("data/(.*).mds$", "reports/\\1.mds.txt", out.txt)

        write.table(mdspca, file=out.txt)
    }
}
