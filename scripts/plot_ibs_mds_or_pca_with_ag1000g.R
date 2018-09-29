#!/usr/bin/Rscript

options(stringsAsFactors = FALSE)

library(ggplot2)
library(RColorBrewer)
library(colorspace)
library(scales)
library(grid)
library(gridExtra)
library(lattice)
library(plyr)

args = commandArgs(trailingOnly = TRUE)
mdspca.in = args[1] # e.g. "data/ssese_with_ag1000g/all.pass.snp.phased.ag1000g.noinv.mds"
   # or "data/ssese_with_ag1000g/chr3.pass.snp.phased.ag1000g.noinv.eigenvec"
   # or "data/ssese_with_ag1000g/chr3.pass.snp.phased.ag1000g.strict.noinv.withM.eigenvec"
ind.info = read.table("data/ssese_individual_info_simple.txt")

# Bring in plate info to infer sequencing ID (e.g. ssese42.PE)
plates = read.table("data/sample_to_seq_id_mapping.txt")
names(plates) = c("mosquito.id", "seq.id")

# Reduce individual info to just sequenced animals
ind.info = ind.info[ind.info$V1 %in% plates$mosquito.id,]

ind.info.ag1000g = read.table("data/ag1000g.phase1.ar3/samples.all.txt",
    sep="\t", quote="", header=TRUE)

# Figure out if we're plotting MDS or PCA
is.mds = FALSE
if (grepl("mds", mdspca.in)) {
    is.mds = TRUE
    num.components = 10
}

if (is.mds) {
    mdspca = read.table(mdspca.in, header = T)
} else {
    mdspca = read.table(mdspca.in)
    # Name as C to be parallel with MDS column names
    num.components = ncol(mdspca) - 2
    names(mdspca) = c("IID", "FID", paste0("C", 1:num.components))
}

mdspca$pops = do.call(c, lapply(mdspca$IID, function(x) {
    pop.ss = as.character(ind.info[ind.info$V2 == x, 3])
    pop.ag = paste(as.character(ind.info.ag1000g[ind.info.ag1000g$ox_code == x,]$country),
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

mdspca$pops[mdspca$pops %in% c("BANDA", "BUKASA", "BUGALA", "NSADZI", "SSERINYA")] =
    "LVB - Ssese Island"
mdspca$pops[mdspca$pops %in% c("BUWAMA", "KAZZI", "KIYINDI", "MITYANA")] =
    "LVB - Mainland"

# Change Uganda to "Uganda - Ag1000G"
mdspca$pops[mdspca$pops == "Uganda_S"] = "Uganda - Ag1000G [gambiae]"

# Add gambiae or coluzzi
mdspca$pops = gsub("_S$",   " [gambiae]",  mdspca$pops)
mdspca$pops = gsub("_M/S$", " [hybrid]",   mdspca$pops)
mdspca$pops = gsub("_M$",   " [coluzzii]", mdspca$pops)

# Lump Guinea hybrid (N=1) into Guinea
mdspca$pops[mdspca$pops == "Guinea [hybrid]"] = "Guinea [gambiae]"

# Lump all Guinea-Bissau into just Guinea-Bissau
mdspca$pops[grepl("Guinea-Bissau", mdspca$pops)] = "Guinea-Bissau"

# Remove species info for Kenya a la Ag1000G paper
mdspca$pops[mdspca$pops == "Kenya [gambiae]"] = "Kenya"

ordered.pops = names(table(mdspca$pops))

mdspca$pops = factor(mdspca$pops,
                  levels = c(grep("LVB", ordered.pops, invert=TRUE, value=TRUE),
                             grep("LVB", ordered.pops, value=TRUE)))

# Re-order to match now that rearranging done
ordered.pops = levels(mdspca$pops)

mdspca$is.ssese = grepl("LVB", mdspca$pops)

pal = brewer.pal(n = length(table(mdspca$pops)) - 2, "Set1")
new.pal = rep(NA, length(ordered.pops))

# Fill in non-LVB colors
new.pal[grep ("LVB", ordered.pops, invert=TRUE)] = pal
# Set Ssese islands to blue and mainland to red
new.pal[grep ("Ssese Island", ordered.pops)] = "#00BFC4"
new.pal[grep ("Mainland", ordered.pops)]     = "#F8766D"

# Set Kenya to grey and Uganda to mainland LVB color
new.pal[grep ("Kenya",  ordered.pops)] = "#999999"
new.pal[grep ("Uganda", ordered.pops)] = "#F8766D"

if (is.mds) {
    component.name = "Component "
} else {
    component.name = "PC"
}

# Add proportion of variance explained, if PCA
prop.explained = rep("", num.components)
if (grepl("eigenvec$", mdspca.in)) {
    val.in = gsub("vec$", "val", mdspca.in)
    val = read.table(val.in)
    prop.explained = paste0(" (", percent(val$V1 / sum(val$V1)), ")")
}

# Figure out output prefix
if (is.mds) {
    # For MDS, since only 1 vs. 2 is plotted, output prefix is full filename
    out.prefix = gsub("results/(.*).mds$", "reports/\\1.ibs_mds_plot.pdf", mdspca.in)
    out.prefix = gsub("data/(.*).mds$", "reports/\\1.ibs_mds_plot.pdf", out.prefix)
} else {
    out.prefix = gsub("data/(.*).eigenvec$", "reports/\\1.pca_plot.", mdspca.in)
}

out.prefix = gsub("/ssese_with_ag1000g.*/", "/", out.prefix)

if (grepl("missingtoref", mdspca.in)) {
    out.prefix = gsub("(chr[^\\.]*\\.)", "\\1missingtoref.", out.prefix)
}

darken = function(color, factor=1.4){
    col = col2rgb(color)
    col = col/factor
    col = rgb(t(col), maxColorValue=255)
    col
}

# Italicize species names
#ordered.pops = gsub("gambiae", expression(italic("gambiae")), ordered.pops)
#ordered.pops = gsub("coluzzi", expression(italic("coluzzi")), ordered.pops)

make.italic <- function(x) {
    as.expression(lapply(x, function(y) {
        if (grepl("gambiae", y) | grepl("coluzzii", y)) {
            part1 = gsub("(.*\\[).*\\]", "\\1", y)
            part2 = gsub(".*\\[(.*)\\]", "\\1", y)
            bquote(.(part1)*italic(.(part2))*"]")
        } else {
            y
        }
    }
    ))
}

# Actually, combine all Ugandan mainland samples for better readability and tweak names
mdspca$pops = revalue(mdspca$pops,
                        c("Uganda - Ag1000G [gambiae]" = "Mainland Uganda [gambiae]",
                          "LVB - Mainland"             = "Mainland Uganda [gambiae]"))
mdspca$pops = revalue(mdspca$pops,
                        c("LVB - Ssese Island" = "Ssese Islands [gambiae]"))

ordered.pops = c(ordered.pops[1:7], "Mainland Uganda [gambiae]", ordered.pops[10])
ordered.pops[ordered.pops == "LVB - Ssese Island"] = "Ssese Islands [gambiae]"
new.pal = new.pal[-c(9)]

p = ggplot(mdspca, aes(C1, C2, col=pops, fill=pops, shape=is.ssese)) +
        #=geom_point(col='black', cex=2.25) +
        geom_point(cex=1.5) +
        #geom_point(aes(alpha=is.ssese), col='black', cex=0.25, pch=1) +
        xlab(paste0(component.name, "1", prop.explained[1])) +
        ylab(paste0(component.name, "2", prop.explained[2])) +
        scale_color_manual(values = darken(new.pal), labels=make.italic(ordered.pops)) +
        scale_fill_manual(values = new.pal) +
        #scale_alpha_manual(values = c(0, 0.75),
        #                   breaks=c(FALSE, TRUE),
        #                   labels=c("", "LVB Sample")) +
        scale_shape_manual(values = c(21, 21)) +
        coord_fixed() +
        theme_bw() +
        theme(legend.text.align = 0) +
        guides(shape=FALSE, fill=FALSE,
               color = guide_legend(override.aes = list(
                                        col=darken(new.pal),
                                        fill=new.pal,
                                        size=3,
                                        shape=c(21,21)[1 + grepl("LVB", ordered.pops)]),
                                    title="Populations") #,
               #alpha = guide_legend(override.aes = list(
               #                         size=0.5,
               #                         shape=21),
               #                     title="")
        )

# Only save for MDS, which is just component 1 vs 2
if (is.mds) {
    ggsave(file=out.prefix)
}

# --- Graph with inset for full genome plots for both MDS and PCA (dimension 1 vs 2 only)

if (grepl("noinv", mdspca.in)) {

    inset.x.range = range(mdspca[mdspca$is.ssese,]$C1)
    inset.y.range = range(mdspca[mdspca$is.ssese,]$C2)

    if (is.mds) {
        inset.x.range = inset.x.range + c(-0.001, 0.001)
        inset.y.range = inset.y.range + c(-0.001, 0.001)
    } else {
        inset.x.range = inset.x.range + c(-0.003, 0.003)
        inset.y.range = inset.y.range + c(-0.003, 0.003)
    }

    p.noleg = p + geom_rect(xmin=inset.x.range[1], xmax=inset.x.range[2],
                        ymin=inset.y.range[1], ymax=inset.y.range[2],
                        fill=NA, col="lightgrey", lty=2) +
                    guides(shape = FALSE, color=FALSE, alpha=FALSE)

    p.inset = p + guides(shape = FALSE, color=FALSE, alpha=FALSE) +
                    scale_x_continuous(limits=inset.x.range, expand = c(0, 0)) +
                    scale_y_continuous(limits=inset.y.range, expand = c(0, 0))

    # Function to extract legend
    # From: http://stackoverflow.com/questions/11883844/inserting-a-table-under-the-legend-in-a-ggplot2-histogram
    get.legend = function(mygplot) {
        mygplot = mygplot + theme(legend.box.margin=margin(0,1,0,0, unit="cm"))
        tmp = ggplot_gtable(ggplot_build(mygplot))
        leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend = tmp$grobs[[leg]]
        return(legend)
    }

    legend = get.legend(p)

    lay = rbind(c(1,1,1,1,2,2),
                c(1,1,1,1,2,2),
                c(4,4,4,4,3,3))

    p.full = grid.arrange(p.noleg, p.inset, legend, grob(), layout_matrix = lay)
    ggsave(filename=paste0(gsub("pdf$", "", out.prefix), "inset.pdf"), p.full,
        width=8, height=5)
}

# --- For PCA, plot additional dimensions (Only 1-6)

if (!is.mds) {

    for (pc.a in 1:5) {
         for (pc.b in (pc.a + 1):6) {

            p = ggplot(mdspca, aes_string(paste0("C", pc.a), paste0("C", pc.b),
                        col="pops", fill="pops", shape="is.ssese")) +
                    geom_point(col='black', cex=2.25) +
                    geom_point(cex=1.5) +
                    #geom_point(aes(alpha=is.ssese), col='black', cex=0.25, pch=1) +
                    xlab(paste0(component.name, "1", prop.explained[1])) +
                    ylab(paste0(component.name, "2", prop.explained[2])) +
                    scale_color_manual(values = darken(new.pal),
                                       labels=make.italic(ordered.pops)) +
                    scale_fill_manual( values = new.pal) +
                    #scale_alpha_manual(values = c(0, 0.75),
                    #                   breaks=c(FALSE, TRUE),
                    #                   labels=c("", "LVB Sample")) +
                    scale_shape_manual(values = c(21, 21)) +
                    coord_fixed() +
                    theme_bw() +
                    theme(legend.text.align = 0) +
                    guides(shape=FALSE, fill=FALSE,
                           color = guide_legend(override.aes = list(
                                                    col=darken(new.pal),
                                                    fill=new.pal,
                                                    size=3,
                                                    shape=c(21,21)[1 + grepl("LVB", ordered.pops)]),
                                                title="Populations") #,
                           #alpha = guide_legend(override.aes = list(
                           #                         size=0.5,
                           #                         shape=21),
                           #                     title="")
                    )

            ggsave(file=paste0(out.prefix, pc.a, "vs", pc.b, ".pca_plot.pdf"))
        }
    }
}

# --- Make little version for text

pc.a = 1
pc.b = 2

p = ggplot(mdspca, aes_string(paste0("C", pc.a), paste0("C", pc.b),
            col="pops", fill="pops", shape="is.ssese")) +
        geom_point(col='black', cex=2.25) +
        geom_point(cex=1.5) +
        #geom_point(aes(alpha=is.ssese), col='black', cex=0.25, pch=1) +
        xlab(paste0(component.name, "1", prop.explained[1])) +
        ylab(paste0(component.name, "2", prop.explained[2])) +
        scale_color_manual(values = darken(new.pal),
                           labels=make.italic(ordered.pops)) +
        scale_fill_manual( values = new.pal) +
        #scale_alpha_manual(values = c(0, 0.75),
        #                   breaks=c(FALSE, TRUE),
        #                   labels=c("", "LVB Sample")) +
        scale_shape_manual(values = c(21, 21)) +
        coord_fixed() +
        theme_bw(base_size=15) +
        theme(legend.text.align = 0,
              panel.border = element_rect(colour = "black", fill=NA, size=2)) +
        guides(shape=FALSE, fill=FALSE,
               color = guide_legend(override.aes = list(
                                        col=darken(new.pal),
                                        fill=new.pal,
                                        size=3,
                                        shape=c(21,21)[1 + grepl("LVB", ordered.pops)]),
                                    title="Populations") #,
               #alpha = guide_legend(override.aes = list(
               #                         size=0.5,
               #                         shape=21),
               #                     title="")
        )

ggsave(file=paste0(out.prefix, pc.a, "vs", pc.b, ".pca_plot.fortext.pdf"),
    height=7, width=7)

# --- Make scree plot

if (grepl("eigenvec$", mdspca.in)) {
    val.in = gsub("vec$", "val", mdspca.in)
    val = read.table(val.in)
    val$prop.var.exp = val$V1 / sum(val$V1)

    p = ggplot(val, aes(seq_along(prop.var.exp), prop.var.exp)) +
        geom_bar(stat="identity", col='black') +
        ylab("Proportion of Variation Explained") +
        xlab("Principal Component") +
        theme_bw() +
        scale_x_continuous(breaks = 1:num.components) +
        scale_fill_manual(values=c("white", "black"),
                          breaks=c(FALSE, TRUE),
                          guide=FALSE)
    ggsave(filename=paste0(gsub("pdf$", "", out.prefix), "scree.pdf"), p,
        width=4.5, height=6)
}
