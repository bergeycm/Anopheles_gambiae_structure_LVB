#!/usr/bin/env Rscript

# module load R/3.3.0

# ========================================================================================
# --- Compute NJ, etc. around regions of interest
# ========================================================================================

options(stringsAsFactors=FALSE)

args = commandArgs(trailingOnly=TRUE)

chr          = args[1]
center       = as.numeric(args[2])
window       = as.numeric(args[3])
name         = args[4]
do.diploid   = as.logical(args[5])
with.ag1000g = as.logical(args[6])
do.snpeff    = as.logical(args[7])

# New sweeps:
# chr="2L"; center=34044820; window=10000; name="2L_34Mb";    do.diploid=FALSE; with.ag1000g=TRUE; do.snpeff=FALSE
# chr="X";  center=9222217;  window=10000; name="chrX.X_9Mb"; do.diploid=FALSE; with.ag1000g=TRUE; do.snpeff=FALSE

# Small window for testing:
# chr="X";  center=9222217;  window=1000;  name="chrX.X_9Mb"; do.diploid=TRUE; with.ag1000g=TRUE; do.snpeff=FALSE

# Known insecticide genes
# chr="2R"; center=28501972; window=10000; name="2R_CYP6P2";  do.diploid=FALSE; with.ag1000g=TRUE; do.snpeff=FALSE
# chr="X";  center=15241718; window=10000; name="X_CYP9K1";   do.diploid=FALSE; with.ag1000g=TRUE; do.snpeff=FALSE
# chr="3R"; center=28598038; window=10000; name="3R_GSTE";    do.diploid=FALSE; with.ag1000g=TRUE; do.snpeff=FALSE

# Only include HIGH SNPs in snpEff data?
just.high = FALSE
if (window > 10000) {
    just.high = TRUE
}

include.cluster.bar = FALSE

library(adegenet)
library(ape)
library(RColorBrewer)
library(dendextend)

# ----------------------------------------------------------------------------------------
# --- Read in data, extract region of interest, annotate with snpEff
# ----------------------------------------------------------------------------------------

# --- Extract region of interest

chrs = c("2L", "2R", "3L", "3R", "X")
chr.num = which (chr == chrs)
if (chr.num == 5) {
    chr.num = "X"
}

start = center - window / 2
end   = center + window / 2

# Suffix that covers if diploid or haploid and if Ag1000G is included
if (do.diploid) {
    class.suffix = "diploid"
} else {
    class.suffix = "haploid"
}

if (with.ag1000g) {
    class.suffix = paste(class.suffix, "with-ag1000g", sep=".")
} else {
    class.suffix = paste(class.suffix, "ssese-only", sep=".")
}

if (do.snpeff) {
    class.suffix = paste(class.suffix, "with-snpeff", sep=".")
} else {
    class.suffix = paste(class.suffix, "sans-snpeff", sep=".")
}

# Prefix for temporary reduced VCF
in.prefix = paste0("data/for_dendrograms/chr", chr, ".", start, "-", end,
    ".", class.suffix)

# Output prefix
dir.create("data/for_dendrograms", showWarnings=FALSE)
out.prefix = paste0("results/for_dendrograms/chr", chr, ".", name, ".",
    window, ".", class.suffix, ".haps")

# Reduce VCF
if (with.ag1000g) {
    #in.vcf = paste0("data/ssese_with_ag1000g/chr", chr, ".pass.snp.phased.ag1000g.vcf")
    in.vcf = paste0("data/ssese_with_ag1000g/ssese_with_ag1000g.", chr, ".flt.strict.vcf")
} else {
    # Note: Unphased VCF. Only used for snpEff. Phase info read from haps files directly
    in.vcf = paste0("data/chr", chr, ".pass.snp.vcf")
}

vcftools.cmd = paste0("vcftools ",
    "--vcf ", in.vcf, " ",
    "--chr ", chr, " --from-bp ", start, " --to-bp ", end, " ",
    "--max-missing 1 ",
    "--remove data/ag1000g.phase1.ar3/excluded.individuals.notM.txt ",
    "--remove-indels --recode --recode-INFO-all ",
    "--out ", in.prefix)

system(vcftools.cmd)

# Check to see if any output remains
if (file.exists(paste0(in.prefix, ".recode.vcf"))) {
    oldw = getOption("warn")
    options(warn = -1)
    grep.result = as.numeric(system(paste0("grep -v '^#' ", in.prefix, ".recode.vcf | wc -l"),
        intern=TRUE))
    options(warn = oldw)
} else {
    grep.result = 0
}

if (as.numeric(grep.result) == 0) {
    cat(paste("ERROR: No data in this region given window size. Aborting.\n"))
    # Fake output files
    system(paste0("touch ", in.prefix, ".clusters.txt"))
    system(paste0("touch ", out.prefix, ".dendrogram.pdf"))
    system(paste0("touch ", out.prefix, ".fancydendro.pdf"))
    quit()
}

# # Replace chromosome number with actual chromosome in VCF
# sed.cmd = paste0('sed -e "s/^', chr.num, "/", chr, '/" ', in.prefix, ".recode.vcf ",
#     '| sed -e "s/^\\(\\S*\\)\t\\(\\S*\\)\t\\./\\1\t\\2\t\\1:\\2/"',
#     "> ", in.prefix, ".recode.fix.vcf")
#
# system(sed.cmd)

# Run snpEff
if (do.snpeff) {
    snpeff.cmd = paste0("java -Xmx4g -jar ~/bin/snpEff/snpEff.jar ",
        "-v -v -v AgamP4.29 -dataDir ~/snpEff_DB ",
        in.prefix, ".recode.vcf ",
        "> ", in.prefix, ".snp.eff.vcf")

    system(snpeff.cmd)
}

# Replace actual chromosome with chromosome number in VCF to be input into PLINK
# Also add SNP name to third column in form chr:pos
sed.cmd = paste0('sed -e "s/^', chr, "/", chr.num, '/" ', in.prefix, ".recode.vcf ",
    '| sed -e "s/^\\(\\S*\\)\t\\(\\S*\\)\t\\./\\1\t\\2\t\\1:\\2/"',
    "> ", in.prefix, ".recode.fix.vcf")

system(sed.cmd)

# Convert from VCF to PLINK RAW format so that it can be read by read.PLINK()
# as a genlight object
# and remove undesireable individuals (like lab crosses and M form mosquitoes)

# First command is just to make a MAP file by going through PED intermediary
plink.cmd.1 = paste0("plink --vcf ", in.prefix, ".recode.fix.vcf ",
    "--remove-fam data/ag1000g.phase1.ar3/excluded.individuals.notM.txt ",
    "--recode ",
    "--out ", in.prefix)

system(plink.cmd.1)

plink.cmd.2 = paste0("plink --file ", in.prefix, " ",
    "--recode A ",
    "--out ", in.prefix)

system(plink.cmd.2)

# --- Read in file

gen = read.PLINK(paste0(in.prefix, ".raw"),
    map.file = paste0(in.prefix, ".map"),
    parallel = TRUE, n.cores = 8)

# Remove ssese99.PE from gen object
if ("ssese99.PE" %in% gen$ind.names) {
    gen = gen[-which(gen$ind.names == "ssese99.PE"),]
}

# Fix SNP names
map.file = read.table(paste0(in.prefix, ".map"))
snp.names = paste(chr, map.file$V4, sep=":")

# gen$loc.names = gsub("\\.", "", gen$loc.names)
# gen$loc.names = paste0(snp.names, gen$loc.names)

# ----------------------------------------------------------------------------------------

# --- Read in functional data from snpEff output

if (do.snpeff) {

    # This skips multi-allelic SNPs and boring intronic stuff

    #awk.cmd = paste0("awk 'BEGIN {OFS = \"\t\"} { if ($5 !~ \",\" && $2 > ", start, " && $2 < ", end, ") ",
    #                 "print $1,$2,$4,$5,$8 }' ", in.prefix, ".snp.eff.vcf | ",
    #                 "grep -v \"intron_variant|MODIFIER\" | grep -v \"WARNING\" | tr \"\\|\" \"\t\"")

    awk.cmd = paste0("awk 'BEGIN {OFS = \"\t\"} { if ($5 !~ \",\" && $2 > ", start, " && $2 < ", end, ") ",
                     "print $1,$2,$4,$5,$8 }' ", in.prefix, ".snp.eff.vcf | ",
                     "grep -v \"intron_variant|MODIFIER\" | grep -v \"WARNING\"")

    if (as.numeric(system(paste(awk.cmd, "| wc -l"), intern=TRUE)) != 0) {
        snpeff = read.table(pipe(awk.cmd), header=FALSE, fill=TRUE)
    } else {
        # No suitable SNPs with functional info
        snpeff = data.frame(chr=character(), pos=integer(), ref=character(), alt=character())
    }

    # The anno column is split into multiple annotations (that could exist) by
    # splitting on 'ANN', but only the first is used. Then the annotations are
    # be further split into the annotation parts by splitting on "|".

    names(snpeff) = c("chr", "pos", "ref", "alt", "anno")

    get.snpeff.bar.data = function (snpeff.locus.num) {

        if (snpeff.locus.num %% 10 == 0) {
            write(paste("Processing line", snpeff.locus.num, "of", nrow(snpeff)), stderr())
        }

        snpeff.loc = snpeff[snpeff.locus.num,]

        # Reduce snpEff data further.
        anno.split = strsplit(gsub(".*ANN=(.*)", "\\1", snpeff.loc$anno), split=',')[[1]]
        # To just top annotation
        anno.split = anno.split[1]
        anno.df = strsplit(anno.split, "\\|")[[1]]

        loc = snpeff.loc$pos
        alt = snpeff.loc$alt

        gen.idx = which(gen$loc.names == paste0(chr.num, ":", loc, "_", alt))

        snps = do.call(c, lapply(gen$gen, function (x) { as.integer(x[gen.idx]) }))

        if (length(table(snps)) == 1) {
            return (NA)
        }

        snp.df = data.frame(ind=gen$ind.names, genotype=snps)

        # Remove upstream gene variants
        if (anno.df[2] == "upstream_gene_variant") {
            return (NA)
        }

        if (just.high) {
            if (anno.df[3] != "MODERATE" & anno.df[3] != "HIGH") {
                return(NA)
            }
        }

        snpeff.info.string = paste0(anno.df[4], " ", gsub("_", " ", anno.df[2]), ": ",
            gsub("p.", "", anno.df[11]), " (", gsub("c.", "", anno.df[10]), ")")

        result = list(snp.df)
        names(result) = snpeff.info.string

        return (result)
    }

    snpeff.data = lapply(1:nrow(snpeff), get.snpeff.bar.data)

    snpeff.data = snpeff.data[
        do.call(c, lapply(1:length(snpeff.data),
            function (x) { !is.na(snpeff.data[[x]][1]) })
        )
    ]

    # length(snpeff.data)
}

# ----------------------------------------------------------------------------------------
# --- Make dendrograms
# ----------------------------------------------------------------------------------------

dir.create("results/for_dendrograms", showWarnings=FALSE)

if (do.diploid) {

    # --- Compute distance matrix

    x = seploc(gen, n.block=10, parallel=FALSE)

    lD = lapply(x, function(e) dist(as.matrix(e)))

    D = Reduce("+", lD)

    # Plot NJ tree

    pdf(paste0(out.prefix, ".njtree.pdf"), height=20, width=20)
        plot(nj(D), type="unrooted")
    dev.off()

} else {

    # --- Make reduced hap file and convert to FASTA format

    if (with.ag1000g) {
        orig.haps = paste0("data/ssese_with_ag1000g/chr", chr,
            ".pass.snp.phased.ag1000g.haps")
        in.haps = gsub("data/ssese_with_ag1000g", "data/for_dendrograms",
            gsub("haps$", paste0(chr, ".", start, "-", end, ".haps"), orig.haps))
    } else {
        orig.haps = paste0("data/chr", chr, ".pass.snp.phased.haps")
        in.haps = gsub("data", "data/for_dendrograms",
            gsub("haps$", paste0(chr, ".", start, "-", end, ".haps"), orig.haps))
    }

    awk.cmd = paste0("awk -v astart=", start, " -v aend=", end,
        " '{ if ($3 > astart && $3 < aend) print $0 }' ",
        orig.haps, " > ", in.haps)
    system(awk.cmd)

    cp.cmd = paste0("cp ", gsub("haps$", "sample", orig.haps), " ",
        gsub(".haps$", ".sample", in.haps))
    system(cp.cmd)

    hap.to.fa.cmd = paste0("python scripts/hap_to_fa.py ", in.haps)
    paste0("Running command [", hap.to.fa.cmd, "]")
    system("python --version")
    system(paste("ls -lhrt", in.haps))
    system(hap.to.fa.cmd)
    paste0("Finished conversion to FASTA.")

    # --- Read in data and compute distance matrix

    fa.in = gsub(".haps$", ".fasta", in.haps)

    fa = read.FASTA(fa.in)

    # Remove excluded individuals (SY, etc.)
    bad.inds = read.table("data/ag1000g.phase1.ar3/excluded.individuals.notM.txt")$V1
    bad.inds = c(paste0(bad.inds, "_A"), paste0(bad.inds, "_B"))
    bad.indices = match(bad.inds, names(fa))
    bad.indices = bad.indices[!is.na(bad.indices)]
    if (length(bad.indices) > 0) {
        fa = fa[-c(bad.indices)]
    }

    D = dist.dna(fa, model='raw') * length(as.character(fa[[1]]))

    pdf(paste0(out.prefix, ".njtree.pdf"), height=20, width=20)
        plot(nj(D), type="unrooted")
    dev.off()

}

# --- Perform hierarchical cluster analysis

hc = hclust(D)

# --- Figure out population

ind.info = read.table("data/ssese_individual_info_simple_bugala_split.txt")

# Ag1000G info is read even if it is never used
ind.info.ag1000g = read.table("data/ag1000g.phase1.ar3/samples.all.txt",
    sep="\t", quote="", header=TRUE)

pops = do.call(c, lapply(gsub("_.*", "", hc$labels), function(x) {
    pop.ss = as.character(ind.info[ind.info$V2 == x, 3])
    pop.ag = paste(as.character(ind.info.ag1000g[ind.info.ag1000g$ox_code == x,]$country),
                   as.character(ind.info.ag1000g[ind.info.ag1000g$ox_code == x,]$m_s),
                    sep="_")
    if (length(pop.ss)) {
        return(pop.ss)
    } else if (length(pop.ag)) {
        return(pop.ag)
    } else if (length(grepl("SY", x)) > 0) {
        return("SY")
    } else {
        return(x)
    }
}))

# Change Uganda to "Uganda - Ag1000G"
pops[pops == "Uganda_S"] = "Uganda - Ag1000G [gambiae]"

# Add gambiae or coluzzi
pops = gsub("_S$",   " [gambiae]",  pops)
pops = gsub("_M/S$", " [hybrid]",   pops)
pops = gsub("_M$",   " [coluzzii]", pops)

# Lump Guinea hybrid (N=1) into Guinea
pops[pops == "Guinea [hybrid]"] = "Guinea [gambiae]"

# Lump all Guinea-Bissau into just Guinea-Bissau
pops[grepl("Guinea-Bissau", pops)] = "Guinea-Bissau"

# Remove species info for Kenya a la Ag1000G paper
pops[pops == "Kenya [gambiae]"] = "Kenya"

# --- Figure out geographic class of site

geo.class = pops
geo.class[pops %in% c("BUGALAML", "BUWAMA", "KAZZI", "KIYINDI", "MITYANA")] = "mainland"
geo.class[pops %in% c("BANDA", "BUGALAIS", "BUKASA", "NSADZI", "SSERINYA")] = "island"

# Set Ag1000G Ugandan samples to "mainland"
geo.class[geo.class == "Uganda - Ag1000G [gambiae]"] = "mainland"

# --- Make pop list that lumps all Ugandan animals into one group

pops.simple = pops
pops.simple[pops == toupper(pops) | grepl("Uganda", pops)] = "Uganda"

# --- Color by geographic class or population or country

if (with.ag1000g) {
    categories = factor(pops.simple)
} else {
    categories = factor(pops)
}
# Alternative:
# categories = factor(geo.class)

num.cats = length(levels(categories))

pal = brewer.pal(num.cats, "Set1")
getPal = colorRampPalette(pal)

if (with.ag1000g) {
    cols = getPal(num.cats)[factor(pops.simple)]
} else {
    cols = getPal(num.cats)[factor(pops)]
}
# Alternative:
# cols = getPal(num.cats)[factor(geo.class)]

get.color = function (label) {
    cols[which(hc$labels == label)]
}

color.labels = function(n) {
    if (is.leaf(n)) {
        # Get label
        label = attr(n, "label")

        # Figure out label color
        label.color = get.color(label)
        attr(n, "nodePar") = list(lab.col = label.color)
        attr(n, "label") = pops[which(hc$labels == label)]
    }
    n
}

hc.d.col = dendrapply(as.dendrogram(hc), color.labels)

pdf(paste0(out.prefix, ".dendrogram.pdf"), height=20, width=60)
    plot(hc.d.col,
        main = paste0(name, " - chr", chr, ":", start, "-", end),
        type = "rectangle")
dev.off()

# --- Checkpoint

save.image(file=paste0(in.prefix, ".Rdata"))

# ----------------------------------------------------------------------------------------
# --- Plot fancy dendrogram
# ----------------------------------------------------------------------------------------

# devtools::install_github('talgalili/dendextend')

# Scale height so tree looks OK with bars
# Get max height, see http://stackoverflow.com/a/25665440
max.height = max(cophenetic(hc))
hc$height = (hc$height / max.height) * 10

hc.d = as.dendrogram(hc)

hc.d.inds = labels(hc.d)

# Cut tree at height equivalent of:
# 0.0004 SNP differences per bp = 0.4 SNP differences per kb
cut.height = (window / 1000) * 0.4
# Also scale it to new scale (out of 10)
cut.height = (cut.height / max.height) * 10
hc.d.cut = cutree(hc.d, h=cut.height, order_clusters_as_data=FALSE)

# Now just use same haplotypes from cutting lower
hc.d.cut.k10 = hc.d.cut
### # Also cut into 10 clusters for map of haplotypes (map o' haps)
### hc.d.cut.k10 = cutree(hc.d, k=10, order_clusters_as_data=FALSE)

# Grab clusters containing certain number of individuals
if (with.ag1000g) {
    clust.ind.cutoff = 4
} else {
    clust.ind.cutoff = 2
}
clusts = table(hc.d.cut)[table(hc.d.cut) >= clust.ind.cutoff]
# Set boring clusters-of-one to NA
hc.d.cut[!hc.d.cut %in% names(clusts)] = NA

# Set labels to population (includes countries for Ag1000G and islands/sites for Ssese)
pops.hc.d = do.call(rbind, lapply(hc.d.inds, function(x) {
    pops[which(hc$labels == x)] }))
labels(hc.d) = pops.hc.d

# Make data.frame of individuals with haplotypes for map
haps.k10 = data.frame(ind = names(hc.d.cut.k10),
                      hap = gsub("_.*", "", names(hc.d.cut.k10)),
                      hap.index = hc.d.cut.k10,
                      pop = pops.hc.d)
write.table(haps.k10, file=paste0(out.prefix, ".haplotypes.txt"))

# Turn off labels (throws warning)
labels(hc.d) = NA

# Color labels
cols.hc.d = do.call(rbind, lapply(hc.d.inds, get.color))

if (with.ag1000g) {
    # Figure out Uganda's color
    uganda.col = getPal(num.cats)[which(names(table(pops.simple)) == "Uganda")]
    # Replace Uganda's color with black
    cols.hc.d[which(cols.hc.d == uganda.col)] = "#000000"

    # Set Kenya to grey
    kenya.col = getPal(num.cats)[which(names(table(pops.simple)) == "Kenya")]
    cols.hc.d[which(cols.hc.d == kenya.col)] = "#999999"
}

labels_colors(hc.d) = cols.hc.d

# Color branches
hc.d = color_branches(hc.d, k=length(hc.d.inds), h = 0, col=cols.hc.d)

# Get Ssese-specific info to add in color bar
geo.class.hc.d = do.call(rbind, lapply(hc.d.inds, function(x) {
    geo.class[which(hc$labels == x)] }))
geo.class.hc.d[!geo.class.hc.d %in% c("mainland", "island")] = NA
#                      island     mainland
col.geo.class.hc.d = c("#00BFC4", "#F8766D")[factor(geo.class.hc.d)]

# Color clusters in grey
#cluster.pal = sample(brewer.pal(length(levels(factor(hc.d.cut))), "Greys"))
cluster.cols = c("darkgrey", NA)[as.numeric(is.na(hc.d.cut)) + 1]

# Add bars to color individuals using snpEff functional data
if (do.snpeff) {
    if (length(snpeff.data) > 0) {
        snpeff.row.data = lapply(1:length(snpeff.data), function (snp.idx) {
            do.call(c, lapply(hc.d.inds, function (x) {
                # Remove _A or _B from individual name in hc.d.inds
                x = gsub("_[AB]", "", x)
                se = snpeff.data[[snp.idx]][[1]]
                se[se$ind == x,]$genotype
            }))
        })

        all.snpeff.row.data = do.call(cbind, snpeff.row.data)

        all.snpeff.row.data[all.snpeff.row.data == 1] = "lightblue"
        all.snpeff.row.data[all.snpeff.row.data == 2] = "darkblue"

        all.snpeff.row.data.clust = do.call(cbind,
            lapply(1:ncol(all.snpeff.row.data), function (x) {
                pasted.col = paste(all.snpeff.row.data[,x], cluster.cols)
                pasted.col = gsub("0 ", "", pasted.col)
                pasted.col = gsub("lightblue .*", "lightblue", pasted.col)
                pasted.col = gsub("darkblue .*",  "darkblue",  pasted.col)
            })
        )
    }
}

if (include.cluster.bar) {
    if (do.snpeff == FALSE || length(snpeff.data) == 0) {
        # If snpEff inclusion is off, this just includes cluster cols
        all.snpeff.row.data.clust = cluster.cols
    }
}

# ----------------------------------------------------------------------------------------
# --- Print out cluster membership info
# ----------------------------------------------------------------------------------------

sink(file=paste0(in.prefix, ".clusters.txt"))

for (cl in levels(factor(hc.d.cut))) {

    write(paste0("Cluster ", cl, ":"), stdout())
    this.clust.inds = names(hc.d.cut[which(hc.d.cut == cl)])

    this.clust.pops = unlist(lapply(this.clust.inds, function(ind) {
        pop.ss = as.character(ind.info[ind.info$V2 == ind, 3])
        pop.ag = paste(
            as.character(ind.info.ag1000g[ind.info.ag1000g$ox_code == ind,]$country),
            as.character(ind.info.ag1000g[ind.info.ag1000g$ox_code == ind,]$m_s),
                sep="_")
        if (length(pop.ss)) {
            if (pop.ss %in% c("BUGALAML", "BUWAMA", "KAZZI", "KIYINDI", "MITYANA")) {
                return ("LVB - MAINLAND")
            } else if (pop.ss %in% c("BANDA", "BUGALAIS", "BUKASA", "NSADZI", "SSERINYA")) {
                return ("LVB - ISLAND")
            } else {
                return(pop.ss)
            }
        } else if (length(pop.ag)) {
            return(pop.ag)
        } else if (length(grepl("SY", ind)) > 0) {
            return("SY")
        } else {
            return(ind)
        }
    }))

    for (cl.pop in names(table(this.clust.pops))) {
        write(print(paste0(cl.pop, " - ", sum(this.clust.pops == cl.pop))), stdout())
    }

}

sink()

# ----------------------------------------------------------------------------------------

# Gather up color bar data

col.bars = cbind(cols.hc.d, col.geo.class.hc.d)
# This may only contain cluster info, if snpEff inclusion is off
if (include.cluster.bar) {
    col.bars = cbind(col.bars, all.snpeff.row.data.clust)
}
num.col.bars = ncol(col.bars)

# Reverse color bars
col.bars = col.bars[,num.col.bars:1]

# Empty spaces are for dupliated rows
row.labels = c("Site", "Ugandan Site Type")
if (do.snpeff && length(snpeff.data) > 0) {
    row.labels = c(row.labels, do.call(c, lapply(snpeff.data, names)))
} else {
    if (include.cluster.bar) {
        row.labels = c(row.labels, "Clusters")
    }
}

# Replace warning about lack of stop codon with "?"
row.labels = gsub("WARNING_TRANSCRIPT_NO_START_CODON", "?", row.labels)

# Reverse row.labels
row.labels = rev(row.labels)

# Get number of individuals for placing legends
num.inds = length(cols.hc.d)

# ----------------------------------------------------------------------------------------

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

# ----------------------------------------------------------------------------------------

# Do plot
fancy.dendro.pdf = paste0(out.prefix, ".fancydendro.pdf")
pdf(fancy.dendro.pdf, height=12, width=20)

    par(mar = c(30,20,3,0.5))
    plot(hc.d,
        # Turn off title:
        # main = paste0(name, " - chr", chr, ":", start, "-", end),
        leaflab="none", yaxs='i', type="rectangle",
        yaxt='n', ann=FALSE)
    colored_bars(col.bars, hc.d,
        rowLabels = row.labels, sort_by_labels_order = FALSE, y_scale=2)

    # Add legend for pop colors
    par(xpd=TRUE)
    col.adj = getPal(num.cats)[1:num.cats]
    if (with.ag1000g) {
        col.adj[which(col.adj == uganda.col)] = "#000000"
        col.adj[which(col.adj == kenya.col)]  = "#999999"
        legend.ncol = 3
    } else {
        legend.ncol = 5
    }
    legend(0, -2.75, legend = make.italic(levels(categories)),
        fill = col.adj, ncol=legend.ncol,
        title="Site")

    # Add legend for LVB site type
    legend(0.46 * num.inds, -2.75, legend = c("Island", "Mainland"),
        fill = c("#00BFC4", "#F8766D"), ncol=2,
        title="Ugandan Site Type")

    # Add legend for clusters and SNPs
    if (include.cluster.bar) {
        if (do.snpeff) {
            snpeff.clust.legend = c("0/1 SNP", "1/1 SNP", "Cluster")
            snpeff.clust.fill = c("lightblue", "darkblue", "grey")
        } else {
            snpeff.clust.legend = c("Cluster")
            snpeff.clust.fill = c("grey")
        }
        legend(0.625 * num.inds, -2.75, legend = snpeff.clust.legend,
            fill = snpeff.clust.fill, ncol=2,
            title="Low Variation Clusters")
    }

dev.off()

write(paste("Wrote fancy dendrogram to", fancy.dendro.pdf), stderr())

# ----------------------------------------------------------------------------------------
# --- PCA
# ----------------------------------------------------------------------------------------

pca = glPca(gen, nf=10)

pops.pca = do.call(c, lapply(row.names(pca$scores), function(x) {
    as.character(ind.info[ind.info$V2 == x, 3])
}))

# scatter(pca, posi="topright")

pdf(paste0(out.prefix, ".pca.pdf"), height=20, width=20)
    plot(pca$scores[,1] ~ pca$scores[,2], col=pal[factor(geo.class)], pch=16,
        ylab="PC1", xlab="PC2")
dev.off()

# --- Plot PCA and tree with colors to represent position

pdf(paste0(out.prefix, ".pca.colorplot.pdf"), height=20, width=20)

    pca.cols = colorplot(pca$scores,pca$scores, transp=TRUE, cex=4)
    abline(h=0, v=0, col="grey")
    # add.scatter.eig(pca$eig[1:40],2,1,2, posi="topright", inset=.05, ratio=.3)

    tre = nj(dist(as.matrix(gen)))

    plot.new()

    par(mar=c(0, 0, 0, 0))
    plot(tre, no.margin=TRUE, direction='upwards', cex=0.5)
    tiplabels(pch=20, col=pca.cols, cex=1)

dev.off()

# dapc1 = dapc(gen, n.pca=10, n.da=1)

# loadingplot(abs(pca$loadings), threshold=quantile(abs(pca$loadings), 0.9995), axis=1)
