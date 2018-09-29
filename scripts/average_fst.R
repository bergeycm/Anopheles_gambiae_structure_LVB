#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Average Fst (within species) between islands
# ----------------------------------------------------------------------------------------

options(stringsAsFactors=FALSE)

island.sites   = c("BANDA", "BUGALA", "BUGALAIS",
                   "BUKASA", "KOOME", "NSADZI", "SSERINYA")
mainland.sites = c("BUGALAML", "BUWAMA", "ENTEBBE",
                   "KAZZI", "KIYINDI", "MITYANA")

all.sites = sort(c(island.sites, mainland.sites))

# Between island results are set up like this:
# fst.avg[[species]][[island_1]][[island_2]] = value
fst.avg.l = list()
fst.med.l = list()

# Should we compute Fst between inversions?
do.inversions = FALSE

# Print header
cat(paste("Group1", "Group2", "species", "filter", "stat", "value", sep="\t"), "\n")

# ----------------------------------------------------------------------------------------

het = read.table("data/heterochromatin.bed")

names(het) = c("chr", "start", "end")

# ----------------------------------------------------------------------------------------

inv = read.table("data/inversion_sites.bed")

names(inv) = c("chr", "start", "end", "name")

inv$inversion = gsub("_.*", "", inv$name)

# Make a guess at where the distal breakpoint of 2Rb is
distal.2Rb.guess = c("2R", 18575300, 18575300, "2Rb_dist", "2Rb")
inv = rbind(inv, distal.2Rb.guess)

inv$start = as.numeric(inv$start)
inv$end   = as.numeric(inv$end)

# ----------------------------------------------------------------------------------------

full.coord = list()

for (inversion in unique(inv$inversion)) {

    this.chr   = inv$chr  [inv$inversion == inversion][1]
    this.start = inv$start[inv$inversion == inversion]
    this.end   = inv$end  [inv$inversion == inversion]

    full.coord[[inversion]] = c(min(this.start, this.end),
                                max(this.start, this.end))
}

# ----------------------------------------------------------------------------------------
# --- Get within-species, between-island average Fst
# --- ignoring heterochromatic stuff and inversions in the case of AG
# ----------------------------------------------------------------------------------------

fst.island.files = grep("window", Sys.glob("results/chr*.*.weir.fst"),
    value=TRUE, invert=TRUE)

# Now, just use chromosome 3
fst.island.files = fst.island.files[grepl("3", fst.island.files)]

isl.pairs = unique(gsub("results/chr.*\\.([^\\.]+)\\.weir\\.fst",
    "\\1", fst.island.files))

fst.avg.l = list()
fst.med.l = list()

for (this.isl.pair in isl.pairs) {

    islands = strsplit(this.isl.pair, "_")[[1]]

    if (! islands[1] %in% names(fst.avg.l)) {
        fst.avg.l[[islands[1]]] = list()
        fst.med.l[[islands[1]]] = list()
    }

    # Now, just use chromosome 3
    fst.files = Sys.glob(paste0("results/chr3*.",
        this.isl.pair, ".weir.fst"))

    fst = do.call(rbind, lapply(fst.files,
        function(x) { read.table(x, header=TRUE) }))

    # --- Mask out heterochromatic stuff

    fst$is.het = FALSE

    for (this.het.idx in 1:nrow(het)) {

        this.het = het[this.het.idx,]

        het.bool = fst$CHROM == this.het$chr &
            fst$POS >= this.het$start &
            fst$POS <= this.het$end

        if (sum(het.bool) > 0) {
            fst[het.bool,]$is.het = TRUE
        }
    }

    # --- Mask out inversion regions

    fst$inv.state = "collinear"

    for (this.inv.idx in 1:length(full.coord)) {

        this.inv = full.coord[this.inv.idx]
        this.chr = substr(names(this.inv), start=0, stop=2)

        inv.bool = fst$CHROM == this.chr &
            fst$POS >= this.inv[[1]][1] &
            fst$POS <= this.inv[[1]][2]

        if (sum(inv.bool) > 0) {
            fst[inv.bool,]$inv.state = names(this.inv)
        }
    }

    # --- Get mean of island pair Fst without heterochromatic stuff
    # --- and inversion regions

    bool = fst$inv.state == "collinear" & (!fst$is.het)
    fst.avg = mean(fst[bool,]$WEIR_AND_COCKERHAM_FST,   na.rm=TRUE)
    fst.med = median(fst[bool,]$WEIR_AND_COCKERHAM_FST, na.rm=TRUE)

    cat(paste(islands[1], islands[2], "AG", "no_het_or_inv", "mean",
        fst.avg, sep="\t"), "\n")
    cat(paste(islands[1], islands[2], "AG", "no_het_or_inv", "median",
        fst.med, sep="\t"), "\n")

    fst.avg.l[[islands[1]]][[islands[2]]] = fst.avg
    fst.med.l[[islands[1]]][[islands[2]]] = fst.med

    # --- Also compute average Fst in inversions, for kicks.
    # --- Skip if indicated.
    for (this.inv in unique(inv$inversion)) {

        if (do.inversions) {
            inv.bool = fst$inv.state == this.inv & (!fst$is.het)
            fst.avg = mean(fst[inv.bool,]$WEIR_AND_COCKERHAM_FST,   na.rm=TRUE)
            fst.med = median(fst[inv.bool,]$WEIR_AND_COCKERHAM_FST, na.rm=TRUE)

            cat(paste(islands[1], islands[2], "AG",
                paste0("no_het_", this.inv), "mean",
                fst.avg, sep="\t"), "\n")
            cat(paste(islands[1], islands[2], "AG",
                paste0("no_het_", this.inv), "median",
                fst.med, sep="\t"), "\n")
        }
    }
}

# ----------------------------------------------------------------------------------------
# --- Compute crappy estimate of Nm from Fst for island pairs
# ----------------------------------------------------------------------------------------

# Nm = (1 - Fst) / 4Fst

filter = "no_het_no_inv"

for (isl1 in names(fst.avg.l)) {

    for (isl2 in names(fst.avg.l[[isl1]])) {

        this.fst = fst.avg.l[[isl1]][[isl2]]
        Nm = (1 - this.fst) / (4 * this.fst)

        cat(paste(isl1, isl2, "AG", filter, "Nm", Nm, sep="\t"), "\n")
    }
}
