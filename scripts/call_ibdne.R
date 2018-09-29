#!/usr/bin/env Rscript

# ========================================================================================
# --- Split inferred IBD by site and infer site pop history with IBDNe
# ========================================================================================

options(stringsAsFactors=FALSE)

args = commandArgs(trailingOnly=TRUE)
site = args[1]
chr  = args[2]    # E.g. 3

ibdne.jar = Sys.getenv("IBDNE_JAR")

# ----------------------------------------------------------------------------------------
# --- Function to split IBD file
# ----------------------------------------------------------------------------------------

ibd = rbind(
    read.table(paste0("results/chr", chr, "L.pass.snp.flt.ibdseq.ibd")),
    read.table(paste0("results/chr", chr, "R.pass.snp.flt.ibdseq.ibd"))
)

ind.info = read.table("data/ssese_individual_info_simple_bugala_split.txt")

sites = unique(ind.info$V3)

do.split = function (site) {

    site.inds = ind.info[ind.info$V3 == site, "V2"]
    site.ibd = ibd[which(ibd$V1 %in% site.inds & ibd$V3 %in% site.inds),]

    site.ibd
}

# ----------------------------------------------------------------------------------------
# --- Infer pop history with IBDNe
# ----------------------------------------------------------------------------------------

# --- Write IBD info to temporary file
ibd.tmp = tempfile(pattern = paste0("TMP.chr", chr, ".pass.snp.flt.ibdseq.", site, "."),
    tmpdir = tempdir(),
    fileext = ".ibd")

write.table(do.split(site), file=ibd.tmp,
    sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

ibd.cmd = paste0("cat ", ibd.tmp, " | ",
                      "java -jar ", ibdne.jar, " ",
                          "map=data/Ag_ALL.plink.map ",
                          "out=results/chr", chr, ".pass.snp.flt.", site, ".ne ",
                          "mincm=0.0005 ",
                          "minregion=0.01 ",
                          "trimcm=0.005 ",
                          "gmax=1000 ",
                          "filtersamples=false ",
                          "nthreads=8 ",
                          "npairs=0 ",
                          "nits=1000 ",
                          "nboots=100 ",
                          "seed=", sample(1:10000, 1))

#   Ag1000G params:
#       "minregion": 10.0
#       "filtersamples": "false"
#       "trim": 0.05
#       "minibd": 0.1

system(ibd.cmd)
