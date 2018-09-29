#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Compute allele frequencies by island for genes of interest
# ----------------------------------------------------------------------------------------

options(stringsAsFactors=FALSE)

args = commandArgs(trailingOnly=TRUE)

in.vcf = args[1]    # E.g. "data/chrX.pass.snp.rdg.1mb.vcf"

samps = system(paste0("grep '^#CHROM' ", in.vcf, " | head -n1"), intern=TRUE)

inds = strsplit(samps, split="\t")[[1]][-c(1:10)]
inds = data.frame(IID=inds)

ind.info = read.table("data/ssese_individual_info_simple.txt")
names(ind.info) = c("seq.id", "IID", "island", "site", "sex", "species")

ind.info.simp = merge(inds, ind.info, by="IID", all=FALSE, sort=FALSE)


geno.counts = list()

con = file(vcf, "r")

line.ct = 0

while (TRUE) {
    line = readLines(con, n = 1)

    write(line, stderr())

    if (length(line) == 0) {
        break
    }
    info = strsplit(line, split="\t")[[1]]

    chr = info[1]
    pos = info[2]

    if (grepl("^#", chr)) {
        next
    }

    genos = info[-c(1:10)]
    genos = gsub(":.*", "", genos)

    geno.isl = data.frame(cbind(genos, ind.info.simp$island))
    geno.isl[,1] = factor(geno.isl[,1], levels=c("0/0", "0/1", "1/0", "1/1"))

    allele.1.ct = colSums(table(geno.isl) * c(0,1,1,2))
    allele.ct = colSums(table(geno.isl)) * 2

    allele.1.prop = allele.1.ct / allele.ct

    res = c(chr, pos, allele.1.prop)

    geno.counts[[pos]] = res
}

close(con)

geno.counts.df = as.data.frame(do.call(rbind, geno.counts))
