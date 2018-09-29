#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Plot IBD sharing results
# ----------------------------------------------------------------------------------------

library(ggplot2)

#  FID1      Family ID for first individual
#  IID1      Individual ID for first individual
#  FID2      Family ID for second individual
#  IID2      Individual ID for second individual
#  RT        Relationship type given PED file
#  EZ        Expected IBD sharing given PED file
#  Z0        P(IBD=0)
#  Z1        P(IBD=1)
#  Z2        P(IBD=2)
#  PI_HAT    P(IBD=2)+0.5*P(IBD=1) ( proportion IBD )
#  PHE       Pairwise phenotypic code (1,0,-1 = AA, AU and UU pairs)
#  DST       IBS distance (IBS2 + 0.5*IBS1) / ( N SNP pairs )
#  PPC       IBS binomial test
#  RATIO     Of HETHET : IBS 0 SNPs (expected value is 2)

ibd = read.table("results/IBD/all.pass.snp.flt.noinv.genome", header=TRUE)

ibd.tmp = ibd
names(ibd.tmp)[1:4] = c("FID2", "IID2", "FID1", "IID1")
ibd.tmp = ibd.tmp[,c(3:4,1:2,5:ncol(ibd.tmp))]

ibd = rbind(ibd, ibd.tmp)

ibd$FID1 = factor(ibd$FID1, levels=unique(ibd$FID1))
ibd$FID2 = factor(ibd$FID2, levels=unique(ibd$FID1))

ggplot(ibd, aes(FID1, FID2, fill=PI_HAT)) + geom_tile()

# ----------------------------------------------------------------------------------------

# Bring in individual info

ind.info = read.csv("data/ssese_individual_info_bugala_split.csv")

# Bring in plate info to infer sequencing ID (e.g. ssese42.PE)
plates = read.table("data/sample_to_seq_id_mapping.txt")
names(plates) = c("mosquito.id", "seq.id")

# Reduce individual info to just sequenced animals
ind.info = ind.info[ind.info$mosquito_id %in% plates$mosquito.id,]

# ----------------------------------------------------------------------------------------

island   = as.character(ibd$IID1)
site     = as.character(ibd$IID1)
sex      = as.character(ibd$IID1)

for (i in 1:dim(ind.info)[1]) {

    this.ind = ind.info$mosquito_id[i]
    this.ind.info = ind.info[ind.info$mosquito_id == this.ind,]

    this.seq.id = plates$seq.id[plates$mosquito.id == this.ind]

    island = replace(island, island==this.seq.id, this.ind.info$island)
    site   = replace(site,   site==this.seq.id,   this.ind.info$site_name)
    sex    = replace(sex,    sex==this.seq.id,    this.ind.info$sex)
}

ibd$IID1.island = island
ibd$IID1.site   = site
ibd$IID1.sex    = sex

# ----------------------------------------------------------------------------------------

island   = as.character(ibd$IID2)
site     = as.character(ibd$IID2)
sex      = as.character(ibd$IID2)

for (i in 1:dim(ind.info)[1]) {

    this.ind = ind.info$mosquito_id[i]
    this.ind.info = ind.info[ind.info$mosquito_id == this.ind,]

    this.seq.id = plates$seq.id[plates$mosquito.id == this.ind]

    island = replace(island, island==this.seq.id, this.ind.info$island)
    site   = replace(site,   site==this.seq.id,   this.ind.info$site_name)
    sex    = replace(sex,    sex==this.seq.id,    this.ind.info$sex)
}

ibd$IID2.island = island
ibd$IID2.site   = site
ibd$IID2.sex    = sex

# ----------------------------------------------------------------------------------------

# Order sites by island / mainland
ibd$IID1.island[ibd$IID1.island == "BUGALAML"] = "BUGALA (M)"
ibd$IID1.island[ibd$IID1.island == "BUGALAIS"] = "BUGALA (I)"
ibd$IID2.island[ibd$IID2.island == "BUGALAML"] = "BUGALA (M)"
ibd$IID2.island[ibd$IID2.island == "BUGALAIS"] = "BUGALA (I)"

# Fix sites
ibd$IID1.island[ibd$IID1.island == "KAZZI"] = "KAAZI"
ibd$IID1.island[ibd$IID1.island == "MITYANA"] = "WAMALA"
ibd$IID2.island[ibd$IID2.island == "KAZZI"] = "KAAZI"
ibd$IID2.island[ibd$IID2.island == "MITYANA"] = "WAMALA"


is.sites = c("BANDA", "BUGALA", "BUGALA (I)", "BUKASA", "NSADZI", "SSERINYA")
ml.sites = c("BUGALA (M)", "BUWAMA", "KAAZI", "KIYINDI", "WAMALA")

ibd$IID1.is.island = ibd$IID1.island %in% is.sites
ibd$IID2.is.island = ibd$IID2.island %in% is.sites

site.labels = sort(unique(paste(!ibd$IID1.is.island, ibd$IID1.island, ibd$FID1)))
site.labels = gsub("(TRUE|FALSE) ", "", site.labels)
site.labels.short = gsub(" ssese.*", "", site.labels)

first = which(c(NA, site.labels.short) != c(site.labels.short, NA)) - 1
not.first = which(c(NA, site.labels.short) == c(site.labels.short, NA))
site.labels.short[not.first] = ""

last.island = first[length(ml.sites)]

simple.cap = function(x) {
    s = strsplit(x, " ")[[1]]
    lstr = paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
        sep="", collapse=" ")

    gsub("\\(i\\)", "(I)", gsub("\\(m\\)", "(M)", lstr))
}

site.labels.short = as.vector(sapply(site.labels.short, simple.cap))

num.labels = length(site.labels.short)

p = ggplot(ibd,
        aes(paste(!IID1.is.island, IID1.island, FID1),
            paste(!IID2.is.island, IID2.island, FID2),
            fill=log(PI_HAT))) +
    geom_tile() +
    geom_vline(xintercept=first + 0.5, color='grey') +
    geom_hline(yintercept=first + 0.5, color='grey') +
    geom_vline(xintercept=last.island + 0.5, color='black') +
    geom_hline(yintercept=last.island + 0.5, color='black') +
    scale_x_discrete(
        labels=c(rep("", 4),
                 site.labels.short[-c((num.labels-3):num.labels)]),
        position="top") +
    scale_y_discrete(labels=site.labels.short, position="right") +
    scale_fill_gradient(low='white', high='black', name="Proportion IBD",
        na.value='white') +
    xlab("") + ylab("") +
    theme_bw(base_size=20) +
    theme(axis.text.x = element_text(angle=90, hjust=0, vjust=1),
          axis.text.y = element_text(vjust=-0.5),
          panel.border = element_rect(colour = "black", fill=NA, size=3))

ggsave(p, filename="reports/ibd.pdf", height=7, width=9)

# ----------------------------------------------------------------------------------------

# --- Check for IBD sharing between BANDA and mainland sites

banda = ibd[(ibd$IID1.island == "BANDA" | ibd$IID2.island == "BANDA") & ibd$PI_HAT != 0,]

table(c(banda$IID1.island, banda$IID2.island))
