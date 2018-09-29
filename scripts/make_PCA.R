#!/usr/bin/Rscript

source("http://bioconductor.org/biocLite.R")
#biocLite("gdsfmt")
#biocLite("SNPRelate")

library("gdsfmt")
library("SNPRelate")

args = commandArgs(trailingOnly = TRUE)

bed.in = args[1]
fam.in = gsub(".bed", ".fam", bed.in)
bim.in = gsub(".bed", ".bim", bed.in)

# Convert
gds.in = gsub(".bed", ".gds", bed.in)
snpgdsBED2GDS(bed.in, fam.in, bim.in, gds.in)

# ----------------------------------------------------------------------------------------

# Open the GDS file
(geno = snpgdsOpen(gds.in))

# Get sample id
sample.id = read.gdsn(index.gdsn(geno, "sample.id"))

# Get population
ind.info = read.csv("data/ssese_individual_info.csv")

# Bring in plate info to infer sequencing ID (e.g. ssese42.PE)
plate1 = read.csv("data/ssese_plate1_individuals.csv", header=FALSE)
plate2 = read.csv("data/ssese_plate2_individuals.csv", header=FALSE)
plates = rbind(plate1, plate2)
plates$seq.id = paste0("ssese", 1:nrow(plates), ".PE")
names(plates)[1] = "mosquito.id"

ind.info.id = merge(ind.info, plates, by.x="mosquito_id", by.y="mosquito.id")
mtch.order = match(ind.info.id$seq.id, sample.id)
#mtch.order = mtch.order[!is.na(mtch.order)]
ind.info.id = ind.info.id[order(mtch.order),]
ind.info.id = head(ind.info.id, n=length(sample.id))

stopifnot(sample.id == ind.info.id$seq.id)
pop_code = ind.info$site_name

# (Assumes the order of sample IDs and of population codes are the same)

# ----------------------------------------------------------------------------------------

plot.pca = function (species) {

    # Figure out AG and AR sets to do PCA separately
    this.sp.inds = ind.info.id[ind.info.id$ID.PCR.Result == species,]$seq.id
    
    # ----------------------------------------------------------------------------------------
    
    # Do PCA (Autosome-only by default)
    pca = snpgdsPCA(geno, num.thread=12, sample.id=this.sp.inds)
    
    # Variance proportion (%)
    pc.percent = pca$varprop*100
    head(round(pc.percent, 2))
    
    # ----------------------------------------------------------------------------------------
    
    # Make a data.frame
    pca.results =  data.frame(sample.id = pca$sample.id,
        pop = factor(pop_code)[match(pca$sample.id, sample.id)],
        EV1 = pca$eigenvect[,1],    # the first eigenvector
        EV2 = pca$eigenvect[,2],    # the second eigenvector
        stringsAsFactors = FALSE)
    
    # ----------------------------------------------------------------------------------------
    
    # Plot PCA
    
    pdf.out = gsub("data", "reports", bed.in)
    pdf.out = gsub(".bed", paste0(".PCA.", species, ".pdf"), pdf.out)
    
    pdf(file=pdf.out)
    	plot(pca.results$EV2, pca.results$EV1, col=as.integer(pca.results$pop), 
    			xlab="Eigenvector 2", ylab="Eigenvector 1")
    	legend("bottomright", legend=levels(pca.results$pop), 
    			pch="o", col=1:nlevels(pca.results$pop))
    dev.off()
    
    # Plot first 4 PCs
    
    multi.pdf.out = gsub(".pdf", ".multi.pdf", pdf.out)
    
    pdf(file=multi.pdf.out)
    	lbls = paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
    	pairs(pca$eigenvect[,1:4], col=pca.results$pop, labels=lbls)
    dev.off()

}

lapply(c("AG", "AR"), plot.pca)
