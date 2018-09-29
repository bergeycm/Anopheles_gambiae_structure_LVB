#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

# ----------------------------------------------------------------------------------------
# --- Make pop list for CLUMPAK
# ----------------------------------------------------------------------------------------

tmp = lapply(c(".", ".withM."), function (m.flag) {

    fam.file = paste0("data/ssese_with_ag1000g/adm_subsets/chr3", m.flag,
        "replicate1.LD.fam")

    adm.inds = data.frame(
        ind = read.table(fam.file)$V1
    )

    ind.info = read.table("data/ssese_individual_info_simple_bugala_split.txt")

    # Bring in plate info to infer sequencing ID (e.g. ssese42.PE)
    plates = read.table("data/sample_to_seq_id_mapping.txt")
    names(plates) = c("mosquito.id", "seq.id")

    # Reduce individual info to just sequenced animals
    ind.info = ind.info[ind.info$V1 %in% plates$mosquito.id,]

    ind.info.ag1000g = read.table("data/ag1000g.phase1.ar3/samples.all.txt",
        sep="\t", quote="", header=TRUE)

    pops = data.frame(rbind(cbind(ind.info$V2,              ind.info$V3),
                            cbind(ind.info.ag1000g$ox_code, ind.info.ag1000g$population))
           )
    names(pops) = c("ind", "pop")

    adm = merge(adm.inds, pops, sort=FALSE)

    adm$pop[adm$pop %in% c("BANDA", "BUKASA", "BUGALA", "BUGALAIS",
        "NSADZI", "SSERINYA")] = "LVB_SseseIslands"
    adm$pop[adm$pop %in% c("BUGALAML", "BUWAMA", "KAZZI", "KIYINDI", "MITYANA")] =
        "LVB_Mainland"

    out.file = paste0("data/ssese_with_ag1000g/adm_subsets/pop.list", m.flag, "txt")

    write.table(adm[,2], file=out.file,
        sep=" ", quote=FALSE, col.names=FALSE, row.names=FALSE)
})
