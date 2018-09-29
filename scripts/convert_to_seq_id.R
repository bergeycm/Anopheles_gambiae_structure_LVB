#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Script that converts mosquito ID lists (BBS-C-M1) to sequencing IDs (ssese21.PE)
# ----------------------------------------------------------------------------------------

options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly = TRUE)
list.in = args[1]
list.out = gsub("samples", "seqids", list.in)

if (file.size(list.in) == 0) {
    file.copy(list.in, list.out)
    write("No lines in input. Skipping file.", stderr())
    quit()
}

id.list = read.table(list.in)
names(id.list) = "mosquito.id"

# Bring in plate info to infer sequencing ID (e.g. ssese42.PE)
plates = read.table("data/sample_to_seq_id_mapping.txt")
names(plates) = c("mosquito.id", "seq.id")

list.matched = merge(id.list, plates)

write.table(list.matched$seq.id, file=list.out,
    quote=FALSE, row.names=FALSE, col.names=FALSE)
