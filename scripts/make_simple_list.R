#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Make simple individual info file with sequence IDs.
# ----------------------------------------------------------------------------------------

options(stringsAsFactors = FALSE)

# Bring in plate info to infer sequencing ID (e.g. ssese42.PE)
plates = read.table("data/sample_to_seq_id_mapping.txt")
names(plates) = c("mosquito.id", "seq.id")

ind.info = read.csv("data/ssese_individual_info.csv")
plates.info = merge(plates, ind.info, by.x="mosquito.id", by.y="mosquito_id")

write.table(plates.info[,c(1:4, 6, 10)], file="data/ssese_individual_info_simple.txt", 
    quote=FALSE, row.names=FALSE, col.names=FALSE)
