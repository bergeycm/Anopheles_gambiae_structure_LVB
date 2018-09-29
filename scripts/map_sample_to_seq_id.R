#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Script that converts mosquito ID lists (BBS-C-M1) to sequencing IDs (ssese21.PE)
# ----------------------------------------------------------------------------------------

options(stringsAsFactors = FALSE)

# Bring in plate info to infer sequencing ID (e.g. ssese42.PE)
plate1 = read.csv("data/ssese_plate1_individuals.csv", header=FALSE)
plate2 = read.csv("data/ssese_plate2_individuals.csv", header=FALSE)
plates = rbind(plate1, plate2)
plates$seq.id = paste0("ssese", 1:nrow(plates), ".PE")
names(plates)[1] = "mosquito.id"

write.table(plates, file="data/sample_to_seq_id_mapping.txt",
    quote=FALSE, row.names=FALSE, col.names=FALSE)
