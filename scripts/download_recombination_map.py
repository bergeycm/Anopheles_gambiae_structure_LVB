#!/usr/bin/env python

# ----------------------------------------------------------------------------------------
# --- Download recombination maps and convert to PLINK format
# ----------------------------------------------------------------------------------------

# Stolen from https://github.com/malariagen/ag1000g/blob/
# 3c13e0bf95f4e262901ac5cbbc3fd9d7b36c9ecf/notebooks/njh/
# 2017-03-15.run_ibdseq.ipynb

import pandas as pd
import numpy as np
import os

chroms = ['2L', '2R', '3L', '3R', 'X']

# --- Download original recombination maps

for chrom in chroms:
    cmd = "wget ftp://ngs.sanger.ac.uk/production/ag1000g/phase2/AR1/recombination_maps/"
    cmd = cmd + "Ag_" + chrom + ".map -O data/Ag_" + chrom + ".map"
    os.system(cmd)

# --- Convert to PLINK format

def convert_map_to_plink(chrom, outf):

    map_df = pd.read_csv("data/Ag_" + chrom + ".map", sep="\t")

    tdf = map_df[["dist", "pos"]]
    tdf.insert(0, "varlab", np.repeat(".", tdf.shape[0]))
    tdf.insert(0, "chrom", np.repeat(chrom, tdf.shape[0]))

    tdf.to_csv(outf, sep=" ", header=False, index=False)

for chrom in chroms:
    convert_map_to_plink(chrom, "data/Ag_" + chrom + ".plink.map")

# --- Combine all maps

os.system("rm -f data/Ag_ALL.plink.map")
os.system("cat data/Ag_*.plink.map > data/Ag_ALL.plink.map")
