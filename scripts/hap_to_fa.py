#!/usr/bin/env python

# ----------------------------------------------------------------------------------------
# --- Create FASTA files from phased haplotypes
# ----------------------------------------------------------------------------------------

# From: http://genometoolbox.blogspot.com/2013/07/create-fasta-sequences-for-phased.html

# Import and Create Files
import sys
haps_file=sys.argv[1]
sample_file=haps_file.strip(".haps")+".sample"
fasta_file=haps_file.strip(".haps")+".fasta"

# Make Dictionary from Haps file with RS Numbers as Keys
rs_dict={}
snps=[]
haps=open(haps_file)
haps=haps.readlines()
for i in range(len(haps)):
	col=haps[i].strip().split()
	if col[3] in ["A","C","G","T"] and col[4] in ["A","C","G","T"]:
		snp_id = col[0] + ":" + col[2]
		for j in range(5,len(col)):
			if col[j]=="0":
				col[j]=col[3]
			elif col[j]=="1":
				col[j]=col[4]
		rs_dict[snp_id]=col[5:]
		snps.append(snp_id)

# Make Dictionary with Sample ID as Keys
id_dict={}
sample=open(sample_file)
sample=sample.readlines()

ids=[]
for i in range(2,len(sample)):
	col=sample[i].strip().split()
	ids.append(col[1]+"_A")
	ids.append(col[1]+"_B")

for id in ids:
	id_dict[id]=[]

for snp in snps:
	for id in ids:
		id_dict[id].append(rs_dict[snp][ids.index(id)])

# Output FASTA file
fasta=open(fasta_file,"w")

for id in ids:
	print >> fasta, ">"+id
	print >> fasta, "".join(id_dict[id])

fasta.close()
