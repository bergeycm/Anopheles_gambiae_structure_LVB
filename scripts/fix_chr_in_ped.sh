#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Edit the MAP file and change 0 to the right chromosome index
# ----------------------------------------------------------------------------------------

MAP_FILE=$1

# Grab chromosome from input file
CHR=$(echo $MAP_FILE | sed -e "s/.*chr\(.*\)\.pass.*/\1/" | sed -e "s/\.missingtoref//g")

# Figure out chromosome index
CHR_INDEX=$(ls data/chr*.flt.vcf.gz | \
    grep -n "${CHR}" | cut -d":" -f1)

# Replace zero or chr name with chromosome index
sed -i -e "s/^[^\t]*\t/${CHR_INDEX}\t/" $MAP_FILE

exit
