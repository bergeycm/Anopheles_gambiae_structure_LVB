#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Create a raster copy of PDFs of XP-EHH results
# ----------------------------------------------------------------------------------------

for pdf in `ls reports/xpehh_and_fst_between_populations_*pdf`; do
    png=${pdf/pdf/png}
    if [[ ! -s $png ]]; then
        echo "Converting $pdf..."
        convert $pdf $png
    fi
done
