#!/bin/bash
# Use perldoc <script_name> for documentation on each pipeline component.
# AUTHOR: Ben Bai (u5205339@anu.edu.au)
# DATE: Nov 2013-Feb 2014

# Include the config file, which is in valid bash format
CONFIG_FILE="$1"
source $CONFIG_FILE

cd "$OUT_DIR"

echo ===== Started gather_exons.sh at $(date) =====

# 9. gatherContigs
echo gathercontigs at $(date)
perl "$SCRIPT_DIR/pl/gathercontigs.pl" "$OUT_DIR" "$SAMPLES_LIST" "$TARGET_PROTEIN_SEQS_LIST" \
                                       "$TARGET_EXON_SEQS_LIST"

# 10. gatherAmbigContigs
echo gatherAmbigContigs at $(date)
perl "$SCRIPT_DIR/pl/gatherAmbigContigs.pl" "$OUT_DIR" "$SAMPLES_LIST" "$TARGET_EXON_SEQS_LIST"

exit        
