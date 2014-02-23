#!/bin/bash
# Use perldoc <script_name> for documentation on each pipeline component.
# AUTHOR: Ben Bai (u5205339@anu.edu.au)
# DATE: Nov 2013-Feb 2014

# Include the config file, which is in valid bash format
CONFIG_FILE="$1"
source $CONFIG_FILE

cd "$OUT_DIR"

echo ===== Started prepare_data.sh at $(date) =====

# 0. prepare_data
echo prepareData at $(date)
perl "$SCRIPT_DIR/pl/prepareData.pl" "$OUT_DIR" \
                                     "$ALL_PROTEIN_SEQS" "$TARGET_PROTEIN_SEQS_LIST" \
                                     "$MAKEBLASTDB_PATH"

exit        
