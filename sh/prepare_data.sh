#!/bin/bash

# Include the config file, which is in valid bash format
CONFIG_FILE="$1"
source $CONFIG_FILE

# 0. prepare_data
# Copy in files
# Prepare the blast dbs

cd "$OUT_DIR"
echo ===== Started prepare_data.sh at $(date) =====

# 0. prepare_data
echo prepareData at $(date)
perl "$SCRIPT_DIR/pl/prepareData.pl" "$OUT_DIR" \
                                     "$ALL_PROTEIN_SEQS" "$TARGET_PROTEIN_SEQS_LIST" \
                                     "$MAKEBLASTDB_PATH"

exit        
