#!/bin/bash
#$ -cwd
#$ -m aeb
#$ -M u5205339@anu.edu.au
#$ -r y
#$ -N job_gather_exons

# Include the config file, which is in valid bash format
CONFIG_FILE="$1"
source $CONFIG_FILE

# 5. gatherContigs
#     for each exon, collect the best exon from each sample

cd "$OUT_DIR"

echo ===== Started gather_exons.sh at $(date) =====

# 5.1 gatherContigs
echo gathercontigs at $(date)
perl "$SCRIPT_DIR/pl/gathercontigs.pl" "$OUT_DIR" "$SAMPLES_LIST" "$TARGET_PROTEIN_SEQS_LIST" \
                                       "$TARGET_EXON_SEQS_LIST"

echo gatherAmbigContigs at $(date)
perl "$SCRIPT_DIR/pl/gatherAmbigContigs.pl" "$OUT_DIR" "$SAMPLES_LIST" "$TARGET_EXON_SEQS_LIST"

exit        
