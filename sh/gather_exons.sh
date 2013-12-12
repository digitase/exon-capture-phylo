#!/bin/bash
#$ -cwd
#$ -m aeb
#$ -M u5205339@anu.edu.au
#$ -r y
#$ -N job_gather_exons

# Include the config file, which is in valid bash format
# source "test.config"
CONFIG_FILE="$1"
source $CONFIG_FILE

# 5. gatherContigs
#     for each exon, collect the best exon from each sample

cd "$OUT_DIR"

# 5. gatherContigs
echo gathercontigs at $(date)
perl "$SCRIPT_DIR/pl/gathercontigs.pl" "$OUT_DIR" "$LIBRARIES_LIST" "$TARGET_PROTEIN_SEQS_LIST" \
                                       "$TARGET_EXON_SEQS_DIR" "$TARGET_EXON_SEQS_LIST" \
                                       "$CONTIG_NUM_FILE" \
    #1>> "$sample_name.out" 2>> "$sample_name.err"

