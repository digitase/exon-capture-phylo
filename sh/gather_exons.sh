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

echo Started gather_exons.sh with at $(date)

# 5.1 gatherContigs
# TODO get rid of contig num file argument
echo gathercontigs at $(date)
perl "$SCRIPT_DIR/pl/gathercontigs.pl" "$OUT_DIR" "$LIBRARIES_LIST" "$TARGET_PROTEIN_SEQS_LIST" \
                                       "$TARGET_EXON_SEQS_LIST"

# 5.2 gather contigs from each exon by sample
echo best2ref at $(date)
perl "$SCRIPT_DIR/pl/best2refs.pl" "$OUT_DIR" "$LIBRARIES_LIST" "$TARGET_EXON_SEQS_LIST"

exit        
