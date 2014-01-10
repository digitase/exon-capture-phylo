#!/bin/bash

# 1. assembleByProt
    # blastx: transcriptome library on a list to Anolis targets 
    # collect reads hit to each target

# 2. callVelvetAssemblies
#     Velvet: assemble exons at 6 k-values

# 3. catContigs
#     CAP3: combine (perfectly) overlapping contigs
#     exonerate: extract exon sequence

# 4. callBestContig
#     blastx: exon sequences to Anolis targets
#     select one exon by reciprocal best hit

# Include the config file, which is in valid bash format
# source "test.config"
sample_name="$1"
CONFIG_FILE="$2"
source "$CONFIG_FILE"
echo Started assemble_exons.sh with "$sample_name" at $(date)

# 1. assembleByProt
# TODO I/O intensive: produces an alignment for each exon for each sample
# Possibly pipe
echo assembleByProt with blastx database name "$TARGET_PROTEIN_BLAST_DB_NAME"
perl "$SCRIPT_DIR/pl/assembleByProtv2.pl" "$sample_name" "$LIBRARIES_DIR" "$OUT_DIR" \
                                          "$TARGET_PROTEIN_SEQS" "$TARGET_PROTEIN_BLAST_DB_NAME" "$TARGET_PROTEIN_SEQS_LIST"

# 2. callVelvetAssemblies
echo callVelvetAssemblies at $(date)
printf "%s\n" "${VELVET_K_VALUES[@]}" | xargs -n 1 -I {} perl "$SCRIPT_DIR/pl/callVelvetAssemblies.pl" \
                                                                  "$sample_name" "$OUT_DIR" "$TARGET_PROTEIN_SEQS_LIST" {} 

# 3. catcontigs
echo catContigs at $(date)
export PATH=$PATH:"$CAP3_DIR"

perl "$SCRIPT_DIR/pl/catcontigs.pl" "$sample_name" "$OUT_DIR" \
                                    "$TARGET_PROTEIN_SEQS_LIST" "$TARGET_PROTEIN_SEQS" \
                                    ${VELVET_K_VALUES[@]}

# 4. callBestContig
echo bestcontig_distrib at $(date)
perl "$SCRIPT_DIR/pl/bestcontig_distrib.pl" "$sample_name" "$OUT_DIR" "$LIBRARIES_LIST" \
                                            "$TARGET_EXON_SEQS" "$TARGET_EXON_SEQS_LIST" \
                                            "$ALL_PROTEIN_SEQS" "$ALL_PROTEIN_BLAST_DB_NAME" \
                                            "$MIN_OVERLAP"

exit
