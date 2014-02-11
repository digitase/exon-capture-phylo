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
sample_name="$1"
CONFIG_FILE="$2"
source "$CONFIG_FILE"
echo ===== Started assemble_exons.sh with "$sample_name" at $(date) =====

# 1. assembleByProt
echo assembleByProt "$sample_name" at $(date)
perl "$SCRIPT_DIR/pl/assembleByProtv2.pl" "$sample_name" "$SAMPLES_DIR" "$OUT_DIR" \
                                          "$TARGET_PROTEIN_BLAST_DB_NAME" "$TARGET_PROTEIN_SEQS_LIST" \
                                          "$BLAST_PROCS_PER_SAMPLE" "$BLASTALL_PATH"

# 2. callVelvetAssemblies
echo callVelvetAssemblies "$sample_name" at $(date)
printf "%s\n" "${VELVET_K_VALUES[@]}" | xargs -n 1 -I {} perl "$SCRIPT_DIR/pl/callVelvetAssemblies.pl" \
                                                                  "$sample_name" "$OUT_DIR" "$TARGET_PROTEIN_SEQS_LIST" \
                                                                      "$VELVETH_PATH" "$VELVETG_PATH" {} 

# 3. catcontigs
echo catContigs "$sample_name" at $(date)
perl "$SCRIPT_DIR/pl/catcontigs.pl" "$sample_name" "$OUT_DIR" \
                                    "$TARGET_PROTEIN_SEQS_LIST" \
                                    "$CAP3_PATH" "$EXONERATE_PATH" \
                                    ${VELVET_K_VALUES[@]}

# 4. callBestContig
echo bestcontig_distrib "$sample_name" at $(date)
perl "$SCRIPT_DIR/pl/bestcontig_distrib.pl" "$sample_name" "$OUT_DIR" \
                                            "$ALL_EXON_SEQS" "$TARGET_EXON_SEQS_LIST" \
                                            "$ALL_PROTEIN_BLAST_DB_NAME" \
                                            "$MIN_OVERLAP" "$EXONERATE_PATH" "$BLASTALL_PATH"

# 5. gather best contigs from each exon
echo best2ref "$sample_name" at $(date)
perl "$SCRIPT_DIR/pl/best2refs.pl" "$sample_name" "$OUT_DIR" "$TARGET_EXON_SEQS_LIST"

exit
