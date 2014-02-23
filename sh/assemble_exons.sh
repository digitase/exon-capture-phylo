#!/bin/bash
# Use perldoc <script_name> for documentation on each pipeline component.
# AUTHOR: Ben Bai (u5205339@anu.edu.au)
# DATE: Nov 2013-Feb 2014

# Include the config file, which is in valid bash format
sample_name="$1"
CONFIG_FILE="$2"
source "$CONFIG_FILE"

echo ===== Started assemble_exons.sh with "[$sample_name]" at $(date) =====

# 1. assembleByProt
echo assembleByProt "[$sample_name]" at $(date)
perl "$SCRIPT_DIR/pl/assembleByProtv2.pl" "$sample_name" "$SAMPLES_DIR" "$OUT_DIR" \
                                          "$FORWARD_READS_SUFFIX" "$REVERSE_READS_SUFFIX" "$UNPAIRED_READS_SUFFIX" \
                                          "$TARGET_PROTEIN_SEQS_LIST" \
                                          "$BLAST_PROCS_PER_SAMPLE" "$BLASTALL_PATH" "$BLAST_EVALUE"

# 2. callVelvetAssemblies
echo callVelvetAssemblies "[$sample_name]" at $(date)
printf "%s\n" "${VELVET_K_VALUES[@]}" | xargs -n 1 -I {} perl "$SCRIPT_DIR/pl/callVelvetAssemblies.pl" \
                                                              "$sample_name" "$OUT_DIR" "$TARGET_PROTEIN_SEQS_LIST" \
                                                              "$VELVETH_PATH" "$VELVETG_PATH" {} 

# 3. catcontigs
echo catContigs "[$sample_name]" at $(date)
perl "$SCRIPT_DIR/pl/catcontigs.pl" "$sample_name" "$OUT_DIR" \
                                    "$TARGET_PROTEIN_SEQS_LIST" \
                                    "$CAP3_PATH" "$EXONERATE_PATH" \
                                    ${VELVET_K_VALUES[@]}

# 4. bestcontig_distrib
echo bestcontig_distrib "[$sample_name]" at $(date)
perl "$SCRIPT_DIR/pl/bestcontig_distrib.pl" "$sample_name" "$OUT_DIR" \
                                            "$ALL_EXON_SEQS" "$TARGET_EXON_SEQS_LIST" \
                                            "$MIN_OVERLAP" "$EXONERATE_PATH" "$BLASTALL_PATH"

# 5. best2ref
echo best2ref "[$sample_name]" at $(date)
perl "$SCRIPT_DIR/pl/best2refs.pl" "$sample_name" "$OUT_DIR" "$TARGET_EXON_SEQS_LIST"

exit
