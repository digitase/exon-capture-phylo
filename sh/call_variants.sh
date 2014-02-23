#!/bin/bash
# Use perldoc <script_name> for documentation on each pipeline component.
# AUTHOR: Ben Bai (u5205339@anu.edu.au)
# DATE: Nov 2013-Feb 2014

# Include the config file, which is in valid bash format
sample_name="$1"
CONFIG_FILE="$2"
source "$CONFIG_FILE"

echo ===== Started call_variants.sh with "[$sample_name]" at $(date) =====

# 6. mapsnp
echo mapsnp "[$sample_name]" at $(date)
perl "$SCRIPT_DIR/pl/mapsnp.pl" "$sample_name" "$SAMPLES_DIR" "$OUT_DIR" \
                                "$FORWARD_READS_SUFFIX" "$REVERSE_READS_SUFFIX" "$UNPAIRED_READS_SUFFIX" \
                                "$BOWTIE2_BUILD_PATH" "$BOWTIE2_PATH" "$SAMTOOLS_PATH"

# 7. gatkSNPcalls
echo gatkSNPcalls "[$sample_name]" at $(date)
perl "$SCRIPT_DIR/pl/gatkSNPcalls.pl" "$sample_name" "$SAMPLES_DIR" "$OUT_DIR" "$PICARD_DIR" "$GATK_DIR" "$JAVA_MAX_HEAP_SIZE" "$SAMTOOLS_PATH"

# 8. vcf2ambigfasta
echo vcf2ambigfasta "[$sample_name]" at $(date)
perl "$SCRIPT_DIR/pl/vcf2ambigfasta.pl" "$sample_name" "$OUT_DIR"

exit
