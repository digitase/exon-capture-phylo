#!/bin/bash

# 6. mapsnp
     # Use bowtie2 to map original library reads onto assemblies

# 7. gatkSNPcalls
#    call SNPs

# Requires
# -l virtual_free=10G,h_vmem=12G
# -q bigmem.q

# Include the config file, which is in valid bash format
sample_name="$1"
CONFIG_FILE="$2"
source "$CONFIG_FILE"
echo ===== Started call_variants.sh with "[$sample_name]" at $(date) =====

echo mapsnp "[$sample_name]" at $(date)
perl "$SCRIPT_DIR/pl/mapsnp.pl" "$sample_name" "$SAMPLES_DIR" "$OUT_DIR" \
                                "$FORWARD_READS_SUFFIX" "$REVERSE_READS_SUFFIX" "$UNPAIRED_READS_SUFFIX" \
                                "$BOWTIE2_BUILD_PATH" "$BOWTIE2_PATH" "$SAMTOOLS_PATH"

echo gatkSNPcalls "[$sample_name]" at $(date)
perl "$SCRIPT_DIR/pl/gatkSNPcalls.pl" "$sample_name" "$SAMPLES_DIR" "$OUT_DIR" "$PICARD_DIR" "$GATK_DIR" "$JAVA_MAX_HEAP_SIZE" "$SAMTOOLS_PATH"

echo vcf2ambigfasta "[$sample_name]" at $(date)
perl "$SCRIPT_DIR/pl/vcf2ambigfasta.pl" "$sample_name" "$OUT_DIR"

exit
