#!/bin/bash

# 6. mapsnp
     # Use bowtie2 to map original library reads onto assemblies

# 7. gatkSNPcalls
#    call SNPs

# Requires
# -l virtual_free=10G,h_vmem=12G
# -q bigmem.q

# Include the config file, which is in valid bash format
# source "test.config"
sample_name="$1"
CONFIG_FILE="$2"
source "$CONFIG_FILE"
echo Started call_variants.sh with "$sample_name" at $(date)

echo mapsnp at $(date) 
perl "$SCRIPT_DIR/pl/mapsnp.pl" "$sample_name" "$LIBRARIES_DIR" "$OUT_DIR"

exit

echo gatkSNPcalls at $(date)
perl "$SCRIPT_DIR/pl/gatkSNPcalls.pl" "$sample_name" "$LIBRARIES_DIR" "$OUT_DIR" \
                                      "$PICARD_DIR" "$GATK_DIR"

exit
