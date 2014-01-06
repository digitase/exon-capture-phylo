#!/bin/bash
#$ -cwd
#$ -m aeb
#$ -M u5205339@anu.edu.au
#$ -r y
#$ -N job_exon_capture_phylo
#$ -R y
#$ -l virtual_free=10G,h_vmem=12G
#$ -q bigmem.q
# -pe threads 5
# -t 1-1
# -tc 20 

# Mail on job abortion/end
# Rerun on abort
# Use 4 threads and reserve slots as they free
# Job arrayed for each library

# Parameters will be in a config file: "source test.config"

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

# 5.1 gatherContigs
#     for each exon, collect the best exon from each sample

# Include the config file, which is in valid bash format
# source "test.config"
echo Started exon_capture_phylo at $(date)
CONFIG_FILE="$1"
source "$CONFIG_FILE"
echo Using config file at "$CONFIG_FILE"

# Create outermost directory level for all script output
if [ ! -d "$OUT_DIR" ]; then
    mkdir "$OUT_DIR"
fi

cd "$OUT_DIR" || exit 

# Use TASK_ID as an index to extract a sample out of library

libs=( $(cat $LIBRARIES_LIST) )
lib_num=${#libs[@]}

# assemble contigs
cat "$LIBRARIES_LIST" | xargs -n 1 --max-procs 20 -I {} "$SCRIPT_DIR/sh/assemble_exons.sh" {} "$SCRIPT_DIR/$CONFIG_FILE"

exit
# gather contigs by exon and by sample
"$SCRIPT_DIR/sh/gather_exons.sh" "$SCRIPT_DIR/$CONFIG_FILE"


# call variants
cat "$LIBRARIES_LIST" | xargs -n 1 --max-procs 20 -I {} "$SCRIPT_DIR/sh/call_variants.sh" {} "$SCRIPT_DIR/$CONFIG_FILE"

exit
