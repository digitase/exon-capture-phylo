#!/bin/bash
#$ -cwd
#$ -m aeb
#$ -M u5205339@anu.edu.au
#$ -r y
#$ -N job_exon_capture_phylo
#$ -R y
#$ -pe orte 8
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

# 5. gatherContigs
#     for each exon, collect the best exon from each sample

# Include the config file, which is in valid bash format
# source "test.config"
echo Started at $(date)
CONFIG_FILE="$1"
source "$CONFIG_FILE"

# Create outermost directory level for all script output
if [ ! -d "$OUT_DIR" ]; then
    mkdir "$OUT_DIR"
fi

# WARNING Clean directory
cd "$OUT_DIR" || exit 
# rm "$OUT_DIR/*" -r

# Use TASK_ID as an index to extract a sample out of library

libs=( $(cat $LIBRARIES_LIST) )
lib_num=${#libs[@]}

# TODO this needs to be throttled to 20 jobs
# http://superuser.com/questions/345447/how-can-i-trigger-a-notification-when-a-job-process-ends
# or xargs
# Using 1-indexing to retain compatability with SGE_TASK_ID
for (( i=1; i<=$lib_num; i++ )); do
    filnum=$((i-1))
    sample_name=${libs[$filnum]}
    echo Performing assemble_exons with SGE_TASK_ID $filnum on sample $sample_name at $(date)
    xargs
    "$SCRIPT_DIR/sh/assemble_exons.sh" "$sample_name" "$SCRIPT_DIR/$CONFIG_FILE" &
done

# This needs to execute after all the array jobs are done
wait

"$SCRIPT_DIR/sh/gather_exons.sh" "$SCRIPT_DIR/$CONFIG_FILE"

