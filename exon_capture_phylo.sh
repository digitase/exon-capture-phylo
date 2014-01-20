#!/bin/bash
#$ -cwd
#$ -m aeb
#$ -M u5205339@anu.edu.au
#$ -r y
#$ -N job_exon_capture_phylo
#$ -R y
#$ -l virtual_free=10G,h_vmem=12G
#$ -q bigmem.q
#$ -tc 20 
#$ -t 1-2

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
echo ===== Started exon_capture_phylo at $(date) =====
CONFIG_FILE="$1"
source "$CONFIG_FILE"
echo Using config file at "$CONFIG_FILE"

# Create outermost directory level for all script output
if [ ! -d "$OUT_DIR" ]; then
    mkdir "$OUT_DIR"
fi

cd "$OUT_DIR" || exit 

# Check if we not using a job array
if [ -z "$SGE_TASK_ID" ] || [ "$SGE_TASK_ID" == "undefined" ]; then

    echo Using xargs to run $MAX_PARALLEL_SAMPLES samples in parallel 

    # prepare data
    "$SCRIPT_DIR/sh/prepare_data.sh" "$SCRIPT_DIR/$CONFIG_FILE"

    # assemble contigs
    < "$LIBRARIES_LIST" xargs -n 1 --max-procs "$MAX_PARALLEL_SAMPLES" -I {} "$SCRIPT_DIR/sh/assemble_exons.sh" "{}" "$SCRIPT_DIR/$CONFIG_FILE"

    # call variants
    < "$LIBRARIES_LIST" xargs -n 1 --max-procs "$MAX_PARALLEL_SAMPLES" -I {} "$SCRIPT_DIR/sh/call_variants.sh" "{}" "$SCRIPT_DIR/$CONFIG_FILE"

    # gather contigs by exon and by sample
    "$SCRIPT_DIR/sh/gather_exons.sh" "$SCRIPT_DIR/$CONFIG_FILE"

# We are using a job array
else

    readarray -t SAMPLE_NAMES_ARRAY < "$LIBRARIES_LIST"
    sample_num=${#SAMPLE_NAMES_ARRAY[@]}
    sample_name=${SAMPLE_NAMES_ARRAY[$SGE_TASK_ID - 1]}
    echo SGE Job Array SGE_TASK_ID=$SGE_TASK_ID detected. Using sample $sample_name out of $sample_num samples

    # data is already prepared continue
    if [ -d "$OUT_DIR/prepare_data_complete.lock" ]; then
        echo already prepped!
        :
    # else check if the data is being prepared. If not, prepare the data, then inform that the data is prepared
    elif mkdir "$OUT_DIR/preparing_data.lock"; then
        echo I will prep
        "$SCRIPT_DIR/sh/prepare_data.sh" "$SCRIPT_DIR/$CONFIG_FILE"
        mv "$OUT_DIR/preparing_data.lock/" "$OUT_DIR/prepare_data_complete.lock/"

    # else poll until prepare_data is complete
    else
        
        echo I will wait
        while ! [ -d "$OUT_DIR/prepare_data_complete.lock/" ]; do
            sleep 10
            echo waiting for prepare data complete
        done

    fi
    
    echo main part now

    # assemble contigs
    "$SCRIPT_DIR/sh/assemble_exons.sh" "$sample_name" "$SCRIPT_DIR/$CONFIG_FILE"

    # call variants
    "$SCRIPT_DIR/sh/call_variants.sh" "$sample_name" "$SCRIPT_DIR/$CONFIG_FILE"

    # inform that the sample is done
    echo I finished main part
    echo "$SGE_TASK_ID" >> "$OUT_DIR/completed_ids.txt"

    # if all samples are done, this job needs to gather
    if [ "$(wc -l < "$OUT_DIR/completed_ids.txt")" -eq "$sample_num" ]; then

        # lock it to make sure this is the only job gathering
        if mkdir "$OUT_DIR/gather_data.lock"; then
            echo I gather!
            # gather data 
            "$SCRIPT_DIR/sh/gather_exons.sh" "$SCRIPT_DIR/$CONFIG_FILE"
            
            # clean up
            rmdir "$OUT_DIR/prepare_data_complete.lock/"
            rm "$OUT_DIR/completed_ids.txt"
            rmdir "$OUT_DIR/gather_data.lock"
        else
            echo race condition avoided. heehee
        fi

    else
        echo I am not the last job. cya, not gathering
    fi

fi

echo ===== Finished at $(date) =====
exit
