#!/bin/bash

# SGE configuration
#$ -cwd
#$ -M u5205339@anu.edu.au
#$ -m bea
#$ -N job_exon_capture_phylo
#$ -R y
#$ -l virtual_free=6G,h_vmem=8G
#$ -tc 4
#$ -t 1-2

echo ===== Started exon_capture_phylo at $(date) =====

# Include the config file, which is in valid bash format
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

    echo Using xargs to run $XARGS_PARALLEL_SAMPLES samples in parallel 

    # prepare data
    "$SCRIPT_DIR/sh/prepare_data.sh" "$SCRIPT_DIR/$CONFIG_FILE"

    # assemble contigs
    < "$SAMPLES_LIST" xargs -n 1 --max-procs "$XARGS_PARALLEL_SAMPLES" -I {} "$SCRIPT_DIR/sh/assemble_exons.sh" "{}" "$SCRIPT_DIR/$CONFIG_FILE"

    # call variants
    < "$SAMPLES_LIST" xargs -n 1 --max-procs "$XARGS_PARALLEL_SAMPLES" -I {} "$SCRIPT_DIR/sh/call_variants.sh" "{}" "$SCRIPT_DIR/$CONFIG_FILE"

    # gather contigs by exon and by sample
    "$SCRIPT_DIR/sh/gather_exons.sh" "$SCRIPT_DIR/$CONFIG_FILE"

# We are using a job array
else

    readarray -t SAMPLE_NAMES_ARRAY < "$SAMPLES_LIST"
    sample_num=${#SAMPLE_NAMES_ARRAY[@]}
    sample_name=${SAMPLE_NAMES_ARRAY[$SGE_TASK_ID - 1]}
    host_name="$(hostname)"
    echo SGE Job Array SGE_TASK_ID="$SGE_TASK_ID" detected. Processing sample "$sample_name" out of "$sample_num" samples on "$host_name"

    # data is already prepared continue
    if [ -d "$OUT_DIR/prepare_data_complete.lock" ]; then
        echo Detected prepared data.
    # else check if the data is being prepared. If not, prepare the data, then inform that the data is prepared
    elif mkdir -m 400 "$OUT_DIR/preparing_data.lock" 2>/dev/null; then
        echo Detected that data is unprepared.
        "$SCRIPT_DIR/sh/prepare_data.sh" "$SCRIPT_DIR/$CONFIG_FILE"
        mv "$OUT_DIR/preparing_data.lock/" "$OUT_DIR/prepare_data_complete.lock/"

    # else poll until prepare_data is complete
    else
        
        while ! [ -d "$OUT_DIR/prepare_data_complete.lock/" ]; do
            sleep 10
            echo Waiting for prepare data completion.
        done

    fi

    # assemble contigs
    "$SCRIPT_DIR/sh/assemble_exons.sh" "$sample_name" "$SCRIPT_DIR/$CONFIG_FILE"

    # call variants
    "$SCRIPT_DIR/sh/call_variants.sh" "$sample_name" "$SCRIPT_DIR/$CONFIG_FILE"

    # inform that the sample is done
    echo "$SGE_TASK_ID" >> "$OUT_DIR/completed_ids.txt"

    # if all samples are done, this job needs to gather
    num_samples_complete="$(wc -l < "$OUT_DIR/completed_ids.txt")"
    if [ "$num_samples_complete" -eq "$sample_num" ]; then

        # lock it to make sure this is the only job gathering
        if mkdir -m 400 "$OUT_DIR/gather_data.lock" 2>/dev/null; then
            echo Final sample complete. Gathering output on "$host_name"
            # gather data 
            "$SCRIPT_DIR/sh/gather_exons.sh" "$SCRIPT_DIR/$CONFIG_FILE"
            
            # clean up
            echo Cleaning up and exiting.
            rmdir "$OUT_DIR/prepare_data_complete.lock/"
            rm "$OUT_DIR/completed_ids.txt"
            rmdir "$OUT_DIR/gather_data.lock"
        else
            # Race condition avoided
            :
        fi

    else
        echo Sample "$sample_name" complete on "$host_name". "$num_samples_complete" of "$sample_num" samples completed. Exiting.
    fi

fi

echo ===== Finished at $(date) =====
exit
