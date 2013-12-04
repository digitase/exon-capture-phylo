#!/bin/bash
#$ -cwd
#$ -m aeb
#$ -M u5205339@anu.edu.au
#$ -r y
#$ -N job_exon_capture_phylo
#$ -t 1-1
#$ -tc 20 
#$ -pe threads 4
#$ -R y

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
# TODO Switch the config file path to a command line argument 
source "test.config"

# Create outermost directory level for all script output
if [ ! -d "$OUT_DIR" ]; then
    mkdir "$OUT_DIR"
fi

cd "$OUT_DIR"

# 1. assembleByProt

# Use TASK_ID as an index to extract a sample out of library
libs=( $(cat $SAMPLE_NAMES) )
filnum=$[SGE_TASK_ID-1]
parfil=${libs[$filnum]}

echo Performing $0 with SGE_TASK ID $SGE_TASK_ID on sample $parfil

# Args
# 1. Sample name
# 2. Samples directory
# 3. Script output directory
# 4. Target protein sequences
# 5. Name for blastx database
perl $SCRIPT_DIR/pl/assembleByProtv2.pl "$parfil" "$SAMPLES_DIR" "$OUT_DIR" "$TARGET_SEQS" "$BLAST_DB_NAME" \
# Create output logs
    1> $parfil.out 2> $parfil.err

perl $SCRIPT_DIR/pl/assembleByProtv2.pl ${parfil} $SAMPLES_DIR $OUT_DIR $BLAST_DB

