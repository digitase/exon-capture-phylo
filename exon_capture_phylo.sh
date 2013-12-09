#!/bin/bash
#$ -cwd
#$ -m aeb
#$ -M u5205339@anu.edu.au
#$ -r y
#$ -N job_exon_capture_phylo
#$ -t 1-1
#$ -tc 20 
#$ -R y
#$ -pe threads 1

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
source "test.config"

# Create outermost directory level for all script output
if [ ! -d "$OUT_DIR" ]; then
    mkdir "$OUT_DIR"
fi

# WARNING Clean directory
cd "$OUT_DIR"
rm ./* -r

# Use TASK_ID as an index to extract a sample out of library
libs=( $(cat $LIBRARIES_LIST) )
filnum=$[SGE_TASK_ID-1]
sample_name=${libs[$filnum]}
echo Performing $0 with SGE_TASK ID $SGE_TASK_ID on sample $sample_name

# 1. assembleByProt
# Args
    # 1. Sample name
    # 2. Samples directory
    # 3. Script output directory
    # 4. Target protein sequences
    # 5. Name for blastx database
# TODO I/O intensive: produces an alignment for each exon for each sample
# Possibly pipe
# TODO Where is the error.log file coming from?
echo assembleByProt with blastx database name "$TARGET_PROTEIN_BLAST_DB_NAME"
perl "$SCRIPT_DIR/pl/assembleByProtv2.pl" "$sample_name" "$LIBRARIES_DIR" "$OUT_DIR" \
                                          "$TARGET_PROTEIN_SEQS" "$TARGET_PROTEIN_BLAST_DB_NAME" \
    #1> "$sample_name.out" 2> "$sample_name.err"

# 2. callVelvetAssemblies
for k_value in ${VELVET_K_VALUES[@]}; do
    echo callVelvetAssemblies at "$k_value"
    perl "$SCRIPT_DIR/pl/callVelvetAssemblies.pl" "$sample_name" "$OUT_DIR" "$k_value" \
        #1>> "$sample_name.out" 2>> "$sample_name.err"
done

# 3. catcontigs
# TODO What is this?
# -l virtual_free=16G,h_vmem=20G
# -q bigmem.q
# TODO TARGET_PROTEIN_SEQS_DIR and TARGET_PROTEIN_SEQS are redundant
echo catContigs
export PATH=$PATH:"$CAP3_LOCATION"

perl "$SCRIPT_DIR/pl/catcontigs.pl" "$sample_name" "$OUT_DIR" \
                                    "$TARGET_PROTEIN_SEQS_LIST" "$TARGET_PROTEIN_SEQS_DIR" \
                                    $VELVET_K_VALUES \
    #1>> "$sample_name.out" 2>> "$sample_name.err"

exit

# 4. callBestContig
echo bestcontig_distrib
perl "$SCRIPT_DIR/pl/bestcontig_distrib.pl" "$sample_name" "$OUT_DIR" "$LIBRARIES_LIST" \
                                            "$TARGET_EXON_SEQS_DIR" "$TARGET_EXON_SEQS_LIST" \
                                            "$TARGET_PROTEIN_SEQS_DIR" "$TARGET_PROTEIN_SEQS_LIST" \
                                            "$ALL_PROTEIN_BLAST_DB_NAME" "$MIN_OVERLAP"
    #1>> "$sample_name.out" 2>> "$sample_name.err"


# 5. gatherContigs
perl "$SCRIPT_DIR/pl/gathercontigs.pl" \
    #1>> "$sample_name.out" 2>> "$sample_name.err"
