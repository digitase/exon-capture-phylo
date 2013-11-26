#!/bin/bash

# Parameters will be in a config file: "source test.config"

# Scripts ordering

# Job arrayed for each library

# 1. assembleByProt
#     blastx: transcriptome library on a list to Anolis targets 
#     collect reads hit to each target

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
    
source test.config

