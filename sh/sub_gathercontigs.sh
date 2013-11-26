#!/bin/bash
#$ -cwd
#$ -m ae
#$ -r y
#$ -N gathercontigs 

DIR=/home2/jgb/assemble/crypto/


cd ${DIR}
nohup perl /home2/jgb/assemble/crypto/scripts/gathercontigs.pl  1>gathercontigs.err 2>gathercontigs.out
