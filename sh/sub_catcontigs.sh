#!/bin/bash
#$ -cwd
#$ -m ae
#$ -r y
#$ -N catcontigs 
#$ -t 1-100:1
#$ -tc 12
#$ -l virtual_free=16G,h_vmem=20G
#$ -q bigmem.q


DIR=/home2/jgb/camaxlibs/
libs=( $(cat ${DIR}camaxalllibs.txt) )
filnum=$[SGE_TASK_ID-1]
parfil=${libs[$filnum]}

echo $SGE_TASK_ID
echo $parfil

cd ${DIR}
nohup perl /home2/jgb/camaxlibs/scripts/catcontigs.pl ${parfil} 1> ${parfil}.err 2> ${parfil}.out
