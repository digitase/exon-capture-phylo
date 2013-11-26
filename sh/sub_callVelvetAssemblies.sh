#!/bin/bash
#$ -cwd
#$ -m ae
#$ -r y
#$ -N assemble
#$ -t 1-100:1
#$ -tc 24 

DIR=/home2/jgb/camaxlibs/
libs=( $(cat ${DIR}camaxalllibs.txt) )
filnum=$[SGE_TASK_ID-1]
parfil=${libs[$filnum]}

# loop over other 5 k values
k=81

echo $SGE_TASK_ID
echo $parfil

cd ${DIR}
nohup perl /home2/jgb/camaxlibs/scripts/callVelvetAssemblies.pl ${parfil} ${k} 1> ${parfil}.err 2>${parfil}.out
