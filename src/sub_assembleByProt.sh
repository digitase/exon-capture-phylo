#!/bin/bash
#$ -cwd
#$ -m ae
#$ -r y
#$ -N assemble
#$ -t 1-100:1
#$ -tc 20 

# job array, throttled to 20
# jobs have same id num
# read library list into array

DIR=/home2/jgb/camaxlibs/

# list of libraries
libs=( $(cat ${DIR}camaxalllibs.txt) )

# task id: numbers from 1-100
filnum=$[SGE_TASK_ID-1]
parfil=${libs[$filnum]}

echo $SGE_TASK_ID
echo $parfil

cd ${DIR}
nohup perl /home2/jgb/camaxlibs/scripts/assembleByProtv2.pl ${parfil} 1> ${parfil}.err 2>${parfil}.out
