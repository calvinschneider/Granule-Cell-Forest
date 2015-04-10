#!/bin/bash
#
#$ -q som
#$ -cwd
#$ -pe openmp 1
#$ -j y
#$ -S /bin/bash
#$ -N dentate
#$ -o /dev/null
#$ -t 1-257
#
module load MATLAB/r2014a
directory="/fast-scratch/calvinjs"
mkdir -p "${directory}/Soma_Position"
mkdir -p "${directory}/Soma_Trimmed"
./dentate_2_p_somataposition $SGE_TASK_ID $directory
