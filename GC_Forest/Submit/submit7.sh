#!/bin/bash
#
#$ -q som
#$ -cwd
#$ -pe openmp 1
#$ -j y
#$ -S /bin/bash
#$ -N dentate
#$ -o /dev/null
#$ -t 1-1186
#
module load MATLAB/r2014a
directory="/fast-scratch/calvinjs"
mkdir -p "${directory}/Trees_Jittered"
./dentate_7_p_jitter $SGE_TASK_ID $directory
