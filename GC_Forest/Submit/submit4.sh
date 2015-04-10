#!/bin/bash
#
#$ -q som
#$ -cwd
#$ -pe openmp 1
#$ -j y
#$ -S /bin/bash
#$ -N dentate
#$ -o /dev/null
#$ -t 1-136
#
module load MATLAB/r2014a
sublayer=1
directory="/fast-scratch/calvinjs"
mkdir -p "${directory}/All_Points/"
mkdir -p "${directory}/All_Points/GCL"
mkdir -p "${directory}/All_Points/IML"
mkdir -p "${directory}/All_Points/MML"
mkdir -p "${directory}/All_Points/OML"
mkdir -p "${directory}/All_Points/OOML"
./dentate_4_p_findallpoints $sublayer $SGE_TASK_ID $directory
