#!/bin/bash
#
#$ -q som
#$ -cwd
#$ -pe openmp 1
#$ -j y
#$ -S /bin/bash
#$ -N dentate
#$ -o /dev/null
#
module load MATLAB/r2014a
directory="/fast-scratch/calvinjs"
./dentate_5_s_choosepoints $directory
