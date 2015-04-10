#!/bin/bash
#
#$ -q som
#$ -cwd
#$ -l mem_free=16G   
#$ -pe openmp 1
#$ -j y
#$ -S /bin/bash
#$ -N dentate
#$ -o /dev/null
#
module load MATLAB/r2014a
mkdir -p Outputs/
./dentate_1_s_fillsomata
