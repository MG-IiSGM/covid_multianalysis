#!/bin/bash
#SBATCH -J dpu
#SBATCH --cpus-per-task=16

python /home/laura/Laura_intel/Desktop/covid_multianalysis/map_sample.py $1 $2 $3 $4 $5 $6 $7 $8 $9
