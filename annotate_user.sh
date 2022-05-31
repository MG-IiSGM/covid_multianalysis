#!/bin/bash
#SBATCH --cpus-per-task=1

python /home/laura/Laura_intel/Desktop/covid_multianalysis/annotate_user.py $1 $2 $3 $4 $5 $6
