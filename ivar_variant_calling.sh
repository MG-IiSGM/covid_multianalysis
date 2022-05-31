#!/bin/bash
#SBATCH --cpus-per-task=1

python /home/laura/Laura_intel/Desktop/covid_multianalysis/ivar_variant_calling.py $1 $2 $3 $4 $5