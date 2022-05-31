#!/bin/bash
#SBATCH --cpus-per-task=1

python /home/laura/Laura_intel/Desktop/covid_multianalysis/variant_filtering.py $1 $2 $3