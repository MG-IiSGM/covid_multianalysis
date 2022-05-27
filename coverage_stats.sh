#!/bin/bash
#SBATCH -J mapper6
#SBATCH --cpus-per-task=1

python /home/laura/Laura_intel/Desktop/covid_multianalysis/coverage_stats.py $1 $2 $3