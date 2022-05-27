#!/bin/bash
#SBATCH -J mapper5
#SBATCH --cpus-per-task=16

python /home/laura/Laura_intel/Desktop/covid_multianalysis/bamstats.py $1 $2 $3 $4