#!/bin/bash
#SBATCH -J useraa
#SBATCH --cpus-per-task=1

python /home/laura/Laura_intel/Desktop/covid_multianalysis/useraa_annotation.py $1 $2 $3 $4 $5
