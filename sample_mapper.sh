#!/bin/bash
#SBATCH -J mapper
#SBATCH --cpus-per-task=16

python /home/laura/Laura_intel/Desktop/covid_multianalysis/sample_mapper.py $1 $2 $3 $4 $5 $6 $7