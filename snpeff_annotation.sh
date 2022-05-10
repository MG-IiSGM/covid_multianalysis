#!/bin/bash
#SBATCH -J dpu_snpeff
#SBATCH --cpus-per-task=1

python /home/laura/Laura_intel/Desktop/covid_multianalysis $1 $2 $3 $4 $5
