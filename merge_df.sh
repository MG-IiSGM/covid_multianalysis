#!/bin/bash
#SBATCH -J merge_df
#SBATCH --cpus-per-task=1

python /home/laura/Laura_intel/Desktop/covid_multianalysis/merge_df.py $1 $2 $3 $4 $5 $6 $7 $8 $9 $10