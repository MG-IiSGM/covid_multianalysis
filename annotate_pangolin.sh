#!/bin/bash
#SBATCH -J pango
#SBATCH --cpus-per-task=2

python /home/laura/Laura_intel/Desktop/covid_multianalysis/annotate_pangolin.py $1 $2 $3 $4 $5
