#!/bin/bash

# Input variables
input_file=$1
[ ! -f $input_file ] && echo "$input_file DOES NOT exists." && exit 1
output_folder=$2
[ ! -d $output_folder ] && echo "Directory $output_folder DOES NOT exists." && exit 1
output_filename=$3
threads=$4
max_ambig=$5

# Annotate consensus fasta files to obtian sample linage using pangolin.
# Output is located in Annotation/pangolin.
pangolin $input_file --outdir $output_folder --outfile $output_filename \
    --threads $threads --max-ambig $max_ambig
