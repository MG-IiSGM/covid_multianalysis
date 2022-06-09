#!/bin/bash

# Input variables
OUTDIR=$1
[ ! -d $OUTDIR ] && echo "Directory $OUTDIR DOES NOT exists." && exit 1
SAMPLE=$2
output_markdup_trimmed_file=$3
[ ! -f $output_markdup_trimmed_file ] && echo "$output_markdup_trimmed_file DOES NOT exists." && exit 1
THREADS=$4

#CREATE Bamstats #######################################
########################################################
# Variable
out_stats_dir=$(echo "$OUTDIR/Stats")
out_stats_bamstats_dir=$(echo "$out_stats_dir/Bamstats")
output_file=$(echo "$out_stats_bamstats_dir/$SAMPLE.bamstats")

samtools flagstat --threads $THREADS $output_markdup_trimmed_file > $output_file

