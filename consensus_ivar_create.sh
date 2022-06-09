#!/bin/bash

# Input variable
OUTDIR=$1
[ ! -d $OUTDIR ] && echo "Directory $OUTDIR DOES NOT exists." && exit 1
SAMPLE=$2
output_markdup_trimmed_file=$3
[ ! -f $output_markdup_trimmed_file ] && echo "$output_markdup_trimmed_file DOES NOT exists." && exit 1
min_quality=20
min_frequency_threshold=0.8
min_depth=20
uncovered_character='N'

#CREATE CONSENSUS with ivar consensus##################
#######################################################
# Variable
out_consensus_dir=$(echo "$OUTDIR/Consensus")
out_consensus_ivar_dir=$(echo "$out_consensus_dir/ivar")
out_ivar_consensus_file=$(echo "$out_consensus_ivar_dir/$SAMPLE.fa")
prefix=$(echo "$out_consensus_ivar_dir/$SAMPLE")

samtools mpileup -aa -A -d 0 -B -Q 0  $output_markdup_trimmed_file | \
    ivar consensus -p $prefix -q $min_quality -t $min_frequency_threshold \
    -m $min_depth -n $uncovered_character

# Replace consensus header
Sequence=$(cat $out_ivar_consensus_file | grep -v ">")
echo ">$SAMPLE" > $out_ivar_consensus_file \
    && echo $Sequence >> $out_ivar_consensus_file