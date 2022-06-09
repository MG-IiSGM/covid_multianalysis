#!/bin/bash

# Input variables
OUTDIR=$1
[ ! -d $OUTDIR ] && echo "Directory $OUTDIR DOES NOT exists." && exit 1
output_markdup_trimmed_file=$2
[ ! -f $output_markdup_trimmed_file ] && echo "$output_markdup_trimmed_file DOES NOT exists." && exit 1
SAMPLE=$3
REF=$4
[ ! -f $REF ] && echo "$REF DOES NOT exists." && exit 1
Annotation=$5
[ ! -f $Annotation ] && echo "$Annotation DOES NOT exists." && exit 1
min_quality=15
min_frequency_threshold=0.01
min_depth=1

#VARIANT CALLING WTIH ivar variants##################
#####################################################
# Variable
out_variant_dir=$(echo "$OUTDIR/Variants")
ivar_folder=$(echo "$out_variant_dir/ivar_raw")
prefix=$(echo "$ivar_folder/$SAMPLE")

samtools mpileup -aa -A -d 0 -B -Q 0 --reference $REF $output_markdup_trimmed_file | \
    ivar variants -p $prefix -q $min_quality -t $min_frequency_threshold \
    -m $min_depth -r $REF -g $Annotation
