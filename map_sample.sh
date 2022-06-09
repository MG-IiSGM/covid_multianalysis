#!/bin/bash

# Input variables
OUTDIR=$1
[ ! -d $OUTDIR ] && echo "Directory $OUTDIR DOES NOT exists." && exit 1
R1=$2
[ ! -f $R1 ] && echo "$R1 DOES NOT exists." && exit 1
R2=$3
[ ! -f $R2 ] && echo "$R2 DOES NOT exists." && exit 1
SAMPLE=$4
REF=$5
[ ! -f $REF ] && echo "$REF DOES NOT exists." && exit 1
THREADS=$6
PRIMERS=$7
[ ! -f $PRIMERS ] && echo "$PRIMERS DOES NOT exists." && exit 1
min_qual=20
window_size=10
min_len=35
ivar_min_len=30
sliding_window_width=4

# QUALITY CHECK in RAW with fastqc
######################################################
# Variables
out_qc_dir=$(echo "$OUTDIR/Quality")
out_qc_pre_dir=$(echo "$out_qc_dir/raw")

# Check quality of input fastq with fastqc and store info in output/Quality/raw
fastqc $R1 $R2 -o $out_qc_pre_dir --threads $THREADS

# QUALITY TRIMMING AND ADAPTER REMOVAL WITH fastp
###################################################
# Variables
out_trim_dir=$(echo "$OUTDIR/Trimmed")
output_trimmed_r1=$(echo "$out_trim_dir/$SAMPLE.trimmed_R1.fastq.gz")
output_trimmed_r2=$(echo "$out_trim_dir/$SAMPLE.trimmed_R2.fastq.gz")
[ ! -d $(echo "$out_trim_dir/html") ] && mkdir $(echo "$out_trim_dir/html")
[ ! -d $(echo "$out_trim_dir/json") ] && mkdir $(echo "$out_trim_dir/json")
json_file=$(echo "$out_trim_dir/json/$SAMPLE""_fastp.json")
html_file=$(echo "$out_trim_dir/html/$SAMPLE""_fastp.html")

# Trim reads by window. Trim regions or reads that does not satisfy min_qual=20, 
# window_size=10, min_len=35 by using fastp. Output is in output/Trimmed
fastp --in1 $R1 --in2 $R2 --out1 $output_trimmed_r1 --out2 $output_trimmed_r2 \
    --detect_adapter_for_pe --cut_tail --cut_window_size $window_size \
    --cut_mean_quality $min_qual --length_required $min_len --json $json_file \
    --html $html_file --thread $THREADS

# QUALITY CHECK in TRIMMED with fastqc
######################################################
# Variables
out_qc_post_dir=$(echo "$out_qc_dir/processed")

# Check quality of trimmed fastq with fastqc and store info in output/Quality/processed
fastqc $output_trimmed_r1 $output_trimmed_r2 -o $out_qc_post_dir --threads $THREADS

# MAPPING WITH BWA - SAM TO SORTED BAM - ADD HEADER SG
#####################################################
# Variables
out_map_dir=$(echo "$OUTDIR/Bam")
sam_file=$(echo "$out_map_dir/$SAMPLE.sam")
bam_file=$(echo "$out_map_dir/$SAMPLE.bam")
sorted_bam_file=$(echo "$out_map_dir/$SAMPLE.sorted.bam")
rg_sorted_bam_file=$(echo "$out_map_dir/$SAMPLE.rg.sorted.bam")

# Map trimmed reads to reference genome. As output a sorted bam file is 
# generated in output/Bam
bwa index $REF
bwa mem -Y -M -t $THREADS -o $sam_file $REF $output_trimmed_r1 $output_trimmed_r2
samtools view -Sb $sam_file --threads $THREADS -o $bam_file && rm $sam_file
samtools sort $bam_file -o $sorted_bam_file && rm $bam_file
python /home/laura/Laura_intel/Desktop/covid_multianalysis/add_SG.py $SAMPLE $sorted_bam_file $rg_sorted_bam_file $output_trimmed_r1 \
    && rm $sorted_bam_file

#MARK DUPLICATES WITH PICARDTOOLS ###################
#####################################################
# Variables
rg_markdup_bam_file=$(echo "$out_map_dir/$SAMPLE.rg.markdup.bam")
rg_markdup_sorted_bam_file=$(echo "$out_map_dir/$SAMPLE.rg.markdup.sorted.bam")
[ ! -d $(echo "$out_map_dir/Stats") ] && mkdir $(echo "$out_map_dir/Stats")
stats_file=$(echo "$out_map_dir/Stats/$SAMPLE.markdup.metrics.txt")

# Mark and remove duplicates in bam file. Previous sorted bam file in output/Bam
# is sustituted by a sorted without duplicates bam file.
picard MarkDuplicates -I $rg_sorted_bam_file -O $rg_markdup_bam_file -M $stats_file \
    && rm $rg_sorted_bam_file
samtools sort $rg_markdup_bam_file -o $rg_markdup_sorted_bam_file \
    && rm $rg_markdup_bam_file

#TRIM PRIMERS WITH ivar trim ########################
#####################################################
# Variable
input_bai=$(echo "$rg_markdup_sorted_bam_file.bai")
prefix=$(echo "$out_map_dir/$SAMPLE.rg.markdup.trimmed")
output_trimmed_bam=$(echo "$prefix.bam")
rg_markdup_trimmed_sorted_bam_file=$(echo "$out_map_dir/$SAMPLE.rg.markdup.trimmed.sorted.bam")

# Trim aligned reads (bam file) using ivar. Output file is in output/Bam.
ivar trim -i $rg_markdup_sorted_bam_file -b $PRIMERS -p $prefix \
    -m $ivar_min_len -q $min_qual -s $sliding_window_width -e \
    && rm $rg_markdup_sorted_bam_file
samtools sort $output_trimmed_bam -o $rg_markdup_trimmed_sorted_bam_file \
    && rm $output_trimmed_bam
samtools index $rg_markdup_trimmed_sorted_bam_file && rm $input_bai
