#!/usr/bin/env python

# Standard library imports
import os
import sys

# Local application imports
from misc import check_file_exists
from preprocessing import fastqc_quality, fastp_trimming
from pe_mapper import bwa_mapping, sam_to_index_bam
from bam_variant import picard_markdup, ivar_trim

def check_quality(r1_file, r2_file, output, name_folder, sub_folder_name, threads, sample):
    """
    Function that check reads quality with fastqc
    Output is located in output/name_folder/sub_folder_name
    """

    # Set name files
    out_dir = os.path.join(output, name_folder)                                 # folder
    out_subdir = os.path.join(out_dir, sub_folder_name)                         # subfolder

    # Check quality with fastqc
    fastqc_quality(r1_file, r2_file, out_subdir, threads)

def trim_read(r1_file, r2_file, sample, output, threads):
    """
    Function that trim reads considering some restrictions
    using fastp

    Output is located in output/Trimmed
    """

    # Set name files
    out_trim_dir = os.path.join(output, "Trimmed")                          # Folder
    out_trim_name_r1 = sample + ".trimmed_R1.fastq.gz"                      # r1 filename
    out_trim_name_r2 = sample + ".trimmed_R2.fastq.gz"                      # r2 filanema
    output_trimming_file_r1 = os.path.join(out_trim_dir, out_trim_name_r1)  # absolute path r1 filename
    output_trimming_file_r2 = os.path.join(out_trim_dir, out_trim_name_r2)  # absolute path r2 filename

    # Use fastp for trimming
    fastp_trimming(r1_file, r2_file, sample, out_trim_dir, threads=threads, 
        min_qual=20, window_size=10, min_len=35)
    
    # Return path of trimmed files
    return (output_trimming_file_r1, output_trimming_file_r2)

def mapping(output_trimming_file_r1, output_trimming_file_r2, sample, output, reference, threads):
    """
    Function that maps reads in reference genome using BWA
    Then, transform sam file to sorted bam file.

    Output is located in Bam directory
    """

    out_map_dir = os.path.join(output, "Bam")                       # folder
    out_map_name = sample + ".rg.sorted.bam"                        # file name
    output_map_file = os.path.join(out_map_dir, out_map_name)       # absolute path file name

        
    # Mapping to reference using BWA. As output a sam file is created
    bwa_mapping(output_trimming_file_r1, output_trimming_file_r2,
                reference, sample, out_map_dir, threads=threads)
    # Converts sam file to bam file using samtools and sort alingments
    sam_to_index_bam(sample, out_map_dir, output_trimming_file_r1, threads=threads)
    
    # Return absolute path mapped file
    return (output_map_file)

def mark_duplicates(output_map_file, sample, output):
    """
    Function that marks and remove duplicates using picard tools

    Sorted bam file is removed and just a sorted bam file without
    duplicates is kept. Output is located in Bam directory.
    """

    out_map_dir = os.path.join(output, "Bam")                           # Folder
    out_markdup_name = sample + ".rg.markdup.sorted.bam"                # Filename
    output_markdup_file = os.path.join(out_map_dir, out_markdup_name)   # absolute path file name


    # Mark and delete duplicates with picard tools
    picard_markdup(output_map_file) 
    
    # Return absolute path marked and remove duplicates file
    return (output_markdup_file)

def trim_ivar(output_markdup_file, sample, primers):
    """
    Function that trim reamining primers of aligned reads using ivar.
    Then uses samtools to sort reads and create an index to access 
    overlapping alignments quickly. 

    Output is located in Bam directory.
    """

    # Trimming of remaining primers is performed with ivar
    ivar_trim(output_markdup_file, primers, sample,
                min_length=30, min_quality=20, sliding_window_width=4)


def sample_mapper(output, r1_file, r2_file, sample, reference, threads, primers):
    """
    Function that maps reads to reference genome.
    """

    # INPUT ARGUMENTS
    ################
    check_file_exists(r1_file)
    check_file_exists(r2_file)

    # QUALITY CHECK in RAW with fastqc
    ######################################################
    # Check quality of input fastq with fastqc and store info in output/Quality/raw
    check_quality(r1_file, r2_file, output, "Quality", "raw", threads, sample)

    # QUALITY TRIMMING AND ADAPTER REMOVAL WITH fastp
    ###################################################
    # Trim reads by window. Trim regions or reads that does not satisfy min_qual=20, window_size=10, min_len=35
    # by using fastp. Output is in output/Trimmed
    output_trimming_file_r1, output_trimming_file_r2 = trim_read(r1_file, r2_file, sample, output, threads)

    # QUALITY CHECK in TRIMMED with fastqc
    ######################################################
    # Check quality of trimmed fastq with fastqc and store info in output/Quality/processed
    check_quality(output_trimming_file_r1, output_trimming_file_r2, output, "Quality", "processed", threads, sample)

    # MAPPING WITH BWA - SAM TO SORTED BAM - ADD HEADER SG
    #####################################################
    # Map trimmed reads to reference genome. As output a sorted bam file is generated in output/Bam
    output_map_file = mapping(output_trimming_file_r1, output_trimming_file_r2, sample, output, reference, threads)

    #MARK DUPLICATES WITH PICARDTOOLS ###################
    #####################################################
    # Mark and remove duplicates in bam file. Previous sorted bam file in output/Bam is sustituted
    # by a sorted without duplicates bam file.
    output_markdup_file = mark_duplicates(output_map_file, sample, output)

    #TRIM PRIMERS WITH ivar trim ########################
    #####################################################
    # Trim aligned reads (bam file) using ivar. Output file is in output/Bam.
    trim_ivar(output_markdup_file, sample, primers)

# Input variables
output = sys.argv[1]
r1_file = sys.argv[2]
r2_file = sys.argv[3]
sample = sys.argv[4]
reference = sys.argv[5]
threads = int(sys.argv[6])
primers = sys.argv[7]

# Function
sample_mapper(output, r1_file, r2_file, sample, reference, threads, primers)