#!/usr/bin/env python


# Standard library imports
from distutils.cmd import Command
from distutils.command import check
import multiprocessing
import os
import sys
import re
import logging
import concurrent.futures

# Third party imports
import argparse
import subprocess
import datetime


# Local application imports
from misc import check_file_exists, extract_sample, check_create_dir, execute_subprocess, \
    extract_read_list, file_to_list, obtain_group_cov_stats, clean_unwanted_files, \
    check_reanalysis, vcf_stats, remove_low_quality, obtain_overal_stats
from preprocessing import fastqc_quality, fastp_trimming, format_html_image
from pe_mapper import bwa_mapping, sam_to_index_bam
from bam_variant import picard_dictionary, samtools_faidx, picard_markdup, ivar_trim, ivar_variants, ivar_consensus, \
    replace_consensus_header, create_bamstat, create_coverage, create_consensus
from vcf_process import filter_tsv_variants
from annotation import annotate_snpeff, annotate_pangolin, user_annotation, user_annotation_aa, annotation_to_html, \
    report_samples_html
from compare_snp import ddtb_add, ddtb_compare, ddbb_create_intermediate, revised_df, remove_position_range


# Local application imports
from covidma import check_quality, trim_read, mapping, mark_duplicates, trim_ivar, ivar_variant_calling \
    , variant_filtering, consensus_create, bamstats, coverage_stats

output = sys.argv[0]
args = sys.argv[1]
logger = sys.argv[2]
r1_file = sys.argv[3]
r2_file = sys.argv[4]
sample_list_F = sys.argv[5]
new_samples = sys.argv[6]
reference = sys.argv[7]


# Extract sample name
sample = extract_sample(r1_file, r2_file)
# Counter for new samples
new_sample_number = 0
# True if samples needs to be analysed
if sample in sample_list_F:

    sample_number = str(sample_list_F.index(sample) + 1)
    sample_total = str(len(sample_list_F))

    out_map_dir = os.path.join(output, "Bam")                                               # Folder
    out_markdup_trimmed_name = sample + ".rg.markdup.trimmed.sorted.bam"                    # Filname
    output_markdup_trimmed_file = os.path.join(out_map_dir, out_markdup_trimmed_name)       # absolute path to filename

    # Check if sample has been already analysed
    # If True, sample does not need to be analysed again
    if sample in new_samples:
        new_sample_number = str(int(new_sample_number) + 1)
        new_sample_total = str(len(new_samples))
        logger.info("\n" + WHITE_BG + "STARTING SAMPLE: " + sample +
                    " (" + sample_number + "/" + sample_total + ")" + " (" + new_sample_number + "/" + new_sample_total + ")" + END_FORMATTING)
    else:
        logger.info("\n" + WHITE_BG + "STARTING SAMPLE: " + sample +
                    " (" + sample_number + "/" + sample_total + ")" + END_FORMATTING)

    # True if not .rg.markdup.trimmed.sorted.bam exits
    # (trimming and mapping of alingned reads is already done)
    if not os.path.isfile(output_markdup_trimmed_file):
        
        # INPUT ARGUMENTS
        ################
        check_file_exists(r1_file)
        check_file_exists(r2_file)

        # QUALITY CHECK in RAW with fastqc
        ######################################################
        # Check quality of input fastq with fastqc and store info in output/Quality/raw
        check_quality(r1_file, r2_file, output, "Quality", "raw", logger, args, sample)

        # QUALITY TRIMMING AND ADAPTER REMOVAL WITH fastp
        ###################################################
        # Trim reads by window. Trim regions or reads that does not satisfy min_qual=20, window_size=10, min_len=35
        # by using fastp. Output is in output/Trimmed
        output_trimming_file_r1, output_trimming_file_r2 = trim_read(r1_file, r2_file, logger, sample, output, args)

        # QUALITY CHECK in TRIMMED with fastqc
        ######################################################
        # Check quality of trimmed fastq with fastqc and store info in output/Quality/processed
        check_quality(output_trimming_file_r1, output_trimming_file_r2, output, "Quality", "processed", logger, args, sample)

        # MAPPING WITH BWA - SAM TO SORTED BAM - ADD HEADER SG
        #####################################################
        # Map trimmed reads to reference genome. As output a sorted bam file is generated in output/Bam
        output_map_file = mapping(logger, output_trimming_file_r1, output_trimming_file_r2, sample, output, reference, args)

        #MARK DUPLICATES WITH PICARDTOOLS ###################
        #####################################################
        # Mark and remove duplicates in bam file. Previous sorted bam file in output/Bam is sustituted
        # by a sorted without duplicates bam file.
        output_markdup_file = mark_duplicates(logger, output_map_file, sample, output)

        #TRIM PRIMERS WITH ivar trim ########################
        #####################################################
        # Trim aligned reads (bam file) using ivar. Output file is in output/Bam.
        trim_ivar(logger, output_markdup_trimmed_file, output_markdup_file, sample, args)

    else:
        logger.info(YELLOW + DIM + output_markdup_trimmed_file +
                " EXIST\nOmmiting BAM mapping and BAM manipulation in sample " + sample + END_FORMATTING)

    ########################END OF MAPPING AND BAM MANIPULATION##########################
    #####################################################################################

    #VARIANT CALLING WTIH ivar variants##################
    #####################################################
    # Variant calling with ivar. Output is located in output/Variants/ivar_raw
    out_ivar_variant_file = ivar_variant_calling(logger, output, output_markdup_trimmed_file, sample, reference, annotation)

    #VARIANT FILTERING ##################################
    #####################################################
    # Filter variants detected with ivar. Output is located in output/Variants/ivar_filtered
    variant_filtering(output, sample, out_ivar_variant_file, logger)

    #CREATE CONSENSUS with ivar consensus##################
    #######################################################
    # Create a consensus fasta file from bam file trimmed and without duplicates.
    # Output is located in output/Consensus/ivar.
    consensus_create(output, sample, output_markdup_trimmed_file, logger)

    ########################CREATE STATS AND QUALITY FILTERS###############################
    #######################################################################################
    #CREATE Bamstats #######################################
    ########################################################
    # Compute metrics from trimmed and without duplicates bam file using samtools.
    # Output is located in output/Stats/Bamstats.
    bamstats(output, sample, output_markdup_trimmed_file, logger, args)

    #CREATE CoverageStats ##################################
    ########################################################
    coverage_stats(output, sample, output_markdup_trimmed_file, logger)