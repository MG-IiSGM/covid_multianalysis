#!/usr/bin/env python

# Standard library imports
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

"""
=============================================================
HEADER
=============================================================

INSTITUTION:IiSGM
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
d^v^b
VERSION=0.1
CREATED: 22 Sep 2020
REVISION:


TODO:
    Adapt check_reanalysis
    Check file with multiple arguments
    Check program is installed (dependencies)
================================================================
END_OF_HEADER
================================================================
"""

# COLORS AND AND FORMATTING

END_FORMATTING = '\033[0m'
WHITE_BG = '\033[0;30;47m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
MAGENTA = '\033[35m'
BLUE = '\033[34m'
CYAN = '\033[36m'
YELLOW = '\033[93m'
DIM = '\033[2m'

logger = logging.getLogger()


def main():
    """
    Create main function to capture code errors: https://stackoverflow.com/questions/6234405/logging-uncaught-exceptions-in-python
    """

    # ARGUMENTS

    def get_arguments():

        parser = argparse.ArgumentParser(
            prog='covidma.py', description='Pipeline to call variants (SNVs) with any non model organism. Specialised in SARS-CoV-2')

        input_group = parser.add_argument_group('Input', 'Input parameters')

        input_group.add_argument('-i', '--input', dest="input_dir", metavar="input_directory",
                                 type=str, required=True, help='REQUIRED.Input directory containing all fast[aq] files')
        input_group.add_argument('-r', '--reference', metavar="reference",
                                 type=str, required=True, help='REQUIRED. File to map against')
        input_group.add_argument('-a', '--annotation', metavar="annotation",
                                 type=str, required=True, help='REQUIRED. gff3 file to annotate variants')
        input_group.add_argument('-s', '--sample', metavar="sample", type=str,
                                 required=False, help='Sample to identify further files')
        input_group.add_argument('-L', '--sample_list', type=str, required=False,
                                 help='Sample names to analyse only in the file supplied')
        input_group.add_argument('-p', '--primers', type=str, default='/home/laura/DATABASES/Anotacion/COVID/primers/nCoV-2019.bed',
                                 required=False, help='Bed file including primers to trim')

        quality_group = parser.add_argument_group(
            'Quality parameters', 'parameters for diferent triming conditions')

        # quality_group.add_argument('-c', '--coverage20', type=int, default=90, required=False,
        #                            help='Minimum percentage of coverage at 20x to clasify as uncovered (Default 90)')
        # quality_group.add_argument('-n', '--min_snp', type=int, required=False,
        #                            default=1, help='SNP number to pass quality threshold')
        quality_group.add_argument('-c', '--min_coverage', type=int, default=70, required=False,
                                    help='Minimum percentage of coverage at 30x to clasify as uncovered (Default 70)')
        quality_group.add_argument('-u', '--max_unmap', type=int, default=10, required=False,
                                    help='Maximum percentage of unmaped to clasify as uncovered (Default 10)')

        output_group = parser.add_argument_group(
            'Output', 'Required parameter to output results')

        output_group.add_argument('-o', '--output', type=str, required=True,
                                  help='REQUIRED. Output directory to extract all results')
        output_group.add_argument('-C', '--noclean', required=False,
                                  action='store_false', help='Clean unwanted files for standard execution')

        params_group = parser.add_argument_group(
            'Parameters', 'parameters for diferent stringent conditions')

        params_group.add_argument('-T', '--threads', type=str, dest="threads",
                                  required=False, default=16, help='Threads to use')
        params_group.add_argument('-M', '--memory', type=str, dest="memory",
                                  required=False, default=32, help='Max memory to use')

        annot_group = parser.add_argument_group(
            'Annotation', 'parameters for variant annotation')

        annot_group.add_argument('-B', '--annot_bed', type=str, default=[],
                                 required=False, action='append', help='bed file to annotate')
        annot_group.add_argument('-V', '--annot_vcf', type=str, default=[],
                                 required=False, action='append', help='vcf file to annotate')
        annot_group.add_argument('-A', '--annot_aa', type=str, default=[],
                                 required=False, action='append', help='aminoacid file to annotate')
        annot_group.add_argument('-R', '--remove_bed', type=str, default=False,
                                 required=False, help='BED file with positions to remove')
        annot_group.add_argument('--mash_database', type=str, required=False,
                                 default=False, help='MASH ncbi annotation containing all species database')
        annot_group.add_argument('--snpeff_database', type=str, required=False,
                                 default='NC_045512.2', help='snpEFF annotation database')

        compare_group = parser.add_argument_group(
            'Compare', 'parameters for compare_snp')

        compare_group.add_argument('-S', '--only_snp', required=False,
                                   action='store_true', help='Use INDELS while comparing')

        arguments = parser.parse_args()

        return arguments

    args = get_arguments()

    ######################################################################
    #####################START PIPELINE###################################
    ######################################################################
    output = os.path.abspath(args.output)
    group_name = output.split("/")[-1]
    reference = os.path.abspath(args.reference)
    annotation = os.path.abspath(args.annotation)

    # LOGGING
    # Create log file with date and time
    right_now = str(datetime.datetime.now())
    right_now_full = "_".join(right_now.split(" "))
    log_filename = group_name + "_" + right_now_full + ".log"
    log_folder = os.path.join(output, 'Logs')
    check_create_dir(log_folder)
    log_full_path = os.path.join(log_folder, log_filename)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s:%(message)s')

    file_handler = logging.FileHandler(log_full_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    # stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    logger.info("\n\n" + BLUE + BOLD +
                "STARTING PIPELINE IN GROUP: " + group_name + END_FORMATTING)

    today = str(datetime.date.today())

    logger.info("ARGUMENTS:")
    logger.info(str(args))

    # Obtain all R1 and R2 from folder
    r1, r2 = extract_read_list(args.input_dir)

    # Check if there are samples to filter out
    sample_list_F = []
    if args.sample_list == None:
        logger.info("\n" + "No samples to filter")
        for r1_file, r2_file in zip(r1, r2):
            sample = extract_sample(r1_file, r2_file)
            sample_list_F.append(sample)
    else:
        logger.info("samples will be filtered")
        sample_list_F = file_to_list(args.sample_list)

    new_samples = check_reanalysis(args.output, sample_list_F)

    logger.info("\n%d samples will be analysed: %s" %
                (len(new_samples), ",".join(new_samples)))

    #PREPARE REFERENCE FOR MAPPING + FAI + DICT #########
    #####################################################

    # picard_dictionary(args)
    samtools_faidx(args)

    #DECLARE FOLDERS CREATED IN PIPELINE ################
    #AND KEY FILES ######################################
    #####################################################
    # Annotation related parameters
    # script_dir = os.path.dirname(os.path.realpath(__file__))

    # Output related
    out_qc_dir = os.path.join(output, "Quality")
    out_qc_pre_dir = os.path.join(out_qc_dir, "raw")  # subfolder
    out_qc_post_dir = os.path.join(out_qc_dir, "processed")  # subfolder
    out_trim_dir = os.path.join(output, "Trimmed")
    out_map_dir = os.path.join(output, "Bam")
    out_variant_dir = os.path.join(output, "Variants")
    out_variant_ivar_dir = os.path.join(
        out_variant_dir, "ivar_raw")  # subfolder
    out_filtered_ivar_dir = os.path.join(
        out_variant_dir, "ivar_filtered")  # subfolder
    out_consensus_dir = os.path.join(output, "Consensus")
    out_consensus_ivar_dir = os.path.join(
        out_consensus_dir, "ivar")  # subfolder

    out_stats_dir = os.path.join(output, "Stats")
    out_stats_bamstats_dir = os.path.join(
        out_stats_dir, "Bamstats")  # subfolder
    out_stats_coverage_dir = os.path.join(
        out_stats_dir, "Coverage")  # subfolder
    out_compare_dir = os.path.join(output, "Compare")

    out_annot_dir = os.path.join(output, "Annotation")
    out_annot_snpeff_dir = os.path.join(out_annot_dir, "snpeff")  # subfolder
    out_annot_pangolin_dir = os.path.join(
        out_annot_dir, "pangolin")  # subfolder
    out_annot_user_dir = os.path.join(out_annot_dir, "user")  # subfolder
    out_annot_user_aa_dir = os.path.join(out_annot_dir, "user_aa")  # subfolder

    new_sample_number = 0

    for r1_file, r2_file in zip(r1, r2):
        # EXtract sample name
        sample = extract_sample(r1_file, r2_file)
        args.sample = sample
        if sample in sample_list_F:

            sample_number = str(sample_list_F.index(sample) + 1)
            sample_total = str(len(sample_list_F))

            out_markdup_trimmed_name = sample + ".rg.markdup.trimmed.sorted.bam"
            output_markdup_trimmed_file = os.path.join(
                out_map_dir, out_markdup_trimmed_name)

            if sample in new_samples:
                new_sample_number = str(int(new_sample_number) + 1)
                new_sample_total = str(len(new_samples))
                logger.info("\n" + WHITE_BG + "STARTING SAMPLE: " + sample +
                            " (" + sample_number + "/" + sample_total + ")" + " (" + new_sample_number + "/" + new_sample_total + ")" + END_FORMATTING)
            else:
                logger.info("\n" + WHITE_BG + "STARTING SAMPLE: " + sample +
                            " (" + sample_number + "/" + sample_total + ")" + END_FORMATTING)

            if not os.path.isfile(output_markdup_trimmed_file):

                args.r1_file = r1_file
                args.r2_file = r2_file

                ##############START PIPELINE#####################
                #################################################

                # INPUT ARGUMENTS
                ################
                check_file_exists(r1_file)
                check_file_exists(r2_file)

                args.output = os.path.abspath(args.output)
                check_create_dir(args.output)

                # QUALITY CHECK in RAW with fastqc
                ######################################################
                check_create_dir(out_qc_dir)

                out_qc_raw_name_r1 = (".").join(r1_file.split(
                    '/')[-1].split('.')[0:-2]) + '_fastqc.html'
                out_qc_raw_name_r2 = (".").join(r2_file.split(
                    '/')[-1].split('.')[0:-2]) + '_fastqc.html'
                output_qc_raw_file_r1 = os.path.join(
                    out_qc_pre_dir, out_qc_raw_name_r1)
                output_qc_raw_file_r2 = os.path.join(
                    out_qc_pre_dir, out_qc_raw_name_r2)

                if os.path.isfile(output_qc_raw_file_r1) and os.path.isfile(output_qc_raw_file_r2):
                    logger.info(YELLOW + DIM + output_qc_raw_file_r1 +
                                " EXIST\nOmmiting QC for sample " + sample + END_FORMATTING)
                else:
                    logger.info(
                        GREEN + "Checking quality in sample " + sample + END_FORMATTING)
                    logger.info("R1: " + r1_file + "\nR2: " + r2_file)
                    fastqc_quality(r1_file, r2_file,
                                   out_qc_pre_dir, args.threads)

                """
                TODO: Human filter
                """

                # QUALITY TRIMMING AND ADAPTER REMOVAL WITH fastp
                ###################################################
                out_trim_name_r1 = sample + ".trimmed_R1.fastq.gz"
                out_trim_name_r2 = sample + ".trimmed_R2.fastq.gz"
                output_trimming_file_r1 = os.path.join(
                    out_trim_dir, out_trim_name_r1)
                output_trimming_file_r2 = os.path.join(
                    out_trim_dir, out_trim_name_r2)

                if os.path.isfile(output_trimming_file_r1) and os.path.isfile(output_trimming_file_r2):
                    logger.info(YELLOW + DIM + output_trimming_file_r1 +
                                " EXIST\nOmmiting Trimming for sample " + sample + END_FORMATTING)
                else:
                    logger.info(GREEN + "Trimming sample " +
                                sample + END_FORMATTING)
                    fastp_trimming(r1_file, r2_file, sample, out_trim_dir,
                                   threads=args.threads, min_qual=20, window_size=10, min_len=35)

                # QUALITY CHECK in TRIMMED with fastqc
                ######################################################
                check_create_dir(out_qc_dir)

                out_qc_pos_r1 = sample + ".trimmed_R1_fastqc.html"
                out_qc_pos_r2 = sample + ".trimmed_R2_fastqc.html"
                output_qc_precessed_file_r1 = os.path.join(
                    out_qc_post_dir, out_qc_pos_r1)
                output_qc_precessed_file_r2 = os.path.join(
                    out_qc_post_dir, out_qc_pos_r2)

                if os.path.isfile(output_qc_precessed_file_r1) and os.path.isfile(output_qc_precessed_file_r2):
                    logger.info(YELLOW + DIM + output_qc_raw_file_r1 +
                                " EXIST\nOmmiting QC for sample " + sample + END_FORMATTING)
                else:
                    logger.info(
                        GREEN + "Checking quality in processed sample " + sample + END_FORMATTING)
                    logger.info("R1: " + output_trimming_file_r1 +
                                "\nR2: " + output_trimming_file_r2)
                    fastqc_quality(
                        output_trimming_file_r1, output_trimming_file_r2, out_qc_post_dir, args.threads)

                # MAPPING WITH BWA - SAM TO SORTED BAM - ADD HEADER SG
                #####################################################
                out_map_name = sample + ".rg.sorted.bam"
                output_map_file = os.path.join(out_map_dir, out_map_name)

                if os.path.isfile(output_map_file):
                    logger.info(YELLOW + DIM + output_map_file +
                                " EXIST\nOmmiting Mapping for sample " + sample + END_FORMATTING)
                else:
                    logger.info(GREEN + "Mapping sample " +
                                sample + END_FORMATTING)
                    logger.info("R1: " + output_trimming_file_r1 + "\nR2: " +
                                output_trimming_file_r2 + "\nReference: " + reference)
                    bwa_mapping(output_trimming_file_r1, output_trimming_file_r2,
                                reference, sample, out_map_dir, threads=args.threads)
                    sam_to_index_bam(
                        sample, out_map_dir, output_trimming_file_r1, threads=args.threads)

                #MARK DUPLICATES WITH PICARDTOOLS ###################
                #####################################################
                out_markdup_name = sample + ".rg.markdup.sorted.bam"
                output_markdup_file = os.path.join(
                    out_map_dir, out_markdup_name)

                if os.path.isfile(output_markdup_file):
                    logger.info(YELLOW + DIM + output_markdup_file +
                                " EXIST\nOmmiting Duplucate Mark for sample " + sample + END_FORMATTING)
                else:
                    logger.info(GREEN + "Marking Dupes in sample " +
                                sample + END_FORMATTING)
                    logger.info("Input Bam: " + output_map_file)
                    picard_markdup(output_map_file)

                #TRIM PRIMERS WITH ivar trim ########################
                #####################################################

                if os.path.isfile(output_markdup_trimmed_file):
                    logger.info(YELLOW + DIM + output_markdup_trimmed_file +
                                " EXIST\nOmmiting Duplucate Mark for sample " + sample + END_FORMATTING)
                else:
                    logger.info(
                        GREEN + "Trimming primers in sample " + sample + END_FORMATTING)
                    logger.info("Input Bam: " + output_markdup_file)
                    ivar_trim(output_markdup_file, args.primers, sample,
                              min_length=30, min_quality=20, sliding_window_width=4)
            else:
                logger.info(YELLOW + DIM + output_markdup_trimmed_file +
                            " EXIST\nOmmiting BAM mapping and BAM manipulation in sample " + sample + END_FORMATTING)

            ########################END OF MAPPING AND BAM MANIPULATION#####################################################################
            ################################################################################################################################

            #VARIANT CALLING WTIH ivar variants##################
            #####################################################
            check_create_dir(out_variant_dir)
            out_ivar_variant_name = sample + ".tsv"
            out_ivar_variant_file = os.path.join(
                out_variant_ivar_dir, out_ivar_variant_name)

            if os.path.isfile(out_ivar_variant_file):
                logger.info(YELLOW + DIM + out_ivar_variant_file +
                            " EXIST\nOmmiting Variant call for  sample " + sample + END_FORMATTING)
            else:
                logger.info(
                    GREEN + "Calling variants with ivar in sample " + sample + END_FORMATTING)
                ivar_variants(reference, output_markdup_trimmed_file, out_variant_dir, sample,
                              annotation, min_quality=15, min_frequency_threshold=0.01, min_depth=1)

            #VARIANT FILTERING ##################################
            #####################################################
            check_create_dir(out_filtered_ivar_dir)
            out_ivar_filtered_file = os.path.join(
                out_filtered_ivar_dir, out_ivar_variant_name)

            if os.path.isfile(out_ivar_filtered_file):
                logger.info(YELLOW + DIM + out_ivar_filtered_file +
                            " EXIST\nOmmiting Variant filtering for  sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + "Filtering variants in sample " +
                            sample + END_FORMATTING)
                filter_tsv_variants(out_ivar_variant_file, out_filtered_ivar_dir, min_frequency=0.7,
                                    min_total_depth=10, min_alt_dp=4, is_pass=True, only_snp=False)

            #CREATE CONSENSUS with ivar consensus##################
            #######################################################
            check_create_dir(out_consensus_dir)
            check_create_dir(out_consensus_ivar_dir)
            out_ivar_consensus_name = sample + ".fa"
            out_ivar_consensus_file = os.path.join(
                out_consensus_ivar_dir, out_ivar_consensus_name)

            if os.path.isfile(out_ivar_consensus_file):
                logger.info(YELLOW + DIM + out_ivar_consensus_file +
                            " EXIST\nOmmiting Consensus for  sample " + sample + END_FORMATTING)
            else:
                logger.info(
                    GREEN + "Creating consensus with ivar in sample " + sample + END_FORMATTING)
                ivar_consensus(output_markdup_trimmed_file, out_consensus_ivar_dir, sample,
                               min_quality=20, min_frequency_threshold=0.8, min_depth=20, uncovered_character='N')
                logger.info(
                    GREEN + "Replacing consensus header in " + sample + END_FORMATTING)
                replace_consensus_header(out_ivar_consensus_file)

            ########################CREATE STATS AND QUALITY FILTERS########################################################################
            ################################################################################################################################
            #CREATE Bamstats#######################################
            #######################################################
            check_create_dir(out_stats_dir)
            check_create_dir(out_stats_bamstats_dir)
            out_bamstats_name = sample + ".bamstats"
            out_bamstats_file = os.path.join(
                out_stats_bamstats_dir, out_bamstats_name)

            if os.path.isfile(out_bamstats_file):
                logger.info(YELLOW + DIM + out_bamstats_file +
                            " EXIST\nOmmiting Bamstats for  sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + "Creating bamstats in sample " +
                            sample + END_FORMATTING)
                create_bamstat(output_markdup_trimmed_file,
                               out_stats_bamstats_dir, sample, threads=args.threads)

            #CREATE Bamstats#######################################
            #######################################################
            check_create_dir(out_stats_coverage_dir)
            out_coverage_name = sample + ".cov"
            out_coverage_file = os.path.join(
                out_stats_coverage_dir, out_coverage_name)

            if os.path.isfile(out_coverage_file):
                logger.info(YELLOW + DIM + out_coverage_file +
                            " EXIST\nOmmiting Bamstats for  sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + "Creating coverage in sample " +
                            sample + END_FORMATTING)
                create_coverage(output_markdup_trimmed_file,
                                out_stats_coverage_dir, sample)

    # fastqc OUTPUT FORMAT FOR COMPARISON
    ######################################################
    logger.info(
        GREEN + "Creating summary report for quality result " + END_FORMATTING)
    # format_html_image(out_qc_dir)

    # coverage OUTPUT SUMMARY
    ######################################################
    logger.info(
        GREEN + "Creating summary report for coverage result " + END_FORMATTING)
    obtain_group_cov_stats(out_stats_coverage_dir, group_name)

    # READS and VARIANTS OUTPUT SUMMARY
    ######################################################
    logger.info(GREEN + "Creating overal summary report " + END_FORMATTING)
    obtain_overal_stats(output, group_name)

    # REMOVE UNCOVERED
    ##############################################################################################################################
    logger.info(GREEN + "Removing low quality samples" + END_FORMATTING)
    remove_low_quality(output, min_coverage=args.min_coverage, max_unmap=args.max_unmap)

    #ANNOTATION WITH SNPEFF, USER INOUT AND PANGOLIN ####
    #####################################################
    logger.info("\n\n" + BLUE + BOLD + "STARTING ANNOTATION IN GROUP: " +
                group_name + END_FORMATTING + "\n")
    check_create_dir(out_annot_dir)
    check_create_dir(out_annot_snpeff_dir)
    check_create_dir(out_annot_pangolin_dir)
    # SNPEFF
    if args.snpeff_database != False:
        # CHANGE FOR RAW/FILTERED ANNOTATION
        for root, _, files in os.walk(out_filtered_ivar_dir):
            if root == out_filtered_ivar_dir:  # CHANGE FOR RAW/FILTERED ANNOTATION
                for name in files:
                    if name.endswith('.tsv'):
                        sample = name.split('.')[0]
                        filename = os.path.join(root, name)
                        out_annot_file = os.path.join(
                            out_annot_snpeff_dir, sample + ".annot")
                        if os.path.isfile(out_annot_file):
                            logger.info(YELLOW + DIM + out_annot_file +
                                        " EXIST\nOmmiting snpEff Annotation for sample " + sample + END_FORMATTING)
                        else:
                            logger.info(
                                GREEN + "Annotating sample with snpEff: " + sample + END_FORMATTING)
                            output_vcf = os.path.join(
                                out_annot_snpeff_dir, sample + '.vcf')
                            annotate_snpeff(
                                filename, output_vcf, out_annot_file, database=args.snpeff_database)
    # USER DEFINED
    if not args.annot_bed and not args.annot_vcf:
        logger.info(
            YELLOW + BOLD + "Ommiting User Annotation, no BED or VCF files supplied" + END_FORMATTING)
    else:
        check_create_dir(out_annot_user_dir)
        # CHANGE FOR RAW/FILTERED ANNOTATION
        for root, _, files in os.walk(out_variant_ivar_dir):
            if root == out_variant_ivar_dir:  # CHANGE FOR RAW/FILTERED ANNOTATION
                for name in files:
                    if name.endswith('.tsv'):
                        sample = name.split('.')[0]
                        logger.info(
                            'User bed/vcf annotation in sample {}'.format(sample))
                        filename = os.path.join(root, name)
                        out_annot_file = os.path.join(
                            out_annot_user_dir, sample + ".tsv")
                        user_annotation(
                            filename, out_annot_file, vcf_files=args.annot_vcf, bed_files=args.annot_bed)

    # USER AA DEFINED
    if not args.annot_aa:
        logger.info(
            YELLOW + BOLD + "Ommiting User aa Annotation, no AA files supplied" + END_FORMATTING)
    else:
        check_create_dir(out_annot_user_aa_dir)
        for root, _, files in os.walk(out_annot_snpeff_dir):
            if root == out_annot_snpeff_dir:
                for name in files:
                    if name.endswith('.annot'):
                        sample = name.split('.')[0]
                        logger.info(
                            'User aa annotation in sample {}'.format(sample))
                        filename = os.path.join(root, name)
                        out_annot_aa_file = os.path.join(
                            out_annot_user_aa_dir, sample + ".tsv")
                        if os.path.isfile(out_annot_aa_file):
                            user_annotation_aa(
                                out_annot_aa_file, out_annot_aa_file, aa_files=args.annot_aa)
                        else:
                            user_annotation_aa(
                                filename, out_annot_aa_file, aa_files=args.annot_aa)

    # PANGOLIN
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures_pangolin = []

        for root, _, files in os.walk(out_consensus_ivar_dir):
            if root == out_consensus_ivar_dir:
                for name in files:
                    if name.endswith('.fa'):
                        sample = name.split('.')[0]
                        filename = os.path.join(root, name)
                        out_pangolin_filename = sample + ".lineage.csv"
                        out_pangolin_file = os.path.join(
                            out_annot_pangolin_dir, out_pangolin_filename)
                        if os.path.isfile(out_pangolin_file):
                            logger.info(YELLOW + DIM + out_pangolin_file +
                                        " EXIST\nOmmiting Lineage for  sample " + sample + END_FORMATTING)
                        else:
                            logger.info(
                                GREEN + "Obtaining Lineage in sample " + sample + END_FORMATTING)
                            future = executor.submit(
                                annotate_pangolin, filename, out_annot_pangolin_dir, out_pangolin_filename, threads=args.threads, max_ambig=0.6)
                            futures_pangolin.append(future)
                for future in concurrent.futures.as_completed(futures_pangolin):
                    logger.info(future.result())
                    # annotate_pangolin(filename, out_annot_pangolin_dir,
                    #                out_pangolin_filename, threads=args.threads, max_ambig=0.6)

    # USER AA TO HTML
    annotated_samples = []
    logger.info('Adapting annotation to html in {}'.format(group_name))
    for root, _, files in os.walk(out_annot_user_aa_dir):
        if root == out_annot_user_aa_dir:
            for name in files:
                if name.endswith('.tsv'):
                    sample = name.split('.')[0]
                    annotated_samples.append(sample)
                    filename = os.path.join(root, name)
                    annotation_to_html(filename, sample)
    annotated_samples = [str(x) for x in annotated_samples]
    report_samples_html_all = report_samples_html.replace(
        'ALLSAMPLES', ('","').join(annotated_samples))  # NEW
    with open(os.path.join(out_annot_user_aa_dir, '00_all_samples.html'), 'w+') as f:
        f.write(report_samples_html_all)

    # SNP COMPARISON using tsv variant files
    ######################################################
    logger.info("\n\n" + BLUE + BOLD + "STARTING COMPARISON IN GROUP: " +
                group_name + END_FORMATTING + "\n")

    check_create_dir(out_compare_dir)
    folder_compare = today + "_" + group_name
    path_compare = os.path.join(out_compare_dir, folder_compare)
    check_create_dir(path_compare)
    full_path_compare = os.path.join(path_compare, group_name)

    # ddtb_add(out_filtered_ivar_dir, full_path_compare)
    compare_snp_matrix_recal = full_path_compare + ".revised.final.tsv"
    compare_snp_matrix_INDEL = full_path_compare + ".revised_INDEL.final.tsv"
    compare_snp_matrix_recal_intermediate = full_path_compare + ".revised_intermediate.tsv"
    compare_snp_matrix_INDEL_intermediate = full_path_compare + \
        ".revised_INDEL_intermediate.tsv"
    recalibrated_snp_matrix_intermediate = ddbb_create_intermediate(
        out_variant_ivar_dir, out_stats_coverage_dir, min_freq_discard=0.1, min_alt_dp=4, only_snp=args.only_snp)
    recalibrated_snp_matrix_intermediate.to_csv(
        compare_snp_matrix_recal_intermediate, sep="\t", index=False)
    compare_snp_matrix_INDEL_intermediate_df = remove_position_range(
        recalibrated_snp_matrix_intermediate)
    compare_snp_matrix_INDEL_intermediate_df.to_csv(
        compare_snp_matrix_INDEL_intermediate, sep="\t", index=False)
    recalibrated_revised_df = revised_df(recalibrated_snp_matrix_intermediate, path_compare, min_freq_include=0.7,
                                         min_threshold_discard_sample=0.07, min_threshold_discard_position=0.4, remove_faulty=True, drop_samples=True, drop_positions=True)
    recalibrated_revised_df.to_csv(
        compare_snp_matrix_recal, sep="\t", index=False)
    recalibrated_revised_INDEL_df = revised_df(compare_snp_matrix_INDEL_intermediate_df, path_compare, min_freq_include=0.7,
                                               min_threshold_discard_sample=0.07, min_threshold_discard_position=0.4, remove_faulty=True, drop_samples=True, drop_positions=True)
    recalibrated_revised_INDEL_df.to_csv(
        compare_snp_matrix_INDEL, sep="\t", index=False)

    ddtb_compare(compare_snp_matrix_recal, distance=0)
    ddtb_compare(compare_snp_matrix_INDEL, distance=0, indel=True)

    logger.info("\n\n" + MAGENTA + BOLD + "COMPARING FINISHED IN GROUP: " +
                group_name + END_FORMATTING + "\n")

    #####################CONSENSUS WITH REFINED CALL######
    ######################################################
    logger.info(GREEN + "Creating refined consensus" + END_FORMATTING)
    create_consensus(reference, compare_snp_matrix_recal,
                     out_stats_coverage_dir, out_consensus_dir)

    logger.info("\n\n" + MAGENTA + BOLD +
                "#####END OF PIPELINE COVID MULTI ANALYSIS#####" + END_FORMATTING + "\n")


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        logger.exception(e)
        raise
