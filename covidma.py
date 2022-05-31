#!/usr/bin/env python

# Standard library imports
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
from compare_snp import ddtb_add, ddtb_compare, ddbb_create_intermediate, ddbb_create_intermediate_ori, revised_df, remove_position_range

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

### ARGUMENTS

def get_arguments():

    parser = argparse.ArgumentParser(
        prog='covidma.py', description='Pipeline to call variants (SNVs) with any non model organism. Specialised in SARS-CoV-2')

    input_group = parser.add_argument_group('Input', 'Input parameters')

    input_group.add_argument('-i', '--input', dest="input_dir", metavar="input_directory",
                                type=str, required=True, help='REQUIRED.Input directory containing all fast[aq] files')
    input_group.add_argument('-name', '--name', dest="name_sbatch", metavar="name sbatch",
                                type=str, required=True, help='REQUIRED.Name for sbatch. It must contain 3 letters, if not, the first three letters are considered')
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

    quality_group.add_argument('-c', '--coverage20', type=int, default=90, required=False,
                                help='Minimum percentage of coverage at 20x to clasify as uncovered (Default 90)')
    quality_group.add_argument('-n', '--min_snp', type=int, required=False,
                                default=1, help='SNP number to pass quality threshold')

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

def create_logFile(group_name, logger, output, args):
    """
    Function that creates a log file to store all steps 
    of the pipeline.

    In output/Logs a .log file is created
    """

    # Name of the filename
    log_filename = group_name + "_" + str(datetime.datetime.now().date()) + "_" + str(datetime.datetime.now().time()) + ".log"
    # Folder Logs
    log_folder = os.path.join(output, 'Logs')
    # Create Logs file
    check_create_dir(log_folder)
    # Absolute path
    log_full_path = os.path.join(log_folder, log_filename)

    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s:%(message)s')

    file_handler = logging.FileHandler(log_full_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)
    logger.info("\n\n" + "STARTING PIPELINE IN GROUP: " + group_name + END_FORMATTING)
    logger.info("ARGUMENTS:")
    logger.info(str(args))

def obtain_reads(args, logger):
    """
    Function that extract name reads files and check wether
    samples need reanalysis or not
    """

    # Obtain all R1 and R2 files from folder
    r1, r2 = extract_read_list(args.input_dir)

    # Check if there are samples to filter out
    sample_list_F = []
    if args.sample_list == None:
        logger.info("\n" + "No samples to filter")
        # As not samples need to be filtered, we extract all samples
        for r1_file, r2_file in zip(r1, r2):
            sample = extract_sample(r1_file, r2_file)
            sample_list_F.append(sample)
    else:
        logger.info("samples will be filtered")
        sample_list_F = file_to_list(args.sample_list)

    # Check if samples have been already analysed (Bam directory is created)
    # new samples contain samples not already analysed second time
    # a run is performed
    new_samples = check_reanalysis(args.output, sample_list_F)
    logger.info("\n%d samples will be analysed: %s" %
                (len(new_samples), ",".join(new_samples)))
    
    # Prepare reference for mapping + fai + dict
    samtools_faidx(args)

    # Returns paired-reads absolute path and samples ID
    return (r1, r2, sample_list_F, new_samples)

def set_folder_structure(output):

    # output/Quality
    out_qc_dir = os.path.join(output, "Quality")
    check_create_dir(out_qc_dir)
    # output/Quality/raw
    out_qc_pre_dir = os.path.join(out_qc_dir, "raw")  # subfolder
    check_create_dir(out_qc_pre_dir)
    # output/Quality/processed
    out_qc_post_dir = os.path.join(out_qc_dir, "processed")  # subfolder
    check_create_dir(out_qc_post_dir)

    # output/Trimmed
    out_trim_dir = os.path.join(output, "Trimmed")
    check_create_dir(out_trim_dir)

    # output/Bam
    out_map_dir = os.path.join(output, "Bam")
    check_create_dir(out_map_dir)

    # output/Variants
    out_variant_dir = os.path.join(output, "Variants")
    check_create_dir(out_variant_dir)
    # output/Variants/ivar_raw
    out_variant_ivar_dir = os.path.join(
        out_variant_dir, "ivar_raw")  # subfolder
    check_create_dir(out_variant_ivar_dir)
    # output/Variants/ivar_filtered
    out_filtered_ivar_dir = os.path.join(
        out_variant_dir, "ivar_filtered")  # subfolder
    check_create_dir(out_filtered_ivar_dir)

    # output/Consensus
    out_consensus_dir = os.path.join(output, "Consensus")
    check_create_dir(out_consensus_dir)
    # output/Consensus/ivar
    out_consensus_ivar_dir = os.path.join(
        out_consensus_dir, "ivar")  # subfolder
    check_create_dir(out_consensus_ivar_dir)

    # output/Stats
    out_stats_dir = os.path.join(output, "Stats")
    check_create_dir(out_stats_dir)
    # output/Stats/Bamstats
    out_stats_bamstats_dir = os.path.join(
        out_stats_dir, "Bamstats")  # subfolder
    check_create_dir(out_stats_bamstats_dir)
    # output/Stats/Coverage
    out_stats_coverage_dir = os.path.join(
        out_stats_dir, "Coverage")  # subfolder
    check_create_dir(out_stats_coverage_dir)

    # output/Compare
    out_compare_dir = os.path.join(output, "Compare")
    check_create_dir(out_compare_dir)

    # output/Annotation
    out_annot_dir = os.path.join(output, "Annotation")
    check_create_dir(out_annot_dir)
    # output/Annotation/snpeff
    out_annot_snpeff_dir = os.path.join(out_annot_dir, "snpeff")  # subfolder
    check_create_dir(out_annot_snpeff_dir)
    # output/Annotation/pangolin
    out_annot_pangolin_dir = os.path.join(
        out_annot_dir, "pangolin")  # subfolder
    check_create_dir(out_annot_pangolin_dir)
    # output/Annotation/user
    out_annot_user_dir = os.path.join(out_annot_dir, "user")  # subfolder
    check_create_dir(out_annot_user_dir)
    # output/Annotation/user_aa
    out_annot_user_aa_dir = os.path.join(out_annot_dir, "user_aa")  # subfolder
    check_create_dir(out_annot_user_aa_dir)

def user_aa_to_html(out_annot_user_aa_dir, group_name):
    """
    Function that converts user_aa annotation to html.

    Html files are generated in Annotation/user_aa.
    """

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

def snp_comparison(name_s, logger, output, group_name, out_variant_ivar_dir, out_stats_coverage_dir):
    """
    Function that performs SNP comparison. As output is obtained a 
    pairwise distance matrix and dendrogram. 

    Output is located in Compare/FOLDER/.
    """

    logger.info("\n\n" + BLUE + BOLD + "STARTING COMPARISON IN GROUP: " +
                group_name + END_FORMATTING + "\n")
    out_consensus_dir = os.path.join(output, "Consensus")                                       # Folder
    

    out_compare_dir = os.path.join(output, "Compare")       # Folder
    folder_compare = str(datetime.date.today()) + "_" + group_name
    path_compare = os.path.join(out_compare_dir, folder_compare)
    check_create_dir(path_compare)
    full_path_compare = os.path.join(path_compare, group_name)

    # ddtb_add(out_filtered_ivar_dir, full_path_compare)
    compare_snp_matrix_recal = full_path_compare + ".revised.final.tsv"
    compare_snp_matrix_INDEL = full_path_compare + ".revised_INDEL.final.tsv"
    compare_snp_matrix_recal_intermediate = full_path_compare + ".revised_intermediate.tsv"
    compare_snp_matrix_INDEL_intermediate = full_path_compare + ".revised_INDEL_intermediate.tsv"

    # Get number proc
    n_files = len([f for f in os.listdir(out_variant_ivar_dir) if f.endswith(".tsv")])
    if n_files > 6000:
        nproc = 96
    elif n_files > 4000:
        nproc = 64
    elif n_files > 2000:
        nproc = 32
    elif n_files > 500:
        nproc = 16
    elif n_files > 100:
        nproc = 8
    elif n_files > 50:
        nproc = 4
    elif n_files > 30:
        nproc = 2
    else:
        nproc = 1

    if nproc > 1:
        logger.info("\n\n" + "USING: " +
                str(nproc) + " cores" + "\n")
        recalibrated_snp_matrix_intermediate = ddbb_create_intermediate(name_s, 
                path_compare, out_variant_ivar_dir, out_stats_coverage_dir, min_freq_discard=0.1, 
                min_alt_dp=4, only_snp=False, nproc=nproc)
    else:
        logger.info("\n\n" + "USING: " +
                str(nproc) + " core" + "\n")
        recalibrated_snp_matrix_intermediate = ddbb_create_intermediate_ori(
        out_variant_ivar_dir, out_stats_coverage_dir, min_freq_discard=0.1, min_alt_dp=4, only_snp=args.only_snp)
    
    recalibrated_snp_matrix_intermediate.to_csv(compare_snp_matrix_recal_intermediate, sep="\t", index=False)

    compare_snp_matrix_INDEL_intermediate_df = remove_position_range(recalibrated_snp_matrix_intermediate)
    compare_snp_matrix_INDEL_intermediate_df.to_csv(compare_snp_matrix_INDEL_intermediate, sep="\t", index=False)

    recalibrated_revised_df = revised_df(recalibrated_snp_matrix_intermediate, path_compare, min_freq_include=0.7,
                                         min_threshold_discard_sample=0.07, min_threshold_discard_position=0.4, remove_faulty=True, drop_samples=True, drop_positions=True)
    recalibrated_revised_df.to_csv(compare_snp_matrix_recal, sep="\t", index=False)
    
    recalibrated_revised_INDEL_df = revised_df(compare_snp_matrix_INDEL_intermediate_df, path_compare, min_freq_include=0.7,
                                               min_threshold_discard_sample=0.07, min_threshold_discard_position=0.4, remove_faulty=True, drop_samples=True, drop_positions=True)
    recalibrated_revised_INDEL_df.to_csv(compare_snp_matrix_INDEL, sep="\t", index=False)

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
        

def covidma(output, args, logger, r1, r2, sample_list_F, new_samples, group_name, reference, annotation):

    # Loop for paralellization
    name_s = args.name_sbatch[:3]      # name sbatch
    l = len(r1) - 1

    # Loop to fix posible problems
    for time in range(3):
        counter = 0
        new_sample_number = 0

        for i in range(len(r1)):
            r1_file = r1[i]
            r2_file = r2[i]
            
            # Extract sample name
            sample = extract_sample(r1_file, r2_file)
            args.sample = sample
            if sample in sample_list_F:
                sample_number = str(sample_list_F.index(sample) + 1)
                sample_total = str(len(sample_list_F))

                out_map_dir = os.path.join(output, "Bam")
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
                
                # Necessary variables
                out_variant_dir = os.path.join(output, "Variants") 
                out_filtered_ivar_dir = os.path.join(out_variant_dir, "ivar_filtered") 
                out_variant_ivar_dir = os.path.join(out_variant_dir, "ivar_raw")
                out_stats_dir = os.path.join(output, "Stats")  
                out_stats_coverage_dir = os.path.join(out_stats_dir, "Coverage")
                out_stats_bamstats_dir = os.path.join(out_stats_dir, "Bamstats")  # subfolder
                out_annot_dir = os.path.join(output, "Annotation")              # Folder
                out_annot_user_aa_dir = os.path.join(out_annot_dir, "user_aa")  # subfolder
                out_consensus_dir = os.path.join(output, "Consensus")               # Folder
                out_consensus_ivar_dir = os.path.join(out_consensus_dir, "ivar")    # subfolder
                out_annot_pangolin_dir = os.path.join(out_annot_dir, "pangolin")    # subfolder
                threads = str(args.threads)
                primers = args.primers

                # create annot_vcf file
                f_annot_vcf = open(output + "/" + "annot_vcf.txt", "w")
                for s in args.annot_vcf:
                    to_write = s.strip() + "\n"
                    f_annot_vcf.write(to_write)
                f_annot_vcf.close()

                # create annot_aa file
                f_annot_aa = open(output + "/" + "annot_aa.txt", "w")
                for s in args.annot_aa:
                    to_write = s.strip() + "\n"
                    f_annot_aa.write(to_write)
                f_annot_aa.close()

                # create annot_bed file
                f_annot_bed = open(output + "/" + "annot_bed.txt", "w")
                for s in args.annot_bed:
                    to_write = s.strip() + "\n"
                    f_annot_bed.write(to_write)
                f_annot_bed.close()

                ######################################################
                # QUALITY CHECK in RAW with fastqc
                # QUALITY TRIMMING AND ADAPTER REMOVAL WITH fastp
                # QUALITY CHECK in TRIMMED with fastqc
                # MAPPING WITH BWA - SAM TO SORTED BAM - ADD HEADER SG
                #MARK DUPLICATES WITH PICARDTOOLS 
                #TRIM PRIMERS WITH ivar trim
                ######################################################
                if os.path.isfile(output_markdup_trimmed_file) and os.stat(output_markdup_trimmed_file).st_size:
                    logger.info(YELLOW + DIM + output_markdup_trimmed_file +
                                " EXIST\nOmmiting BAM mapping and BAM manipulation in sample " + sample + END_FORMATTING) 
                else:
                    os.system("sbatch -J %s /home/laura/Laura_intel/Desktop/covid_multianalysis/sample_mapper.sh %s %s %s %s %s %s %s > jobid.batch" 
                            %(name_s + "1", output, r1_file, r2_file, sample, reference, threads, primers))


                #VARIANT CALLING WTIH ivar variants##################
                #####################################################
                out_ivar_variant_name = sample + ".tsv"
                out_ivar_variant_file = os.path.join(
                    out_variant_ivar_dir, out_ivar_variant_name)
                if os.path.isfile(out_ivar_variant_file) and os.stat(out_ivar_variant_file).st_size > 140:
                    logger.info(YELLOW + DIM + out_ivar_variant_file +
                                " EXIST\nOmmiting Variant call for  sample " + sample + END_FORMATTING)
                else:
                    if os.path.isfile(out_ivar_variant_file) and os.stat(out_ivar_variant_file).st_size:
                        os.remove(out_ivar_variant_file)
                    logger.info(
                        GREEN + "Calling variants with ivar in sample " + sample + END_FORMATTING)
                    if os.path.exists("jobid.batch"):
                        f = open("jobid.batch", "r")
                        jobid = f.readline().strip().split()[-1]
                        f.close()
                        os.remove("jobid.batch")
                        os.system("sbatch -J %s --dependency=afterok:%s /home/laura/Laura_intel/Desktop/covid_multianalysis/ivar_variant_calling.sh %s %s %s %s %s > jobid.batch" 
                            %(name_s + "2", jobid, output, output_markdup_trimmed_file, sample, reference, annotation))
                    else:
                        os.system("sbatch -J %s /home/laura/Laura_intel/Desktop/covid_multianalysis/ivar_variant_calling.sh %s %s %s %s %s > jobid.batch" 
                            %(name_s + "2", output, output_markdup_trimmed_file, sample, reference, annotation))


                #VARIANT FILTERING ##################################
                #####################################################
                out_ivar_filtered_file = os.path.join(
                    out_filtered_ivar_dir, out_ivar_variant_name)
                if os.path.isfile(out_ivar_filtered_file) and os.stat(out_ivar_filtered_file).st_size > 140:
                    logger.info(YELLOW + DIM + out_ivar_filtered_file +
                                " EXIST\nOmmiting Variant filtering for  sample " + sample + END_FORMATTING)
                else:
                    if os.path.isfile(out_ivar_filtered_file) and os.stat(out_ivar_filtered_file).st_size:
                        os.remove(out_ivar_filtered_file)
                    logger.info(GREEN + "Filtering variants in sample " +
                                sample + END_FORMATTING)
                    if os.path.exists("jobid.batch"):
                        f = open("jobid.batch", "r")
                        jobid = f.readline().strip().split()[-1]
                        f.close()
                        os.remove("jobid.batch")
                        os.system("sbatch -J %s --dependency=afterok:%s /home/laura/Laura_intel/Desktop/covid_multianalysis/variant_filtering.sh %s %s %s > jobid.batch" 
                            %(name_s + "3", jobid, output, sample, out_ivar_variant_file))
                    else:
                        os.system("sbatch -J %s /home/laura/Laura_intel/Desktop/covid_multianalysis/variant_filtering.sh %s %s %s > jobid.batch" 
                            %(name_s + "3", output, sample, out_ivar_variant_file))
                
                #CREATE CONSENSUS with ivar consensus##################
                #######################################################
                out_ivar_consensus_name = sample + ".fa"
                out_ivar_consensus_file = os.path.join(
                    out_consensus_ivar_dir, out_ivar_consensus_name)

                if os.path.isfile(out_ivar_consensus_file) and os.stat(out_ivar_consensus_file).st_size > 20000:
                    logger.info(YELLOW + DIM + out_ivar_consensus_file +
                                " EXIST\nOmmiting Consensus for  sample " + sample + END_FORMATTING)
                else:
                    if os.path.isfile(out_ivar_consensus_file) and os.stat(out_ivar_consensus_file).st_size:
                        os.remove(out_ivar_consensus_file)
                    logger.info(
                        GREEN + "Creating consensus with ivar in sample " + sample + END_FORMATTING)
                    if os.path.exists("jobid.batch"):
                        f = open("jobid.batch", "r")
                        jobid = f.readline().strip().split()[-1]
                        f.close()
                        os.remove("jobid.batch")
                        os.system("sbatch -J %s --dependency=afterok:%s /home/laura/Laura_intel/Desktop/covid_multianalysis/consensus_ivar_create.sh %s %s %s > jobid.batch" 
                            %(name_s + "4", jobid, output, sample, output_markdup_trimmed_file))
                    else:
                        os.system("sbatch -J %s /home/laura/Laura_intel/Desktop/covid_multianalysis/consensus_ivar_create.sh %s %s %s > jobid.batch" 
                            %(name_s + "4", output, sample, output_markdup_trimmed_file))
                
                ########################CREATE STATS AND QUALITY FILTERS###############################
                #######################################################################################
                #CREATE Bamstats #######################################
                ########################################################
                out_bamstats_name = sample + ".bamstats"
                out_bamstats_file = os.path.join(
                    out_stats_bamstats_dir, out_bamstats_name)

                if os.path.isfile(out_bamstats_file) and os.stat(out_bamstats_file).st_size:
                    logger.info(YELLOW + DIM + out_bamstats_file +
                                " EXIST\nOmmiting Bamstats for  sample " + sample + END_FORMATTING)
                else:
                    if os.path.isfile(out_bamstats_file) and not os.stat(out_bamstats_file).st_size:
                        os.remove(out_bamstats_file)
                    logger.info(GREEN + "Creating bamstats in sample " +
                                sample + END_FORMATTING)
                    if os.path.exists("jobid.batch"):
                        f = open("jobid.batch", "r")
                        jobid = f.readline().strip().split()[-1]
                        f.close()
                        os.remove("jobid.batch")
                        os.system("sbatch -J %s --dependency=afterok:%s /home/laura/Laura_intel/Desktop/covid_multianalysis/bamstats.sh %s %s %s %s > jobid.batch" 
                            %(name_s + "5", jobid, output, sample, output_markdup_trimmed_file, threads))
                    else:
                        os.system("sbatch -J %s /home/laura/Laura_intel/Desktop/covid_multianalysis/bamstats.sh %s %s %s %s > jobid.batch" 
                            %(name_s + "5", output, sample, output_markdup_trimmed_file, threads))

                #CREATE CoverageStats ##################################
                ########################################################
                out_coverage_name = sample + ".cov"
                out_coverage_file = os.path.join(
                    out_stats_coverage_dir, out_coverage_name)

                if os.path.isfile(out_coverage_file) and os.stat(out_coverage_file).st_size:
                    logger.info(YELLOW + DIM + out_coverage_file +
                                " EXIST\nOmmiting Coverage for sample " + sample + END_FORMATTING)
                else:
                    if os.path.isfile(out_coverage_file) and not os.stat(out_coverage_file).st_size:
                        os.remove(out_coverage_file)
                    logger.info(GREEN + "Creating coverage in sample " +
                                sample + END_FORMATTING)
                    if os.path.exists("jobid.batch"):
                        f = open("jobid.batch", "r")
                        jobid = f.readline().strip().split()[-1]
                        f.close()
                        os.remove("jobid.batch")
                        os.system("sbatch -J %s --dependency=afterok:%s /home/laura/Laura_intel/Desktop/covid_multianalysis/coverage_stats.sh %s %s %s > jobid.batch" 
                            %(name_s + "6", jobid, output, sample, output_markdup_trimmed_file))
                    else:
                        os.system("sbatch -J %s /home/laura/Laura_intel/Desktop/covid_multianalysis/coverage_stats.sh %s %s %s > jobid.batch" 
                            %(name_s + "6", output, output, sample, output_markdup_trimmed_file))
                            
                if os.path.exists("jobid.batch"):
                    os.remove("jobid.batch")
                
                # CHECK PARALELL
                os.system('while [ "$(squeue | grep $USER | grep "%s" | wc -l)" -ge "60" ]; do sleep 0.1; done' %(name_s))
                os.system('if [ %s = %s ]; then while [ $(squeue | grep $USER | grep "%s" | wc -l) != 0 ]; do sleep 0.1; done; fi' %(str(counter), l, name_s))
                counter += 1
        os.system('while [ $(squeue | grep $USER | grep "%s" | wc -l) != 0 ]; do sleep 0.1; done' %(name_s))

    # coverage OUTPUT SUMMARY
    ######################################################
    logger.info(GREEN + "Creating summary report for coverage result " + END_FORMATTING)
    # Output is located in output/Stats/Coverage
    obtain_group_cov_stats(out_stats_coverage_dir, group_name)

    # READS and VARIANTS OUTPUT SUMMARY
    ######################################################
    logger.info(GREEN + "Creating overal summary report " + END_FORMATTING)
    # Output is located in output/Stats
    obtain_overal_stats(output, group_name)

    #ANNOTATION WITH SNPEFF #############################
    #####################################################
    logger.info("\n\n" + BLUE + BOLD + "STARTING ANNOTATION IN GROUP: " +
        group_name + END_FORMATTING + "\n")
    
    # Annotates variants using snpEff. Output is located in output/Annotation/snpeff.
    out_annot_dir = os.path.join(output, "Annotation")              # Folder
    out_annot_snpeff_dir = os.path.join(out_annot_dir, "snpeff")    # subfolder

    if args.snpeff_database:
        snpeff_database = args.snpeff_database
        # CHANGE FOR RAW/FILTERED ANNOTATION
        counter = 0
        l = len([file for file in os.listdir(out_filtered_ivar_dir) if file.endswith(".tsv")]) - 1
        for root, _, files in os.walk(out_filtered_ivar_dir):
            if root == out_filtered_ivar_dir:  # CHANGE FOR RAW/FILTERED ANNOTATION
                for name in files:
                    if name.endswith('.tsv'):
                        logger.info(GREEN + "Creating SNPEFF annotation " + name + END_FORMATTING)
                        os.system("sbatch -J %s /home/laura/Laura_intel/Desktop/covid_multianalysis/snpeff_annotation.sh %s %s %s %s %s > jobid.batch" 
                        %(name_s + "snf", output, name, root, sample, snpeff_database))
                        os.system('while [ "$(squeue | grep $USER | grep "%s" | wc -l)" = "96" ]; do sleep 0.1; done' %(name_s))
                        os.system('if [ %s = %s ]; then while [ $(squeue | grep $USER | grep "%s" | wc -l) != 0 ]; do sleep 0.1; done; fi' %(str(counter), l, name_s))
                        counter += 1

    # USER DEFINED ANNOTATION ###########################
    #####################################################
    # Annotate variants using user files. Output is located in output/Annotation/user.
    if not args.annot_bed and not args.annot_vcf:
        logger.info(
            YELLOW + BOLD + "Ommiting User Annotation, no BED or VCF files supplied" + END_FORMATTING)
    else:
        # CHANGE FOR RAW/FILTERED ANNOTATION
        counter = 0
        l = len([file for file in os.listdir(out_variant_ivar_dir) if file.endswith(".tsv")]) - 1
        for root, _, files in os.walk(out_variant_ivar_dir):
            if root == out_variant_ivar_dir:  # CHANGE FOR RAW/FILTERED ANNOTATION
                for name in files:
                    if name.endswith('.tsv'):
                        logger.info(GREEN + "Creating User annotation " + name + END_FORMATTING)
                        os.system("sbatch -J %s /home/laura/Laura_intel/Desktop/covid_multianalysis/annotate_user.sh %s %s %s %s %s > jobid.batch" 
                        %(name_s + "use", root, output, sample, output + "/" + "annot_vcf.txt", output + "/" + "annot_bed.txt"))
                        os.system('while [ "$(squeue | grep $USER | grep "%s" | wc -l)" = "96" ]; do sleep 0.1; done' %(name_s))
                        os.system('if [ %s = %s ]; then while [ $(squeue | grep $USER | grep "%s" | wc -l) != 0 ]; do sleep 0.1; done; fi' %(str(counter), l, name_s))
                        counter += 1

    # USER AA DEFINED ANNOTATION ########################
    #####################################################
    # Annotate variants using user aa files. Output is located in output/Annotation/user_aa.
    # USER AA DEFINED
    if not args.annot_aa:
        logger.info(
            YELLOW + BOLD + "Ommiting User aa Annotation, no AA files supplied" + END_FORMATTING)
    else:
        counter = 0
        l = len([file for file in os.listdir(out_annot_snpeff_dir) if file.endswith(".annot")]) - 1
        for root, _, files in os.walk(out_annot_snpeff_dir):
            if root == out_annot_snpeff_dir:
                for name in files:
                    if name.endswith('.annot'):
                        logger.info(GREEN + "Creating User aa annotation " + name + END_FORMATTING)
                        os.system("sbatch -J %s /home/laura/Laura_intel/Desktop/covid_multianalysis/useraa_annotation.sh %s %s %s %s %s > jobid.batch" 
                        %(name_s + "aa", name, sample, output, root, output + "/" + "annot_aa.txt"))
                        os.system('while [ "$(squeue | grep $USER | grep %s | wc -l)" = "96" ]; do sleep 0.1; done' %(name_s))
                        os.system('if [ %s = %s ]; then while [ $(squeue | grep $USER | grep "%s" | wc -l) != 0 ]; do sleep 0.1; done; fi' %(str(counter), l, name_s))
                        counter += 1

    os.system("rm %s %s %s" %(output + "/" + "annot_vcf.txt", output + "/" + "annot_aa.txt", output + "/" + "annot_bed.txt"))

    #LINAGE WITH PANGOLIN ###############################
    #####################################################
    # Annotate consensus fasta files to obtian sample linage using pangolin.
    # Output is located in Annotation/pangolin.
    counter = 0
    l = len([file for file in os.listdir(out_consensus_ivar_dir) if file.endswith(".fa")]) - 1
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
                            logger.info(GREEN + "Creating pangolin annotation " + name + END_FORMATTING)
                            # Annotate variants with pangolin to obtain the linage
                            os.system("sbatch -J %s /home/laura/Laura_intel/Desktop/covid_multianalysis/annotate_pangolin.sh %s %s %s %s %s > jobid.batch" 
                            %(name_s + "pan", filename, out_annot_pangolin_dir, out_pangolin_filename, "2", "0.6"))
                            os.system('while [ "$(squeue | grep $USER | grep %s | wc -l)" = "50" ]; do sleep 0.1; done' %(name_s))
                            os.system('if [ %s = %s ]; then while [ $(squeue | grep $USER | grep "%s" | wc -l) != 0 ]; do sleep 0.1; done; fi' %(str(counter), l, name_s))
                            counter += 1
    os.remove("jobid.batch")
    os.system("rm %s" %("slurm-*"))

    # USER AA TO HTML
    # Convert tsv user_aa annotation to html. Files are created in output/Annotation/user_aa.
    user_aa_to_html(out_annot_user_aa_dir, group_name)

    # SNP COMPARISON using tsv variant files
    ######################################################
    snp_comparison(name_s, logger, output, group_name, out_variant_ivar_dir, out_stats_coverage_dir)

# Parse arguments
args = get_arguments()

# Global variables
output = os.path.abspath(args.output)               # Output path (argument)
check_create_dir(output)
group_name = output.split("/")[-1]                  # Group name (output folder name)
reference = os.path.abspath(args.reference)         # Reference path (argument)
annotation = os.path.abspath(args.annotation)       # Annotation path (argument)

# LOGGING
# Create log file with date and time
logger = logging.getLogger()
create_logFile(group_name, logger, output, args)

# Set output folder structre
set_folder_structure(output)

# Extract name reads and already analysed reads
r1, r2, sample_list_F, new_samples = obtain_reads(args, logger)

# Run COVIDMA
covidma(output, args, logger, r1, r2, sample_list_F, new_samples, group_name, reference, annotation)
