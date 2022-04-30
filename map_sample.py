
# Standard library imports
import os
import sys


# Local application imports
from misc import check_file_exists, extract_sample
from preprocessing import fastqc_quality, fastp_trimming
from pe_mapper import bwa_mapping, sam_to_index_bam
from bam_variant import picard_markdup, ivar_trim, ivar_variants, ivar_consensus, \
    replace_consensus_header, create_bamstat, create_coverage
from vcf_process import filter_tsv_variants


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

def check_quality(r1_file, r2_file, output, name_folder, sub_folder_name, sample):
    """
    Function that check reads quality with fastqc
    Output is located in output/name_folder/sub_folder_name
    """

    # Set name files
    out_dir = os.path.join(output, name_folder)                                 # folder
    out_name_r1 = r1_file.replace(".fastq.gz", "_fastqc.html").split("/")[-1]   # fastq_r1 file name
    out_name_r2 = r2_file.replace(".fastq.gz", "_fastqc.html").split("/")[-1]   # fastq_r2 file name
    out_subdir = os.path.join(out_dir, sub_folder_name)                         # subfolder
    output_qc_file_r1 = os.path.join(out_subdir, out_name_r1)                   # absolute path r1
    output_qc_file_r2 = os.path.join(out_subdir, out_name_r2)                   # absolute path r2

    # Check if quality have been already computed
    if os.path.isfile(output_qc_file_r1) and os.path.isfile(output_qc_file_r2):
        print(YELLOW + DIM + output_qc_file_r1 +
                    " EXIST\nOmmiting QC for sample " + sample + END_FORMATTING)
    else:
        print(
            GREEN + "Checking quality in sample " + sample + END_FORMATTING)
        print("R1: " + r1_file + "\nR2: " + r2_file)
        # Check quality with fastqc
        fastqc_quality(r1_file, r2_file, out_subdir, 2)

def trim_read(r1_file, r2_file, sample, output):
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

    # Check if trimming has already been performed
    if os.path.isfile(output_trimming_file_r1) and os.path.isfile(output_trimming_file_r2):
        print(YELLOW + DIM + output_trimming_file_r1 +
                    " EXIST\nOmmiting Trimming for sample " + sample + END_FORMATTING)
    else:
        print(GREEN + "Trimming sample " +
                    sample + END_FORMATTING)
        # Use fastp for trimming
        fastp_trimming(r1_file, r2_file, sample, out_trim_dir, threads=2, 
            min_qual=20, window_size=10, min_len=35)
    
    # Return path of trimmed files
    return (output_trimming_file_r1, output_trimming_file_r2)

def mapping(output_trimming_file_r1, output_trimming_file_r2, sample, output, reference):
    """
    Function that maps reads in reference genome using BWA
    Then, transform sam file to sorted bam file.

    Output is located in Bam directory
    """

    out_map_dir = os.path.join(output, "Bam")                       # folder
    out_map_name = sample + ".rg.sorted.bam"                        # file name
    output_map_file = os.path.join(out_map_dir, out_map_name)       # absolute path file name

    # Check if mapping has already been performed
    if os.path.isfile(output_map_file):
        print(YELLOW + DIM + output_map_file +
                    " EXIST\nOmmiting Mapping for sample " + sample + END_FORMATTING)
    else:
        print(GREEN + "Mapping sample " +
                    sample + END_FORMATTING)
        print("R1: " + output_trimming_file_r1 + "\nR2: " +
                    output_trimming_file_r2 + "\nReference: " + reference)
        
        # Mapping to reference using BWA. As output a sam file is created
        bwa_mapping(output_trimming_file_r1, output_trimming_file_r2,
                    reference, sample, out_map_dir, threads=2)
        # Converts sam file to bam file using samtools and sort alingments
        sam_to_index_bam(sample, out_map_dir, output_trimming_file_r1, threads=2)
    
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

    # Check if mark duplicates has already been performed
    if os.path.isfile(output_markdup_file):
        print(YELLOW + DIM + output_markdup_file +
                    " EXIST\nOmmiting Duplucate Mark for sample " + sample + END_FORMATTING)
    else:
        print(GREEN + "Marking Dupes in sample " +
                    sample + END_FORMATTING)
        print("Input Bam: " + output_map_file)
        # Mark and delete duplicates with picard tools
        picard_markdup(output_map_file) 
    
    # Return absolute path marked and remove duplicates file
    return (output_markdup_file)

def trim_ivar(output_markdup_trimmed_file, output_markdup_file, sample, primers):
    """
    Function that trim reamining primers of aligned reads using ivar.
    Then uses samtools to sort reads and create an index to access 
    overlapping alignments quickly. 

    Output is located in Bam directory.
    """

    # Check if output_markdup_trimmed_file exits
    # If True, trimming is not performed
    if os.path.isfile(output_markdup_trimmed_file):
        print(YELLOW + DIM + output_markdup_trimmed_file +
                    " EXIST\nOmmiting Duplucate Mark for sample " + sample + END_FORMATTING)
    else:
        print(
            GREEN + "Trimming primers in sample " + sample + END_FORMATTING)
        print("Input Bam: " + output_markdup_file)
        # Trimming of remaining primers is performed with ivar
        ivar_trim(output_markdup_file, primers, sample,
                    min_length=30, min_quality=20, sliding_window_width=4)

def ivar_variant_calling(output, output_markdup_trimmed_file, sample, reference, annotation):
    """
    Function that performs variant calling using ivar

    Output is located in Variants/ivar_raw.
    """

    out_variant_dir = os.path.join(output, "Variants")                                      # Folder
    out_ivar_variant_name = sample + ".tsv"                                                 # Variant calling outup file name
    out_variant_ivar_dir = os.path.join(out_variant_dir, "ivar_raw")                        # subfolder
    out_ivar_variant_file = os.path.join(out_variant_ivar_dir, out_ivar_variant_name)

    # Check if variant calling has already been performed
    if os.path.isfile(out_ivar_variant_file):
        print(YELLOW + DIM + out_ivar_variant_file +
                    " EXIST\nOmmiting Variant call for  sample " + sample + END_FORMATTING)
    else:
        print(
            GREEN + "Calling variants with ivar in sample " + sample + END_FORMATTING)
        # Perform variant calling using ivar
        ivar_variants(reference, output_markdup_trimmed_file, out_variant_dir, sample,
                        annotation, min_quality=15, min_frequency_threshold=0.01, min_depth=1)
    
    # Return vcf file name
    return (out_ivar_variant_file)

def variant_filtering(output, sample, out_ivar_variant_file):
    """
    Filter variants that does not satisfy:
        * min_frequency=0.7
        * min_total_depth=10
        * min_alt_dp=4
        * is_pass=True
        * only_snp=False
    
    Output is located in Variants/ivar_filtered.
    """

    out_variant_dir = os.path.join(output, "Variants")                                      # Folder
    out_ivar_variant_name = sample + ".tsv"                                                 # vcf file name
    out_filtered_ivar_dir = os.path.join(out_variant_dir, "ivar_filtered")                  # subfolder
    out_ivar_filtered_file = os.path.join(out_filtered_ivar_dir, out_ivar_variant_name)     # Absolute path file name

    # Check if variant filtering has already been performed
    # If True, filtering is skipped
    if os.path.isfile(out_ivar_filtered_file):
        print(YELLOW + DIM + out_ivar_filtered_file +
                    " EXIST\nOmmiting Variant filtering for  sample " + sample + END_FORMATTING)
    else:
        print(GREEN + "Filtering variants in sample " +
                    sample + END_FORMATTING)
        # Filter variants checking in pandas DataFrame 
        filter_tsv_variants(out_ivar_variant_file, out_filtered_ivar_dir, min_frequency=0.7,
                            min_total_depth=10, min_alt_dp=4, is_pass=True, only_snp=False)

def consensus_create(output, sample, output_markdup_trimmed_file):
    """
    Function that creates a consensus file from bam file after
    trimming and removing duplicates.

    Output consensus file is in Consensus/ivar.
    """

    out_consensus_dir = os.path.join(output, "Consensus")                                       # Folder                                                   # Check if folder exists
    out_consensus_ivar_dir = os.path.join(out_consensus_dir, "ivar")                            # subfolder
    out_ivar_consensus_name = sample + ".fa"                                                    # consensus fasta name file
    out_ivar_consensus_file = os.path.join(out_consensus_ivar_dir, out_ivar_consensus_name)     # Abosolute path consensus file

    # Check if consesnsus file already exists
    # If True, a consensus file is not created
    if os.path.isfile(out_ivar_consensus_file):
        print(YELLOW + DIM + out_ivar_consensus_file +
                    " EXIST\nOmmiting Consensus for  sample " + sample + END_FORMATTING)
    else:
        print(
            GREEN + "Creating consensus with ivar in sample " + sample + END_FORMATTING)
        # Create a consensus file using ivar.
        ivar_consensus(output_markdup_trimmed_file, out_consensus_ivar_dir, sample,
                        min_quality=20, min_frequency_threshold=0.8, min_depth=20, uncovered_character='N')
        print(
            GREEN + "Replacing consensus header in " + sample + END_FORMATTING)
        # As header of fasta file is the sample name
        replace_consensus_header(out_ivar_consensus_file)

def bamstats(output, sample, output_markdup_trimmed_file):
    """
    Function that compute some metrics from bam file
    trimmed and without duplicates, using samtools.

    Output is located in Stats/Bamstats.
    """

    out_stats_dir = os.path.join(output, "Stats")                                       # Folder
    out_stats_bamstats_dir = os.path.join(out_stats_dir, "Bamstats")                    # subfolder
    out_bamstats_name = sample + ".bamstats"                                            # Output filename
    out_bamstats_file = os.path.join(out_stats_bamstats_dir, out_bamstats_name)         # absolute path to filename

    # Check if Bam stats have already been computed
    # If True it is skipped.
    if os.path.isfile(out_bamstats_file):
        print(YELLOW + DIM + out_bamstats_file +
                    " EXIST\nOmmiting Bamstats for  sample " + sample + END_FORMATTING)
    else:
        print(GREEN + "Creating bamstats in sample " +
                    sample + END_FORMATTING)
        # Compute metrics from bam file using samtools
        create_bamstat(output_markdup_trimmed_file,
                        out_stats_bamstats_dir, sample, threads=2)

def coverage_stats(output, sample, output_markdup_trimmed_file):
    """
    Function that compute some metrics realted to coverage
    from bam file trimmed and without duplicates, using samtools.

    Output is located in Stats/Coverage.
    """

    out_stats_dir = os.path.join(output, "Stats")                                   # Folder
    out_stats_coverage_dir = os.path.join(out_stats_dir, "Coverage")                # subfolder
    out_coverage_name = sample + ".cov"                                             # Name output file
    out_coverage_file = os.path.join(out_stats_coverage_dir, out_coverage_name)     # Absolute output path name

    # Check if metrics have already been computed
    # If True, this part is skipped
    if os.path.isfile(out_coverage_file):
        print(YELLOW + DIM + out_coverage_file +
                    " EXIST\nOmmiting Bamstats for  sample " + sample + END_FORMATTING)
    else:
        print(GREEN + "Creating coverage in sample " +
                    sample + END_FORMATTING)
        # Compute coverage metrics using samtools
        create_coverage(output_markdup_trimmed_file,
                        out_stats_coverage_dir, sample)

def map_sample(output, primers, r1_file, r2_file, reference, annotation):
    """
    Function that maps reads to reference genome and performs variant calling.
    """
    sample = extract_sample(r1_file, r2_file)

    out_map_dir = os.path.join(output, "Bam")                                               # Folder
    out_markdup_trimmed_name = sample + ".rg.markdup.trimmed.sorted.bam"                    # Filname
    output_markdup_trimmed_file = os.path.join(out_map_dir, out_markdup_trimmed_name)       # absolute path to filename

    # Check if sample has been already analysed
    # If True, sample does not need to be analysed again
    print("\n" + WHITE_BG + "STARTING SAMPLE: " + sample + END_FORMATTING)

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
        check_quality(r1_file, r2_file, output, "Quality", "raw", sample)

        # QUALITY TRIMMING AND ADAPTER REMOVAL WITH fastp
        ###################################################
        # Trim reads by window. Trim regions or reads that does not satisfy min_qual=20, window_size=10, min_len=35
        # by using fastp. Output is in output/Trimmed
        output_trimming_file_r1, output_trimming_file_r2 = trim_read(r1_file, r2_file, sample, output)

        # QUALITY CHECK in TRIMMED with fastqc
        ######################################################
        # Check quality of trimmed fastq with fastqc and store info in output/Quality/processed
        check_quality(output_trimming_file_r1, output_trimming_file_r2, output, "Quality", "processed", sample)

        # MAPPING WITH BWA - SAM TO SORTED BAM - ADD HEADER SG
        #####################################################
        # Map trimmed reads to reference genome. As output a sorted bam file is generated in output/Bam
        output_map_file = mapping(output_trimming_file_r1, output_trimming_file_r2, sample, output, reference)

        #MARK DUPLICATES WITH PICARDTOOLS ###################
        #####################################################
        # Mark and remove duplicates in bam file. Previous sorted bam file in output/Bam is sustituted
        # by a sorted without duplicates bam file.
        output_markdup_file = mark_duplicates(output_map_file, sample, output)

        #TRIM PRIMERS WITH ivar trim ########################
        #####################################################
        # Trim aligned reads (bam file) using ivar. Output file is in output/Bam.
        trim_ivar(output_markdup_trimmed_file, output_markdup_file, sample, primers)

    else:
        print(YELLOW + DIM + output_markdup_trimmed_file +
                " EXIST\nOmmiting BAM mapping and BAM manipulation in sample " + sample + END_FORMATTING)

    ########################END OF MAPPING AND BAM MANIPULATION##########################
    #####################################################################################

    #VARIANT CALLING WTIH ivar variants##################
    #####################################################
    # Variant calling with ivar. Output is located in output/Variants/ivar_raw
    out_ivar_variant_file = ivar_variant_calling(output, output_markdup_trimmed_file, sample, reference, annotation)

    #VARIANT FILTERING ##################################
    #####################################################
    # Filter variants detected with ivar. Output is located in output/Variants/ivar_filtered
    variant_filtering(output, sample, out_ivar_variant_file)

    #CREATE CONSENSUS with ivar consensus##################
    #######################################################
    # Create a consensus fasta file from bam file trimmed and without duplicates.
    # Output is located in output/Consensus/ivar.
    consensus_create(output, sample, output_markdup_trimmed_file)

    ########################CREATE STATS AND QUALITY FILTERS###############################
    #######################################################################################
    #CREATE Bamstats #######################################
    ########################################################
    # Compute metrics from trimmed and without duplicates bam file using samtools.
    # Output is located in output/Stats/Bamstats.
    bamstats(output, sample, output_markdup_trimmed_file)

    #CREATE CoverageStats ##################################
    ########################################################
    coverage_stats(output, sample, output_markdup_trimmed_file)


output = sys.argv[0]
primers = sys.argv[1]
r1_file = sys.argv[2]
r2_file = sys.argv[3]
reference = sys.argv[4]
annotation = sys.argv[5]

map_sample(output, primers, r1_file, r2_file, reference, annotation)
