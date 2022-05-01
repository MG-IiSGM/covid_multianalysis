#!/usr/bin/env python

import os
import gzip
import argparse
#import argcomplete
import subprocess
from misc import check_file_exists, extract_sample, check_create_dir, execute_subprocess, check_remove_file


"""
=============================================================
HEADER
=============================================================

INSTITUTION:IiSGM
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
VERSION=0.1
CREATED: 27 March 2019
REVISION:
    20200923 - Revisited for covidma pipeline

TODO

================================================================
END_OF_HEADER
================================================================
"""

"""
def bowtie2_mapping(args):
    r1 = os.path.abspath(args.r1_file)
    r2 = os.path.abspath(args.r2_file)
    reference = os.path.abspath(args.reference)

    sample = extract_sample(r1,r2)
    output_dir = obtain_output_dir(args, "Bam")
    sample_name = sample + ".sam"
    output_file = os.path.join(output_dir, sample_name)

    check_create_dir(output_dir)

    if args.extensive_mapping:
        extensive_command = "-a"
    else:
        extensive_command = ""
    #bowtie2 index
    cmd_index = ["bowtie2-build", reference, reference]
    execute_subprocess(cmd_index)
    
    #bowtie map
    cmd_map = ["bowtie2", "-1", r1, "-2", r2, "-S", output_file, "-q", "--very-sensitive-local", "-p", str(args.threads), "-x", reference, extensive_command]
    execute_subprocess(cmd_map)
"""


def bwa_mapping(r1, r2, reference, sample, output_dir, threads=8):
    """
    #Store output in a file when it is outputted in stdout
    https://stackoverflow.com/questions/4965159/how-to-redirect-output-with-subprocess-in-python
    """
    r1 = os.path.abspath(r1)
    r2 = os.path.abspath(r2)
    reference = os.path.abspath(reference)

    sample_name = sample + ".sam"
    output_file = os.path.join(output_dir, sample_name)

    check_create_dir(output_dir)
    
    cmd_index = ["bwa", "index", reference]
    execute_subprocess(cmd_index)
    
    cmd_map = ["bwa", "mem", "-Y", "-M", "-t", str(threads), "-o", output_file, reference, r1, r2]
    execute_subprocess(cmd_map)
    """
    Create file whew it outputs thoug stdout --> Easier with -o param
    with open(output_file, "w") as outfile:
        #map reads and save it in th eoutput file
        subprocess.run(["bwa", "mem", "-t", str(args.threads), reference, r1, r2], 
        stdout=outfile, stderr=subprocess.PIPE, check=True, universal_newlines=True)
    """

def add_SG(sample, input_bam, output_bg_sorted, r1):
    """
    @MN00227:45:000H255J3:1:11102:21214:1110 1:N:0:18
    @NS500454:48:HKG57BGXX:1:11101:17089:1032 2:N:0:TCCTGAGC+TCTTACGC
    @NS500454:27:HJJ32BGXX:1:11101:12392:1099 1:N:0:2

    @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:
    <is filtered>:<control number>:<sample number | barcode1'+barcode2'>
    ID = Read group identifier {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE} 
    PU = Platform Unit #optional
    SM = Sample
    PL = Platform/technology used to produce the read (ILLUMINA, SOLID, LS454, HELICOS and PACBIO)
    LB = DNA preparation library identifier
    """

    with gzip.open(r1) as f:
        first_line = f.readline().strip().decode()
    #print(first_line)
    first_line_list = first_line.split(":")
    if len(first_line_list) > 4:
        rg_id = ".".join([first_line_list[2],first_line_list[3],first_line_list[-1]])
        rg_pu = ".".join([first_line_list[2],first_line_list[3],first_line_list[-1]])
    else:
        rg_id = first_line_list[0]
        rg_pu = first_line_list[0]
    rg_sm = sample
    rg_pl = "ILLUMINA"
    rg_lb = "lib_" + sample

    rg_id_param = "RGID=" + rg_id
    rg_pu_param = "RGPU=" + rg_pu
    rg_sm_param = "RGSM=" + rg_sm
    rg_pl_param = "RGPL=" + rg_pl
    rg_lb_param = "RGLB=" + rg_lb

    #picard_jar = get_picard_path()

    input_param = "INPUT=" + input_bam
    output_param = "OUTPUT=" + output_bg_sorted


    # java -jar picard.jar AddOrReplaceReadGroups \ 
    # INPUT=reads.bam \ OUTPUT=reads_addRG.bam \ RGID=H0164.2 \ #be sure to change from default of 1
    # RGLB= library1 \ RGPL=illumina \ RGPU=H0164ALXX140820.2 \ RGSM=sample1 \ 
    # SORT_ORDER=coordinate \ CREATE_INDEX=true

    cmd = ["picard", "AddOrReplaceReadGroups", 
    input_param, output_param, rg_id_param, rg_lb_param, rg_pl_param, rg_pu_param, rg_sm_param,
    "SORT_ORDER=coordinate"]
    execute_subprocess(cmd)

def sam_to_index_bam(sample, output_dir, r1, threads):
    # input_sam_path = os.path.abspath(input_sam)
    # if output_bam == "inputdir":
    #     output_bam = os.path.dirname(input_sam_path)
    # else:
    #     output_bam = output_bam

    sample_name = sample + ".sam"
    input_sam_path = os.path.join(output_dir, sample_name)

    input_name = (".").join(os.path.basename(input_sam_path).split(".")[:-1])

    output_bam_name = input_name + ".bam"
    output_bam_path = os.path.join(output_dir, output_bam_name)

    output_sorted_name = input_name + ".sorted.bam"
    output_sorted_path = os.path.join(output_dir, output_sorted_name)

    output_bg_sorted_name = input_name + ".rg.sorted.bam"
    output_bg_sorted_path = os.path.join(output_dir, output_bg_sorted_name)

    cmd_view = ["samtools", "view", "-Sb", input_sam_path, "--threads", str(threads), "-o", output_bam_path,]
    execute_subprocess(cmd_view)

    check_remove_file(input_sam_path)
    
    cmd_sort = ["samtools", "sort", output_bam_path, "-o", output_sorted_path]
    execute_subprocess(cmd_sort)

    check_remove_file(output_bam_path)

    add_SG(sample, output_sorted_path, output_bg_sorted_path, r1)

    check_remove_file(output_sorted_path)

    """
    #samtools index: samtools index $output_dir/$sample".sorted.bam"
    subprocess.run(["samtools", "index", output_sorted_path], 
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    """
