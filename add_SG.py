#!/usr/bin/env python

import gzip
import sys
from misc import execute_subprocess

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

sample = sys.argv[1]
output_sorted_path = sys.argv[2]
output_bg_sorted_path = sys.argv[3]
r1 = sys.argv[4]

add_SG(sample, output_sorted_path, output_bg_sorted_path, r1)