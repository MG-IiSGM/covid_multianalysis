#!/usr/bin/env python

# Standard library imports
from fileinput import filename
import os
import sys


# Local application imports
from annotation import user_annotation

def annotate_user(name, root, output, sample, f_annot_vcf, f_annot_bed):
    """
    Function that annotate variants from user files.

    Output is located in Annotation/user.
    """

    # list annot_vcf
    annot_vcf = []
    f = open(f_annot_vcf, "r")
    for s in f:
        annot_vcf.append(s.strip())
    f.close()

    # list annot_bed
    annot_bed = []
    f = open(f_annot_bed, "r")
    for s in f:
        annot_bed.append(s.strip())
    f.close()

    out_annot_dir = os.path.join(output, "Annotation")          # folder
    out_annot_user_dir = os.path.join(out_annot_dir, "user")    # subfolder
    sample = name.split('.')[0]
    print('User bed/vcf annotation in sample {}'.format(sample))
    filename = os.path.join(root, name)
    out_annot_file = os.path.join(
        out_annot_user_dir, sample + ".tsv")
    # Perform user annotation
    user_annotation(filename, out_annot_file, vcf_files=annot_vcf, bed_files=annot_bed)

name = sys.argv[1]
root = sys.argv[2]
output = sys.argv[3]
sample = sys.argv[4]
f_annot_vcf = sys.argv[5]
f_annot_bed = sys.argv[6]

annotate_user(name, root, output, sample, f_annot_vcf, f_annot_bed)