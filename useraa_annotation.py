#!/usr/bin/env python

# Standard library imports
import os
import sys

from annotation import user_annotation_aa

def useraa_annotation(name, sample, output, root, f_annot_aa):
    """
    Function that annotate variants from user aa files.

    Output is located in Annotation/user_aa.
    """

    # list annot_vcf
    annot_aa = []
    f = open(f_annot_aa, "r")
    for s in f:
        annot_aa.append(s.strip())
    f.close()

    out_annot_dir = os.path.join(output, "Annotation")              # Folder
    out_annot_user_aa_dir = os.path.join(out_annot_dir, "user_aa")  # subfolder
    sample = name.split('.')[0]
    print( 'User aa annotation in sample {}'.format(sample))
    filename = os.path.join(root, name)
    out_annot_aa_file = os.path.join(
        out_annot_user_aa_dir, sample + ".tsv")
    # Perform user annotation
    if os.path.isfile(out_annot_aa_file):
        user_annotation_aa(
            out_annot_aa_file, out_annot_aa_file, aa_files=annot_aa)
    else:
        user_annotation_aa(
            filename, out_annot_aa_file, aa_files=annot_aa)

name = sys.argv[1]
sample = sys.argv[2]
output = sys.argv[3]
root = sys.argv[4]
f_annot_aa = sys.argv[5]

useraa_annotation(name, sample, output, root, f_annot_aa)