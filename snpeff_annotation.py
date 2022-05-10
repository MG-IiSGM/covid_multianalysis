#!/usr/bin/env python

# Standard library imports
import os
import sys


# Local application imports
from annotation import annotate_snpeff


def snpeff_annotation(output, name, root, sample, snpeff_database):
    """
    Function that annotates variants using snpEff.

    Annotations are stored in Annotation/snpeff.
    """

    out_annot_dir = os.path.join(output, "Annotation")              # Folder
    out_annot_snpeff_dir = os.path.join(out_annot_dir, "snpeff")    # subfolder
    sample = name.split('.')[0]                                     # sample
    filename = os.path.join(root, name)                             # sample name
    out_annot_file = os.path.join(out_annot_snpeff_dir, sample + ".annot")  # Absolute path file name

    # Check if annotation has already been performed
    # If True, annotation is skipped
    if os.path.isfile(out_annot_file):
        print( out_annot_file +
                    " EXIST\nOmmiting snpEff Annotation for sample " + sample )
    else:
        print( "Annotating sample with snpEff: " + sample )
        output_vcf = os.path.join(
            out_annot_snpeff_dir, sample + '.vcf')
        # Annotates variants using snpEff
        annotate_snpeff(filename, output_vcf, out_annot_file, database=snpeff_database)

output = sys.argv[1]
name = sys.argv[2]
root = sys.argv[3]
sample = sys.argv[4]
snpeff_database = sys.argv[5]

snpeff_annotation(output, name, root, sample, snpeff_database)