#!/usr/bin/env python

# Standard library imports
import os
import sys

# Local application imports
from bam_variant import ivar_variants

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
        print(out_ivar_variant_file +
                    " EXIST\nOmmiting Variant call for  sample " + sample )
    else:
        print( "Calling variants with ivar in sample " + sample )
        # Perform variant calling using ivar
        ivar_variants(reference, output_markdup_trimmed_file, out_variant_dir, sample,
                        annotation, min_quality=15, min_frequency_threshold=0.01, min_depth=1)

# Input variables
output = sys.argv[1]
output_markdup_trimmed_file = sys.argv[2]
sample = sys.argv[3]
reference = sys.argv[4]
annotation = sys.argv[5]

# Function
ivar_variant_calling(output, output_markdup_trimmed_file, sample, reference, annotation)