#!/usr/bin/env python

# Standard library imports
import os
import sys


# Local application imports
from vcf_process import filter_tsv_variants

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
        print( out_ivar_filtered_file +
                    " EXIST\nOmmiting Variant filtering for  sample " + sample )
    else:
        print( "Filtering variants in sample " + sample)
        # Filter variants checking in pandas DataFrame 
        filter_tsv_variants(out_ivar_variant_file, out_filtered_ivar_dir, min_frequency=0.7,
                            min_total_depth=10, min_alt_dp=4, is_pass=True, only_snp=False)

# Input variables
output = sys.argv[1]
sample = sys.argv[2]
out_ivar_variant_file = sys.argv[3]

# Function
variant_filtering(output, sample, out_ivar_variant_file)