#!/usr/bin/env python

# Standard library imports
import os
import sys

# Local application imports
from bam_variant import create_coverage

def coverage_stats(output, sample, output_markdup_trimmed_file):
    """
    Function that compute some metrics realted to coverage
    from bam file trimmed and without duplicates, using samtools.

    Output is located in Stats/Coverage.
    """

    out_stats_dir = os.path.join(output, "Stats")                                   # Folder
    out_stats_coverage_dir = os.path.join(out_stats_dir, "Coverage")                # subfolder
    out_coverage_name = sample + ".cov"                                             # Name output file

    # Compute coverage metrics using samtools
    create_coverage(output_markdup_trimmed_file,
                    out_stats_coverage_dir, sample)

# Input variables
output = sys.argv[1]
sample = sys.argv[2]
output_markdup_trimmed_file = sys.argv[3]

# Function
coverage_stats(output, sample, output_markdup_trimmed_file)