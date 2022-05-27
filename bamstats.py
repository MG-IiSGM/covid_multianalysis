#!/usr/bin/env python

# Standard library imports
import os
import sys


# Local application imports

from bam_variant import create_bamstat


def bamstats(output, sample, output_markdup_trimmed_file, threads):
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
        print( out_bamstats_file +
                    " EXIST\nOmmiting Bamstats for  sample " + sample )
    else:
        print( "Creating bamstats in sample " + sample )
        # Compute metrics from bam file using samtools
        create_bamstat(output_markdup_trimmed_file,
                        out_stats_bamstats_dir, sample, threads=threads)


# Input variables
output = sys.argv[1]
sample = sys.argv[2]
output_markdup_trimmed_file = sys.argv[3]
threads = int(sys.argv[4])

# Function
bamstats(output, sample, output_markdup_trimmed_file, threads)