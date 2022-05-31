#!/usr/bin/env python

# Standard library imports
import os
import sys

# Local application imports
from bam_variant import ivar_consensus, replace_consensus_header

def consensus_ivar_create(output, sample, output_markdup_trimmed_file):
    """
    Function that creates a consensus file from bam file after
    trimming and removing duplicates.

    Output consensus file is in Consensus/ivar.
    """

    out_consensus_dir = os.path.join(output, "Consensus")                                       # Folder
    out_consensus_ivar_dir = os.path.join(out_consensus_dir, "ivar")                            # subfolder
    out_ivar_consensus_name = sample + ".fa"                                                    # consensus fasta name file
    out_ivar_consensus_file = os.path.join(out_consensus_ivar_dir, out_ivar_consensus_name)     # Abosolute path consensus file

    # Create a consensus file using ivar.
    ivar_consensus(output_markdup_trimmed_file, out_consensus_ivar_dir, sample,
                    min_quality=20, min_frequency_threshold=0.8, min_depth=20, uncovered_character='N')
    print( "Replacing consensus header in " + sample )
    # As header of fasta file is the sample name
    replace_consensus_header(out_ivar_consensus_file)

# Input variables
output = sys.argv[1]
sample = sys.argv[2]
output_markdup_trimmed_file = sys.argv[3]

# Function
consensus_ivar_create(output, sample, output_markdup_trimmed_file)