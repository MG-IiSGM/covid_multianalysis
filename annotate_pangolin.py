#!/usr/bin/env python

# Standard library imports
import os
import sys

from annotation import annotate_pangolin


input_file = sys.argv[1]
output_folder = sys.argv[2]
output_filename = sys.argv[3]
threads = int(sys.argv[4])
max_ambig = float(sys.argv[5])

annotate_pangolin(input_file, output_folder, output_filename, threads, max_ambig)