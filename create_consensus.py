#!/usr/bin/env python

from Bio import SeqIO
import pandas as pd
import os 
import sys

def create_consensus(reference, highfreq_df, coverage_folder, out_folder):
    
    # Read csv
    dfsample = pd.read_csv(highfreq_df, sep="\t")
    sample = dfsample.columns[-1]
    coverage_folder = os.path.abspath(coverage_folder)
    cov_file = os.path.join(coverage_folder, sample + ".cov")
    reflist = list(SeqIO.read(reference, "fasta").seq)
    
    for _, row in dfsample.iterrows():
        if str(row[sample]) == '1':
            postition_list = row.Position.split("|")
            ref = postition_list[1]
            pos = int(postition_list[2])
            alt = postition_list[3]
            if reflist[pos - 1] == ref:
                reflist[pos - 1] = alt
    covdf = pd.read_csv(cov_file, sep="\t", names=["#CHROM", "POS", "COV"])
    uncovered = covdf[covdf.COV == 0]
    for _, row in uncovered.iterrows():
        reflist[row.POS - 1] = 'N'
    output_file = os.path.join(out_folder, sample + ".consensus.fasta")
    with open(output_file, 'w+') as fout:
        fout.write('>{}\n{}\n'.format(sample,('').join(reflist)))
    os.remove(highfreq_df)

reference = sys.argv[1]
highfreq_df = sys.argv[2]
coverage_folder = sys.argv[3]
out_folder = sys.argv[4]

create_consensus(reference, highfreq_df, coverage_folder, out_folder)
