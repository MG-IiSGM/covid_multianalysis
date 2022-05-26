#!/usr/bin/env python

# imports
import pandas as pd
import os
import sys

# Auxiliary function
def import_tsv_variants(tsv_file, cov_path,  min_total_depth=4, min_alt_dp=4, only_snp=True):

    base_file = os.path.basename(tsv_file)
    input_file = os.path.abspath(tsv_file)
    sample = base_file.split(".")[0]

    df = pd.read_csv(input_file, sep='\t')
    df = df.drop_duplicates(subset=['POS', 'REF', 'ALT'], keep="first")

    df_lowfreq = df[(df.ALT_DP < 4) &
            (df.ALT_FREQ >= 0.7)]
    
    df_uncover = pd.read_csv(cov_path + "/" + sample + ".cov", sep="\t", header=None)
    df_uncover.columns = ['REGION', 'POS', sample]
    df_uncover = df_uncover[df_uncover[sample] == 0]
    df_uncover = df_uncover.replace(0, '!')

    df = df[((df.TOTAL_DP >= min_total_depth) &
             (df.ALT_DP >= min_alt_dp))]
    
    df = df[['REGION', 'POS', 'REF', 'ALT', 'ALT_FREQ']]
    df_lowfreq = df_lowfreq[['REGION', 'POS', 'REF', 'ALT', 'ALT_FREQ']]
    df_lowfreq['ALT_FREQ'] = '?'

    df = df.rename(columns={'ALT_FREQ': sample})
    df_lowfreq = df_lowfreq.rename(columns={'ALT_FREQ': sample})

    if only_snp == True:
        df = df[~(df.ALT.str.startswith('+') | df.ALT.str.startswith('-'))]
        df_lowfreq = df_lowfreq[~(df_lowfreq.ALT.str.startswith('+') | df_lowfreq.ALT.str.startswith('-'))]
        #print(df_lowfreq)

        return (df, df_lowfreq, df_uncover)
    else:
        #df[sample].update(df[['REGION', 'POS', 'REF', 'ALT']].merge(
        #                df_lowfreq, on=['REGION', 'POS', 'REF', 'ALT'], how='left')[sample])
        return (df, df_lowfreq, df_uncover)

def merge_df(path, tsv_files_f, start, end, old_tag, new_tag, flag, cov_path, out_compare_dir):

    c = 0
    tsv_files = []
    f = open(tsv_files_f, "r")
    for l in f:
        tsv_files.append(l.strip())
    f.close()
    part = tsv_files[start:end]
    for file in part:

        if not c and flag == 0:
            df, df_lowfreq, df_uncover = import_tsv_variants(path + "/" + file, cov_path, only_snp=False)
            #df = import_tsv_variants(path + "/" + file, only_snp=False)

            if df_lowfreq.shape[0]: # If it is not empty
                df_lowfreq.to_csv(out_compare_dir + "/" + file.split(".")[0] + ".lf", index=False, sep="\t")
            df_uncover.to_csv(out_compare_dir + "/" + file.split(".")[0] + ".ucov", index=False, sep="\t")
            c += 1
            continue

        elif not c and flag:
            df = pd.read_csv(path + "/" + file, sep="\t")
            c += 1
            continue

        if flag == 0 and file.endswith(old_tag):
            dfv, df_lowfreq, df_uncover = import_tsv_variants(path + "/" + file, cov_path, only_snp=False)
            #dfv = import_tsv_variants(path + "/" + file, only_snp=False)
            if df_lowfreq.shape[0]: # If it is not empty
                df_lowfreq.to_csv(out_compare_dir + "/" + file.split(".")[0] + ".lf", index=False, sep="\t")
            df_uncover.to_csv(out_compare_dir + "/" + file.split(".")[0] + ".ucov", index=False, sep="\t")
            df = df.merge(dfv, how="outer")

        elif flag and file.endswith(old_tag):
            dfv = pd.read_csv(path + "/" + file, sep = "\t")
            df = df.merge(dfv, how="outer")

    df.to_csv(out_compare_dir + "/" + "%i-%i%s" %(start, end, new_tag), index=False, sep='\t')
    if old_tag != ".tsv":
        for t in part:
            os.system("rm %s" %(path + "/" + t))

path = sys.argv[1]
tsv_files_f = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])
old_tag = sys.argv[5]
new_tag = sys.argv[6]
flag = int(sys.argv[7])
cov_path = sys.argv[8]
out_compare_dir = sys.argv[9]

merge_df(path, tsv_files_f, start, end, old_tag, new_tag, flag, cov_path, out_compare_dir)