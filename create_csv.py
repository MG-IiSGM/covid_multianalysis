#!/usr/bin/env python

# imports
import pandas as pd
import sys

def create_csv(df_file, path, csv_files_f, start, end):
    csv_files = []
    f = open(csv_files_f, "r")
    for l in f:
        csv_files.append(l.strip())
    f.close()
    part = csv_files[start:end]

    # df = pd.read_csv(df_file, sep = "\t")
    df = pd.read_hdf(df_file, mode="r")

    for csv in part:
        df_sample = df[['REGION', 'POS', 'REF', 'ALT', csv]]
        # df_sample.to_csv(path + "/" + csv + ".csv", index=False, sep="\t")
        df_sample.to_hdf(path + "/" + csv + ".csv", "hdf", mode="w", format="fixed", index=False)

df_file = sys.argv[1]
path = sys.argv[2]
csv_files_f = sys.argv[3]
start = int(sys.argv[4])
end = int(sys.argv[5])

create_csv(df_file, path, csv_files_f, start, end)
