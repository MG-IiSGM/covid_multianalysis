import sys
import os
import pandas as pd

def extract_lowfreq(path, csv_files_f, start, end):
    csv_files = []
    f = open(csv_files_f, "r")
    for l in f:
        csv_files.append(l.strip())
    f.close()
    part = csv_files[start:end]

    for csv in part:
        sample = csv.split(".")[0]
        if not os.path.exists(path + "/" + sample + ".lf"):
            os.system("mv %s %s" %(path + "/" + csv, path + "/" + sample + ".slf"))
            continue
        df_csv = pd.read_csv(path + "/" + csv, sep='\t')
        df_lf = pd.read_csv(path + "/" + sample + ".lf", sep="\t")
        df_csv[sample].update(df_csv[['REGION', 'POS', 'REF', 'ALT']].merge(
                        df_lf, on=['REGION', 'POS', 'REF', 'ALT'], how='left')[sample])
        df_csv.to_csv(path + "/" + sample + ".slf", index=False, sep="\t")
        os.system("rm %s %s" %(path + "/" + csv, path + "/" + sample + ".lf"))

path = sys.argv[1]
csv_files_f = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])

extract_lowfreq(path, csv_files_f, start, end)
