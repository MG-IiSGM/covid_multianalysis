import sys
import os
import pandas as pd

def extract_uncovered(path, slf_files_f, indel_position_final_f, start, end):
    slf_files = []
    f = open(slf_files_f, "r")
    for l in f:
        slf_files.append(l.strip())
    f.close()

    indel_position_final = []
    f2 = open(indel_position_final_f, "r")
    for l in f2:
        indel_position_final.append(int(l.strip()))
    f2.close()
    part = slf_files[start:end]

    for slf in part:
        sample = slf.split(".")[0]
        df_slf = pd.read_csv(path + "/" + slf, sep='\t')
        df_uc = pd.read_csv(path + "/" + sample + ".ucov", sep="\t")
        df_uc = df_uc[~df_uc.POS.isin(indel_position_final)]
        df_slf[sample].update(df_slf[['REGION', 'POS']].merge(
                        df_uc, on=['REGION', 'POS'], how='left')[sample])
        df_slf.to_csv(path + "/" + sample + ".uslf", index=False, sep="\t")
        os.system("rm %s %s" %(path + "/" + slf, path + "/" + sample + ".ucov"))

path = sys.argv[1]
csv_files_f = sys.argv[2]
indel_position_final_f = sys.argv[3]
start = int(sys.argv[4])
end = int(sys.argv[5])

extract_uncovered(path, csv_files_f, indel_position_final_f, start, end)
