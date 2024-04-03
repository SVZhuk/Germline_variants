#!/usr/bin/env /kuacc/users/szhuk/.conda/envs/ds/bin/python

import pandas as pd

bed_file_path = "/kuttam_fg/refdata/zhuk/Probes/SureSelectAllExonV6r2_100bp_hg38.sorted.bed"
df = pd.read_csv(bed_file_path, 
                 sep="\t", 
                 header=None)

genes = ["ALS2", "ANG", "C19orf12", "CCNF", "CHCHD10", "CHMP2B", "CHRNA3","CHRNA4",
         "CHRNB4", "CREST", "DAO", "DCTN1", "ELP3", "ERBB4", "EWSR1", "FIG4", "FUS",
         "hnRNPA1", "hnRNPA2B1", "MATR3", "NEFH", "NEK1", "OPTN", "PFN1", "PNPLA6", 
         "PON1", "PON3", "PRPH", "SETX", "SIGMAR1", "SOD1", "SPG11", "SQSTM1", "TAF15", 
         "TARDBP", "TBK1", "TUBA4A", "UBQLN2", "VAPB", "VCP", "CFAP410", "ANXA11", 
         "KIF5A", "ERLIN1", "GLT8D1", "DNAJC7", "SPTLC1", "RPRD2"]

df[4] = df[3].str.split(",")
df = df[df[4].apply(lambda x: any(gene.split("|")[-1] in genes for gene in x))]
df = df.iloc[:, :4]

df_sub = df[3].str.split(",")
first_elements = df_sub.apply(lambda x: x[0].split("|")[-1])
num_unique = first_elements.nunique()

df_length = df[2] - df[1]
total_coverage = df_length.sum()

print(f"Number of genes in the list: {len(genes)}")
print(f"Number of unique first elements: {num_unique}")
print(f"Probes coverage: {total_coverage}bp")

print(df.head())

df.to_csv("/kuttam_fg/refdata/zhuk/Probes/ATX_screen.subset.bed", 
          sep="\t", 
          header=None, 
          index=False)
