#!/usr/bin/env python3

import sys
import os
import shutil
from pathlib import Path

accession = sys.argv[1]
ipd_avrl_dict = sys.argv[2]

genotype = accession + ".genotypes.csv"
allele_diag_list = accession + ".avrl_allele_list.tsv"
allele_fl_list = accession + ".allele_list.tsv"
with open(ipd_avrl_dict) as f_in:
	ipd_diag_dict = json.load(f_in)
import pandas as pd
df_diag = pd.read_csv(allele_diag_list, sep='\t')
df =pd.read_csv(allele_fl_list, sep='\t')
diag_allele_list = list(df_diag['allele'].unique())
ipd_allele_list = list(df['allele'].unique())
missing_diag_list = []
for diag_i in diag_allele_list:
	should_present_allele_list = ipd_diag_dict[diag_i]
	at_least_one_allele = any(item in should_present_allele_list for item in ipd_allele_list)
	if not at_least_one_allele:
		missing_diag_list.append(diag_i)
df.rename(columns={'depth':'read_ct'}, inplace=True)
df['read_ct'] = df['read_ct'].round(decimals=0)
df['accession'] = accession
if len(missing_diag_list) > 0 :
	df_diag_missing = pd.DataFrame({'allele':missing_diag_list})
	df_diag_missing = df_diag.merge(df_diag_missing, on=['allele'], how='inner')
	df_diag_missing.rename(columns={'depth':'read_ct'}, inplace=True)
	df_diag_missing['read_ct'] = df_diag_missing['read_ct'].round(decimals=0)
	df_diag_missing['read_ct'] = df_diag_missing['read_ct'] - .01
	df_diag_missing['accession'] = accession
	df = pd.concat([df, df_diag_missing], ignore_index=True)
df.to_csv(genotype, index=False)
