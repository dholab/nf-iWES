#!/usr/bin/env python3

import sys
import pandas as pd
import json
import os
import shutil
from pathlib import Path

ipd_avrl_dict = sys.argv[1]
ipd_num_lookup = sys.argv[2]
accession = sys.argv[3]
missing_alleles = sys.argv[4]

with open(ipd_avrl_dict) as f_in:
	ipd_diag_dict = json.load(f_in)
with open(ipd_avrl_dict) as f_in:
	ipd_num_dict = json.load(f_in)

sample = accession
# open the list of missing alleles do to no corresponding diag region to the fl.
df_missing = pd.read_csv(missing_alleles,sep='\t',header=None,names=['allele'])
missing_allele_list = list(df_missing['allele'].unique())
# ipd_num_dict
included = False
# open the list of alleles that past the depth/computational chimera filter
df = pd.read_csv(accession + ".allele_list.tsv",sep='\t')
# get a list of the unique alleles
diag_present_list = list(df['allele'].unique())
ipd_allele_list = []
# convert the list to the corresponding fl sequences, as many of them multi-map.
for diag_i in diag_present_list:
	if diag_i in ipd_diag_dict.keys():
		ipd_allele_list = ipd_allele_list + ipd_diag_dict[diag_i]

#  os.makedirs(exhaustive_result,exist_ok=True)
# add on the missing allele list
ipd_allele_list = ipd_allele_list + missing_allele_list
# write the two lists of alleles and corresponding number list for the fasta files.
with open(accession + ".allele_list_fl.txt",'w') as f:
	f.write('\n'.join(ipd_allele_list))
ipd_num_list = [ipd_num_dict[x] for x in ipd_allele_list]
with open(accession + ".allele_list_fl_num.txt",'w') as f:
	f.write('\n'.join(ipd_num_list))
