#!/usr/bin/env python
from argparse import ArgumentParser, ArgumentTypeError
import pandas as pd
import os
import numpy as np
import json
import math
import os


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ArgumentTypeError('Boolean value expected.')


def parse_process_arrays_args(parser: ArgumentParser):
    """Parses the python script arguments from bash and makes sure files/inputs are valid"""
    parser.add_argument('--out_dir',
                        type=str,
                        help='directory to output the results',
                        required=True)
    parser.add_argument('--project_name',
                        type=str,
                        help='Project name for file prefix',
                        required=True)
    parser.add_argument('--config_dir',
                        type=str,
                        help='where the config files are stored',
                        default=None,
                        required=False)
    parser.add_argument('--animal_lookup_path',
                        type=str,
                        help='2 column csv file each row is a file prefix to corresponding sample id',
                        required=True)


    parser.add_argument('--bait_fasta',
                        type=str,
                        help='Where you ipd-diag combined fasta exists',
                        default=None,
                        required=False)
    parser.add_argument('--haplotype_lookup',
                        type=str,
                        help='you haplotype look up created with the program to create the summerized haplotypes',
                        default=None,
                        required=False)

    parser.add_argument('--diag_to_ipd_json',
                        type=str,
                        help='diag to ipd look up json created before',
                        default=None,
                        required=False)


def get_process_arrays_args():
    """	Inputs arguments from bash
    Gets the arguments, checks requirements, returns a dictionary of arguments
    Return: args - Arguments as a dictionary
    """
    parser = ArgumentParser()
    parse_process_arrays_args(parser)
    args = parser.parse_args()
    return args


args = get_process_arrays_args()
# same arguments to a local variable by same name as the argument
out_dir = args.out_dir
project_name = args.project_name
config_dir = args.config_dir

bait_fasta = args.bait_fasta
haplotype_lookup = args.haplotype_lookup
animal_lookup_path = args.animal_lookup_path
diag_to_ipd_json = args.diag_to_ipd_json

if bait_fasta is None:
    bait_fasta = os.path.join(config_dir, 'bait.fasta')

if haplotype_lookup is None:
    haplotype_lookup = os.path.join(config_dir, 'haplotype_lookup.csv')
if diag_to_ipd_json is None:
    diag_to_ipd_json = os.path.join(config_dir, 'diag_to_ipd_lookup.json')
#


os.makedirs(out_dir, exist_ok=True)


df_norm_median = pd.read_csv(os.path.join(out_dir, '{0}_norm_median.csv'.format(project_name)))
df_norm_median_miseq = pd.read_csv(os.path.join(out_dir, '{0}_norm_median_miseq.csv'.format(project_name)))
df_read_ct = pd.read_csv(os.path.join(out_dir, '{0}_read_ct.csv'.format(project_name)))

xlsx_filepath = None
os.makedirs(out_dir, exist_ok=True)


def error_code(x, y, z):
    if y > 2:
        return 'TMH:({0})'.format(z)
    if y == 2:
        return x
    if y == 1:
        # add a dash for a second option, as it needs to be looked at to some degree
        return [x, '-']
    return "NO HAPLO"


def highlight_cells(k, x):
    z = list(x)
    z = z[18:]

    def style_cell(y):
        if isinstance(y, str):
            return None
        if math.isnan(y):
            return None
        attribute_list = []
        # ambig = int(((y*100 - math.floor(y*100))*100) + .3)
        Y_str = '{:.4f}'.format(y)

        unique_map = int(Y_str.split('.')[1][0:2])
        ambig = int(Y_str.split('.')[1][2:4])
        # print(unique_map)
        if unique_map < 1.01:
            attribute_list.append("font-weight: bold")
            # attribute_list.append("comment: None")
        elif unique_map > 98:
            if ambig >= 1.5:
                attribute_list.append("color: red")
            else:
                attribute_list.append("font-weight: bold")
                attribute_list.append("color: #006400")
            return ';'.join(attribute_list)
        else:
            attribute_list.append("color: #9900FF")
        if ambig >= 1.5:
            attribute_list.append('background-color: #F2DFA9')
        if len(attribute_list) > 0:
            return ';'.join(attribute_list)
        return None

    return [None] * 18 + [style_cell(y) for y in z]
# this is a bit complicated:
# First, we need unique names for columns that will/will not be merged

print(diag_to_ipd_json)
with open(diag_to_ipd_json) as f_in:
    diag_to_ipd_dict = json.load(f_in)

df = df_norm_median.copy()
df.rename(columns={'ALLELE': 'allele', 'SAMPLE_NUM': 'accession', 'DEPTH_ADJ': 'read_ct'}, inplace=True)

df_miseq = df_norm_median_miseq.copy()
df_miseq.rename(columns={'ALLELE': 'miseq_allele', 'SAMPLE_NUM': 'accession',
                         'DEPTH_ADJ': 'read_ct',
                         'unique_maps_per_allele': 'unique_maps_per_allele_miseq'}, inplace=True)

df_dict = pd.DataFrame()
for key, value in diag_to_ipd_dict.items():
    df_temp = pd.DataFrame({'miseq_allele': key, 'allele': value})
    df_dict = pd.concat([df_dict, df_temp], ignore_index=True)
df_miseq_m = df_miseq.merge(df_dict, on=['miseq_allele'], how='inner')
df_miseq_m['accession'] = df_miseq_m['accession'].astype(str)

df_j = df[['allele', 'accession', 'unique_maps_per_allele']]
df_j['accession'] = df_j['accession'].astype(str)
df_miseq_m_2 = df_miseq_m.merge(df_j, on=['allele', 'accession'], how='left')
df_miseq_m_2 = df_miseq_m_2.fillna(0)
df_miseq = df_miseq_m_2.groupby(['miseq_allele', 'read_ct', 'accession', 'unique_maps_per_allele_miseq'])[
    'unique_maps_per_allele'].max().reset_index().rename(columns={'miseq_allele': 'allele'})
df_miseq = df_miseq[df_miseq['unique_maps_per_allele'] == 0]
df_miseq.drop(columns=['unique_maps_per_allele'], inplace=True)
df_miseq.rename(columns={'unique_maps_per_allele_miseq': 'unique_maps_per_allele'}, inplace=True)
df_miseq.sort_values(by=['accession', 'allele'], inplace=True)

df = pd.concat([df_miseq, df], ignore_index=True)

df_read_ct_pivot = df_read_ct.pivot_table(values='sample_read_ct',
                                          columns=['gs_id'],
                                          aggfunc=np.max,
                                          fill_value=0).reset_index()
df_read_ct_pivot.rename_axis(['index2'], inplace=True, axis=1)
df_read_ct_pivot.rename(columns={'index': 'gs_id'}, inplace=True)
df_read_ct_pivot.rename_axis(['index'], inplace=True, axis=1)
df_median_counts = df.groupby(['accession'])['read_ct'].median().reset_index().rename(
    columns={'read_ct': 'normalized_median_allele_count'})
df_median_counts['normalized_median_allele_count'] = df_median_counts['normalized_median_allele_count'].round(0)

df_median_counts_pivot = df_median_counts.pivot_table(values='normalized_median_allele_count',
                                                      columns=['accession'],
                                                      aggfunc=np.max,
                                                      fill_value=0).reset_index()
df_median_counts_pivot.rename_axis(['index2'], inplace=True, axis=1)
df_median_counts_pivot.rename(columns={'index': 'gs_id'}, inplace=True)
df_median_counts_pivot.rename_axis(['index'], inplace=True, axis=1)

df.rename(columns={'accession': 'gs_id'}, inplace=True)
df['gs_id'] = df['gs_id'].astype(str)

df_haplo_ready = df.copy()
df_haplo_ready.rename(columns={'allele_ipd': 'allele'}, inplace=True)
df_haplo_ready = df_haplo_ready[['allele', 'read_ct', 'gs_id']]
df_haplo_ready.to_csv('/Volumes/T7/ipd/{0}_genotypes.csv'.format(project_name), index=False)
df_read_ct['gs_id'] = df_read_ct['gs_id'].astype(str)

if os.path.exists(animal_lookup_path):
    animal_lookup = pd.read_csv(animal_lookup_path)
    # print(animal_lookup)
    animal_lookup = animal_lookup[['gs_id']]
    animal_lookup_dict = {}
    animal_lookup['gs_id'] = animal_lookup['gs_id'].astype(str)
    df = df.merge(animal_lookup[['gs_id']], how='inner', on=['gs_id'])
    df_read_ct = df_read_ct.merge(animal_lookup[['gs_id']], how='inner', on=['gs_id'])
# df['allele'] = ['-'.join(x.split('-')[1:])  if else for x in df['allele']]


df['allele_ipd'] = ['-'.join(x.split('-')[1:]) if len(x.split('-')) > 1 else x for x in df['allele']]
df['MHC_TYPE'] = [x.split('_')[1] if int((y * 100) % 100) in [49, 99] else x for x, y in
                  zip(df['allele'], df['read_ct'])]

# print(df)
# df['allele_ipd'] = ['-'.join(x.split('-')[1:])  for x in df['allele']]
# df['MHC_TYPE'] = [x.split('_')[1] if int((y * 100) % 100) in [49, 99] else x.split('_')[0] for x, y in
#                   zip(df['allele'], df['read_ct'])]

# print(df)
type_dict = {'Mamu-A1': 'MHC_A_HAPLOTYPES',
             'Mamu-A2': 'MHC_A_HAPLOTYPES',
             'Mamu-B': 'MHC_B_HAPLOTYPES',
             'Mamu-DRB': 'MHC_DRB_HAPLOTYPES',
             'Mamu-DRB1': 'MHC_DRB_HAPLOTYPES',
             'Mamu-DRB3': 'MHC_DRB_HAPLOTYPES',
             'Mamu-DRB4': 'MHC_DRB_HAPLOTYPES',
             'Mamu-DQA1': 'MHC_DQA_HAPLOTYPES',
             'Mamu-DQB1': 'MHC_DQB_HAPLOTYPES',
             'Mamu-DPA1': 'MHC_DPA_HAPLOTYPES',
             'Mamu-DPB1': 'MHC_DPB_HAPLOTYPES'
             }
#  'Mamu-A3':'MHC_A_HAPLOTYPES', #'Mamu-A4', 'Mamu-A6', 'Mamu-AG1',
# 'Mamu-AG2', 'Mamu-AG3', 'Mamu-AG4', 'Mamu-AG5', 'Mamu-AG6',
# 'Mamu-B02Ps', 'Mamu-B10Ps', 'Mamu-B11L', 'Mamu-B14Ps', 'Mamu-B17',
# 'Mamu-B21Ps',
# 'Mamu-E', 'Mamu-G', 'Mamu-I', 'Mamu-J'}
print(list(df['MHC_TYPE'].unique()))
df['TYPE'] = [type_dict[x] if x in type_dict.keys() else x for x in df['MHC_TYPE']]

df.to_csv(os.path.join(out_dir, 'concat_haplotype.csv'), index=False)

# if xlsx_filepath is None:
xlsx_filepath = os.path.join(out_dir, '{0}.pivot.xlsx'.format(project_name))
df['# Obs'] = df.groupby('allele_ipd')['read_ct'].transform('count')
# change to get normalized values
# max_value = df_read_ct['sample_read_ct'].max()
# df = df.merge(df_read_ct, on='gs_id', how='left')
# df['norm_read_ct'] = df['read_ct'] / df['sample_read_ct'] * max_value
# df['norm_read_ct'] = df['norm_read_ct'].round()
# df['read_ct'] = [x - 0.01 if int((y * 100) % 100) in [49, 99] else x for x, y in zip(df['norm_read_ct'],
# df['read_ct'])]
gs_id_list = list(df['gs_id'].unique())
gs_id_list.sort()
df_gs_id = pd.DataFrame([gs_id_list], columns=gs_id_list)

genotype_pivot = df.pivot_table(values='read_ct',
                                index=['allele_ipd', '# Obs'],
                                columns=['gs_id'],
                                aggfunc=np.max,
                                fill_value='').reset_index()
genotype_pivot.rename_axis(['index'], inplace=True, axis=1)
genotype_pivot.rename(columns={'allele_ipd': 'gs_id'}, inplace=True)
print('df_read_ct_pivot')
print(list(genotype_pivot['gs_id']))

df_allele_count_summary = df.groupby('gs_id')['read_ct'].count().reset_index()
df_allele_count_summary.rename(columns={'read_ct': '# Alleles Identified'}, inplace=True)
df_count_summary_pivot = df_allele_count_summary.pivot_table(values='# Alleles Identified',
                                                             columns=['gs_id'],
                                                             aggfunc=np.max,
                                                             fill_value=0).reset_index()
df_count_summary_pivot.rename_axis(['index2'], inplace=True, axis=1)
df_count_summary_pivot.rename(columns={'index': 'gs_id'}, inplace=True)
df_count_summary_pivot.rename_axis(['index'], inplace=True, axis=1)

diag_dict = {}
for diag, ipd_list in diag_to_ipd_dict.items():
    ipd_list = ['-'.join(x.split('-')[1:]) for x in ipd_list]
    diag = '-'.join(diag.split('-')[1:])
    diag_dict[diag] = ','.join(ipd_list)
# print(diag_dict)
genotype_pivot['ambiguous_alleles'] = [x if x not in diag_dict.keys() else diag_dict[x] for x in
                                       genotype_pivot['gs_id']]
genotype_pivot['header_name'] = genotype_pivot['gs_id']
genotype_pivot['gs_id'] = genotype_pivot['gs_id'].str.replace('__', '*')
genotype_pivot['gs_id'] = genotype_pivot['gs_id'].str.replace('_', ':')
genotype_pivot['gs_id'] = genotype_pivot['gs_id'].str.replace('-diag', '-miseq')
genotype_pivot['gs_id'] = genotype_pivot['gs_id'].str.replace('-nuc', '-cdna')
genotype_pivot['gs_id'] = genotype_pivot['gs_id'].str.replace('A1:028g1-miseq', 'AG3:02g1_A028-miseq')
genotype_pivot['gs_id'] = genotype_pivot['gs_id'].str.replace('A1:110g1-miseq', 'A1:110g1-cdna-miseq')
genotype_pivot['gs_id'] = genotype_pivot['gs_id'].str.replace('A1:119g1-miseq', 'AG3:02g2_A119-miseq')
genotype_pivot['gs_id'] = [x.replace(':', '*', 1) if x.endswith('-miseq') else x for x in genotype_pivot['gs_id']]
genotype_pivot['gs_id'] = genotype_pivot['gs_id'].str.replace('-gen', '')
# print(genotype_pivot)
# diag_dict
# df_count_summary_pivot.merge()

df_list = [df_gs_id, df_count_summary_pivot, df_read_ct_pivot, df_median_counts_pivot, genotype_pivot]
# df_merge = df.copy()
df_haplo_pivot = pd.DataFrame({'gs_id': ['MHC A HAPLOTYPES 1', 'MHC A HAPLOTYPES 2',
                                         'MHC B HAPLOTYPES 1', 'MHC B HAPLOTYPES 2',
                                         'MHC DRB HAPLOTYPES 1', 'MHC DRB HAPLOTYPES 2',
                                         'MHC DQA HAPLOTYPES 1', 'MHC DQA HAPLOTYPES 2',
                                         'MHC DQB HAPLOTYPES 1', 'MHC DQB HAPLOTYPES 2',
                                         'MHC DPA HAPLOTYPES 1', 'MHC DPA HAPLOTYPES 2',
                                         'MHC DPB HAPLOTYPES 1', 'MHC DPB HAPLOTYPES 2']})
if os.path.exists(haplotype_lookup):
    def error_code(x, y, z):
        if y > 2:
            return 'TMH:({0})'.format(z)
        if y == 2:
            return x
        if y == 1:
            # add a dash for a second option, as it needs to be looked at to some degree
            return [x, '-']
        return "NO HAPLO"


    df_haplo_call = pd.read_csv(haplotype_lookup)
    print(df_haplo_call)
    df_haplo_ready.rename(columns={'allele': 'allele_ipd'}, inplace=True)
    df_haplo_merge = df_haplo_call.merge(df_haplo_ready[['allele_ipd', 'gs_id']], on=['allele_ipd'], how='inner')
    df_haplo_merge = df_haplo_merge[
        ['HAPLOTYPE_CALL', 'PREFIX', 'TYPE', 'allele', 'gs_id', 'CALL_COUNT']].drop_duplicates()
    df_haplo_merge['CALL_GSID_COUNT'] = df_haplo_merge.groupby(['HAPLOTYPE_CALL', 'PREFIX', 'TYPE', 'gs_id'])[
        'allele'].transform('count')
    df_haplo_pass = df_haplo_merge[df_haplo_merge['CALL_COUNT'] == df_haplo_merge['CALL_GSID_COUNT']]
    df_haplo_pass = df_haplo_pass[['HAPLOTYPE_CALL', 'PREFIX', 'TYPE', 'gs_id', 'CALL_COUNT']].drop_duplicates()
    df_haplo_pass['PASS_COUNT'] = df_haplo_pass.groupby(['gs_id', 'TYPE'])['HAPLOTYPE_CALL'].transform('count')
    df_haplo_pass['rank'] = df_haplo_pass.groupby(['gs_id', 'TYPE'])['CALL_COUNT'].transform('rank', method='first')
    df_haplo_pass['TMH'] = df_haplo_pass.groupby(['gs_id', 'TYPE'])['HAPLOTYPE_CALL'].transform('; '.join)
    df_haplo_pass['HAPLOTYPE_FINAL'] = [error_code(x, y, z) for x, y, z in
                                        zip(df_haplo_pass['HAPLOTYPE_CALL'], df_haplo_pass['PASS_COUNT'],
                                            df_haplo_pass['TMH'])]
    # if it has only one haplotype, it will make a list [haplotype_call, '-'] which then can be xploded
    df_haplo_pass = df_haplo_pass.explode('HAPLOTYPE_FINAL')
    df_haplo_pass['rank'] = df_haplo_pass.groupby(['gs_id', 'TYPE'])['CALL_COUNT'].transform('rank', method='first')
    df_haplo_pass['TYPE'] = ['{0} {1}'.format(x.replace('_', ' '), int(y)) for x, y in
                             zip(df_haplo_pass['TYPE'], df_haplo_pass['rank'])]
    df_haplo_pass = df_haplo_pass[df_haplo_pass['rank'] < 3]
    df_haplo_pivot = df_haplo_pass.pivot_table(values='HAPLOTYPE_FINAL',
                                               index=['TYPE'],
                                               columns=['gs_id'],
                                               aggfunc=np.max,
                                               fill_value='NO HAPLO').reset_index()
    df_haplo_pivot.rename_axis(['index'], inplace=True, axis=1)
    df_haplo_pivot.rename(columns={'TYPE': 'gs_id'}, inplace=True)

df_xlx_pivot = pd.concat(
    [df_gs_id, df_count_summary_pivot, df_haplo_pivot, df_read_ct_pivot, df_median_counts_pivot, genotype_pivot],
    ignore_index=True)

df_xlx_pivot.rename(columns={'gs_id': 'Animal IDs'}, inplace=True)
if (len(animal_lookup_path) > 0) and os.path.exists(animal_lookup_path):
    animal_lookup = pd.read_csv(animal_lookup_path)
    animal_lookup_dict = {}
    animal_lookup['gs_id'] = animal_lookup['gs_id'].astype(str)
    for idx, row in animal_lookup.iterrows():
        animal_lookup_dict[row['gs_id']] = row['animal_id']
    gs_id_list = list(df_xlx_pivot.columns)
    df_xlx_pivot.rename(columns=animal_lookup_dict, inplace=True)
    animal_list = list(animal_lookup['animal_id'])
    animal_list2 = []

    for gs_id_i in gs_id_list:
        if gs_id_i in animal_lookup_dict.keys():
            animal_list2.append(animal_lookup_dict[gs_id_i])
    animal_list2.sort()
    col_order = ['Animal IDs', '# Obs'] + animal_list2 + ['header_name', 'ambiguous_alleles']
else:
    col_names = list(df_xlx_pivot.columns)
    col_names.remove('Animal IDs')
    col_names.remove('# Obs')
    col_names.remove('ambiguous_alleles')
    col_order = ['Animal IDs', '# Obs'] + col_names + ['header_name', 'ambiguous_alleles']
    animal_list2 = col_names

df_xlx_pivot = df_xlx_pivot[col_order]

df_xlx_pivot_2 = df_xlx_pivot.copy()
for animal_i in animal_list2:
    df_xlx_pivot_2[animal_i] = [int(x) if isinstance(x, float) and not math.isnan(x) else x for x in
                                df_xlx_pivot[animal_i]]

first_item = True
for animal_i in animal_list2:
    if first_item:
        df_xlx_pivot_2 = df_xlx_pivot_2.style.apply(highlight_cells, x=list(df_xlx_pivot[animal_i]), subset=animal_i)
        first_item = False
        continue
    df_xlx_pivot_2 = df_xlx_pivot_2.apply(highlight_cells, x=list(df_xlx_pivot[animal_i]), subset=animal_i)
# for animal_i in animal_list2:
#     df_xlx_pivot[animal_i] = [int(x) if isinstance(x, float) else x for x in df_xlx_pivot[animal_i]]


df_xlx_pivot_2.to_excel(xlsx_filepath, sheet_name='Sheetname_1', index=False)
print(xlsx_filepath)
