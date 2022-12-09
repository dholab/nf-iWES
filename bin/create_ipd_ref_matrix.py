from argparse import ArgumentParser, ArgumentTypeError
import pandas as pd
import numpy as np
import os
from Bio import SeqIO


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
    parser.add_argument('--bait_fasta',
                        type=str,
                        help='Where you ipd-miseq database combined fasta exists',
                        default=None,
                        required=False)
    parser.add_argument('--config_dir',
                        type=str,
                        help='where the config files are stored',
                        default=None,
                        required=False)

    parser.add_argument('--ipd_ref_matrix_dir',
                        type=str,
                        help='Directory Where your ipd reference matrix files are stored',
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
ipd_ref_matrix_dir = args.ipd_ref_matrix_dir
bait_fasta = args.bait_fasta
config_dir = args.config_dir
if (config_dir is None) and ((bait_fasta is None) or (ipd_ref_matrix_dir is None)):
    print("A config_dir or  bait_fast and ipd_ref_matrix_dir must be declared")
    exit()
if ipd_ref_matrix_dir is None:
    ipd_ref_matrix_dir = os.path.join(config_dir, 'ipd_ref_matrix')
if bait_fasta is None:
    bait_fasta = os.path.join(config_dir, 'bait.fasta')


def fasta_to_df(fasta_path=None, header_name='ALLELE', sequence_name='SEQUENCE', as_df=False):
    fasta_sequences = SeqIO.parse(open(fasta_path), 'fasta')
    fasta_dict = {}
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        fasta_dict[name] = sequence
    fasta_sequences.close()
    if as_df:
        return pd.DataFrame(fasta_dict.items(), columns=[header_name, sequence_name])
    return fasta_dict


def make_lists(seq):
    seq_l = len(seq)
    seq_list = []
    prev_seq = ''
    for x in range(-75, seq_l - 75):
        y = x + 151
        if x < 0:
            x = 0
        if y > seq_l:
            y = seq_l
        # if seq[x:y] != prev_seq:
        seq_list.append(seq[x:y])
        # prev_seq = seq[x:y]
    return seq_list


def compare_semi(df_exp_compare, df_exp, start=True):
    def create_range_for(y, z, k):
        if k + 74 > z:
            return list(range(y, z - k + y + 1))
        return list(range(y, 74 + y + 1))

    def filter_dyn_for(x, y):
        return x[-y:]

    def create_range_rev(y, z, k):
        if k - 74 < 0:
            return list(range(y - k, y + 1))
        return list(range(y - 74, y + 1))

    def filter_dyn_rev(x, y):
        return x[:y]

    filter_dyn = {True: filter_dyn_for, False: filter_dyn_rev}
    create_range = {True: create_range_for, False: create_range_rev}

    df_exp_short = df_exp[(df_exp['LEN'] < df_exp['REF_SEQ_LEN']) & (df_exp['LEN'] < 151)]
    if start:
        df_exp_short_start = df_exp_short[(df_exp_short['START'] < 1)]
        filter_exact = (-76, None)
    else:
        df_exp_short_start = df_exp_short[(df_exp_short['START'] > 0)]
        filter_exact = (None, 76)

        # Lets Find trim in half no sense checking any thing that is not a partial
    df_short_seg_start = df_exp_short_start[df_exp_short_start['LEN'] == 76]
    ### different
    df_exp_compare['SEQ_LIST'] = [x[filter_exact[0]:filter_exact[1]] for x in df_exp_compare['SEQ_LIST_2']]
    df_short_seg_start_merge = df_short_seg_start.merge(df_exp_compare, on=['SEQ_LIST'], how='inner')
    df_start_merge_same = df_short_seg_start_merge[
        df_short_seg_start_merge['ALLELE'] == df_short_seg_start_merge['ALLELE_2']]

    df_start_merge_diff = df_short_seg_start_merge[
        df_short_seg_start_merge['ALLELE'] != df_short_seg_start_merge['ALLELE_2']]
    df_start_merge_diff = df_start_merge_diff[df_start_merge_diff['LEN'] != df_start_merge_diff['LEN_2']]
    keep_alleles = list(df_start_merge_diff['ALLELE_2'].unique())
    if start:
        df_start_merge_diff['index_2'] = [create_range[start](y, z, k) for y, z, k in
                                          zip(df_start_merge_diff['index_2'],
                                              df_start_merge_diff['REF_SEQ_LEN_2'],
                                              df_start_merge_diff['END_2'])]
        df_start_merge_diff['index_1'] = [create_range[start](y, z, k) for y, z, k in
                                          zip(df_start_merge_diff['index_1'],
                                              df_start_merge_diff['REF_SEQ_LEN_2'],
                                              df_start_merge_diff['END_2'])]
    else:
        df_start_merge_diff['index_2'] = [create_range[start](y, z, k) for y, z, k in
                                          zip(df_start_merge_diff['index_2'],
                                              df_start_merge_diff['REF_SEQ_LEN_2'],
                                              df_start_merge_diff['START_2'])]
        df_start_merge_diff['index_1'] = [create_range[start](y, z, k) for y, z, k in
                                          zip(df_start_merge_diff['index_1'],
                                              df_start_merge_diff['REF_SEQ_LEN_2'],
                                              df_start_merge_diff['START_2'])]

    # print(df_start_merge_diff[df_start_merge_diff['ALLELE']=='Mamu-B_055g1-diag'])

    df_start_merge_index1 = df_start_merge_diff[['index_1']]
    df_start_merge_index2 = df_start_merge_diff[['index_2']]
    df_start_merge_index1 = df_start_merge_index1.explode('index_1')
    df_start_merge_index2 = df_start_merge_index2.explode('index_2')

    # Reset the index
    df_start_merge_index2 = df_start_merge_index2.reset_index(drop=True).reset_index()
    df_start_merge_index1 = df_start_merge_index1.reset_index(drop=True).reset_index()
    df_start_merge_index2_all = df_start_merge_index2.merge(df_exp_compare, on=['index_2'], how='left')
    df_start_merge_index2_all.sort_values(by=['index'], inplace=True)
    df_start_merge_index1_all = df_start_merge_index1.merge(df_exp, on=['index_1'], how='left')
    df_start_merge_index1_all.sort_values(by=['index'], inplace=True)
    df_start_merge_index2_all.drop(columns=['SEQ_LIST'], inplace=True)
    df_start_all = df_start_merge_index1_all.merge(df_start_merge_index2_all, on='index', how='inner')
    #
    df_start_all.dropna(subset=['LEN', 'LEN_2'], inplace=True)
    # df_start_all.dropna(subset=['LEN_2'], inplace=True)
    df_start_all['LEN'] = df_start_all['LEN'].astype(int)
    # print(df_start_all[df_start_all['SEQ_LIST_2'].isna()])
    df_start_all['SEQ_2_TRIM'] = [filter_dyn[start](x, y) for x, y in
                                  zip(df_start_all['SEQ_LIST_2'], df_start_all['LEN'])]
    # df_start_all['SEQ_2_TRIM_LEN'] = [len(x) for x in df_start_all['SEQ_LIST']]
    df_start_all_filtered = df_start_all[df_start_all['SEQ_2_TRIM'] == df_start_all['SEQ_LIST']]
    return df_start_all_filtered


def compare_shorter(df_exp_compare, df_exp):
    df_exp_short_ref_seq = df_exp[(df_exp['REF_SEQ_LEN'] < 149.5) & (df_exp['LEN'] == df_exp['REF_SEQ_LEN'])]
    df_short_seq_len = pd.DataFrame()
    for idx, row in df_exp_short_ref_seq.iterrows():
        print(row['ALLELE'])
        df_exp_i = df_exp_compare[df_exp_compare['SEQ_LIST_2'].str.contains(row['SEQ_LIST'])]
        index_list = [idx] * len(df_exp_i)
        df_c = pd.DataFrame(row.to_dict(), index=index_list)
        df_concat = pd.concat([df_c.reset_index(drop=True), df_exp_i.reset_index(drop=True)], axis=1,
                              ignore_index=False)
        df_short_seq_len = pd.concat([df_short_seq_len, df_concat], ignore_index=True)
    return df_short_seq_len


def compare_full(df_exp_compare, df_exp, out_dir, file_count=10):
    # shoten the columns to save space and change type from float to int
    df_exp_compare_full = df_exp_compare[['ALLELE_2', 'SEQ_LIST_2', 'START_2', 'END_2']]
    df_exp_compare_full.rename(columns={'SEQ_LIST_2': 'SEQ_LIST'}, inplace=True)
    df_exp_compare_full['START_2', 'END_2'] = df_exp_compare_full['START_2', 'END_2'].astype(int)

    df_exp_full = df_exp[['ALLELE', 'SEQ_LIST', 'START', 'END']]
    df_exp_full['START', 'END'] = df_exp_full['START', 'END'].astype(int)

    for idx, chunk in enumerate(np.array_split(df_exp_compare_full, file_count)):
        print(idx, file_count)
        df_exp_full_merge = df_exp_full.merge(chunk, on=['SEQ_LIST'], how='inner')
        # df_exp_full_merge = df_exp_full_merge[df_exp_full_merge['allele'] != df_exp_full_merge['allele_2']]
        df_exp_full_merge.drop(columns=['SEQ_LIST'], inplace=True)
        print('start gzip')
        df_exp_full_merge.to_csv(os.path.join(out_dir, '{0}_{1}_ipd_full_matrix.csv.gz'.format(idx, file_count)),
                                 compression='gzip',
                                 index=False)

    del df_exp_full_merge
    del chunk
    del df_exp_compare_full
    del df_exp_full
    return


os.makedirs(ipd_ref_matrix_dir, exist_ok=True)
# Read the fasta as a dataframe
df = fasta_to_df(fasta_path=bait_fasta, header_name='ALLELE', sequence_name='SEQUENCE', as_df=True)
df['SEQ_LIST'] = [make_lists(x) for x in df['SEQUENCE']]
# create a sequence of every possibile segment of the reference sequence from 76 to 151 in length
df_exp = df.explode('SEQ_LIST')
df_exp = df_exp.reset_index().rename(columns={'index': 'allele_num'})
df_exp['rank'] = df_exp.groupby('ALLELE')['allele_num'].rank('first')
df_exp['rank'] = df_exp['rank'].astype(int)
df_exp['LEN'] = [len(x) for x in df_exp['SEQ_LIST']]
df_exp['REF_SEQ_LEN'] = [len(x) for x in df_exp['SEQUENCE']]
# needs to consider ref length
df_exp['START'] = [0 if x < 76 else x - 76 for x in df_exp['rank']]
df_exp['END'] = [x + y for x, y in zip(df_exp['START'], df_exp['LEN'])]
df_exp = df_exp[['ALLELE',
                 'REF_SEQ_LEN',
                 'SEQ_LIST',
                 'START',
                 'END',
                 'LEN']].drop_duplicates(keep='first')

df_exp = df_exp.reset_index()
print('create a copy of the fasta converted and expanded dataframe to compare')
df_exp_compare = df_exp.copy()
# make a dataframe that we will join to later for comparisons
# df_exp_compare = df_exp[['index','allele','SEQ_LIST','START','END','LEN','REF_SEQ_LEN']]
df_exp_compare.rename(columns={'index': 'index_2',
                               'ALLELE': 'ALLELE_2',
                               'REF_SEQ_LEN': 'REF_SEQ_LEN_2',
                               'SEQ_LIST': 'SEQ_LIST_2',
                               'START': 'START_2',
                               'END': 'END_2',
                               'LEN': 'LEN_2'}, inplace=True)
# df_exp_compare = df_exp_compare.reset_index().rename(columns={'index':'index_2'})
df_exp.rename(columns={'index': 'index_1'}, inplace=True)

print('create a start based semi compared matrix')
df_start_all_filtered = compare_semi(df_exp_compare, df_exp, start=True)

print('create an end based semi compared matrix')
df_end_all_filtered = compare_semi(df_exp_compare, df_exp, start=False)
print('create shorter than reference sequences lookup')
df_short_seq_len = compare_shorter(df_exp_compare, df_exp)
print('concat shorter comparison dataframes')
df_ipd_matrix = pd.concat([df_end_all_filtered[['ALLELE', 'START', 'END', 'ALLELE_2', 'START_2', 'END_2']],
                           df_start_all_filtered[['ALLELE', 'START', 'END', 'ALLELE_2', 'START_2', 'END_2']],
                           df_short_seq_len[['ALLELE', 'START', 'END', 'ALLELE_2', 'START_2', 'END_2']]],
                          ignore_index=True)
print('export dataframe to gzip')
df_ipd_matrix.drop_duplicates(inplace=True)
df_ipd_matrix.to_csv(os.path.join(ipd_ref_matrix_dir, 'ipd_short_matrix.csv.gz'), index=False, compression='gzip')
print('create compare full lengths and export to files .gz')
compare_full(df_exp_compare, df_exp, ipd_ref_matrix_dir, file_count=10)
print('collecting excess variables')
del df
# n = gc.collect()
print("Number of unreachable objects collected by GC:")
