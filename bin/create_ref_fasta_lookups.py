#!/usr/bin/env python3
from argparse import ArgumentParser, ArgumentTypeError
import pandas as pd
import numpy as np
import os
from Bio import SeqIO
import json
import pysam


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
    parser.add_argument('--config_dir',
                        type=str,
                        help='where the config files are stored',
                        required=True)
    parser.add_argument('--miseq_legacy_db_path',
                        type=str,
                        help='the miseq legacy database fasta filepath',
                        required=True)

    parser.add_argument('--gen_db_path',
                        type=str,
                        help='the genomic database fasta filepath',
                        required=True)
    parser.add_argument('--exon_db_path',
                        type=str,
                        help='the exon2 database fasta filepath',
                        required=True)

    parser.add_argument('--haplotype_json_path',
                        type=str,
                        help='haplotype_json_path',
                        required=True)
    parser.add_argument('--species',
                        type=str,
                        help='species',
                        default='Mamu',
                        required=False)

    parser.add_argument('--cp_path',
                        type=str,
                        help='Directory Where your bbmap executables are stored',
                        required=True)
    parser.add_argument('--threads',
                        type=int,
                        help='Number of threads to run bbmap',
                        default=1,
                        required=False)
    parser.add_argument('--ram',
                        type=int,
                        help='Directory ram to dedicate to bbmap java',
                        default=8000,
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
# cp_path = '~/anaconda3/bin/current'
config_dir = args.config_dir
miseq_legacy_db_path = args.miseq_legacy_db_path
exon_db_path = args.exon_db_path
gen_db_path = args.gen_db_path
haplotype_json_path = args.haplotype_json_path
cp_path = args.cp_path
threads = args.threads
ram = args.ram
species = args.species
raw_fasta = {'gen': gen_db_path, 'exon': exon_db_path}

# change to where your bbmap is to launch java version. (add '/current' to the directory bbmap executables are located)
# what delimiter to split the header names to shorten them: '' is no split
split_char_ipd = ' '
split_char_diag = ''

os.makedirs(config_dir, exist_ok=True)

# ipd_fasta_dir = os.path.join(config_dir, 'ipd')
# diag_fasta_dir = os.path.join(config_dir, 'miseq')
all_filepath = os.path.join(config_dir, 'bait.fasta')
ipd_all_filepath = os.path.join(config_dir, 'ipd.fasta')
miseq_all_filepath = os.path.join(config_dir, 'miseq.fasta')
miseq_long_all_filepath = os.path.join(config_dir, 'miseq_long.fasta')
miseq_ungrouped_all_filepath = os.path.join(config_dir, 'miseq_ungrouped.fasta')
ipd_to_miseq_lookup_path = os.path.join(config_dir, 'ipd_to_gen_lookup.json')
missing_ipd_alleles_path = os.path.join(config_dir, 'missing_ipd_alleles.csv')
miseq_to_ipd_lookup_path = os.path.join(config_dir, 'miseq_to_ipd_lookup.json')
missing_miseq_alleles_path = os.path.join(config_dir, 'missing_diag_alleles.csv')
bam_inpath = os.path.join(config_dir, 'miseq_to_ipd.sam')
bam_outpath = os.path.join(config_dir, 'miseq_to_ipd_filtered.sam')
bam_u_outpath = os.path.join(config_dir, 'miseq_to_ipd_unmapped.sam')
miseq_legacy_short_path = os.path.join(config_dir, 'miseq_legacy_short.fasta')
ipd_miseq_allele_haplotype_path = os.path.join(config_dir, 'haplotype_lookup.csv')
miseq_haplo_sam_path = os.path.join(config_dir, 'miseq_haplo.sam')
ipd_haplo_sam_path = os.path.join(config_dir, 'ipd_haplo.sam')
ipd_to_exon_lookup_path = os.path.join(config_dir, 'gen_to_exon_lookup.json')
exon_to_ipd_lookup_path = os.path.join(config_dir, 'exon_to_gen_lookup.json')
missing_ipd_exom_alleles_path = os.path.join(config_dir, 'missing_ipd_exon_alleles.csv')
missing_exp_ipd_alleles_path = os.path.join(config_dir, 'missing_exon_ipd_alleles.csv')


def rev_comp_st(seq):
    seq = seq.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
    seq = seq.upper()
    seq = seq[::-1]
    return seq


def trim_read_for_diag(seq, cigartupple):
    i = 0
    trim_start = 0
    trim_end = 0
    final_tupple = len(cigartupple) - 1
    for cigar_i in cigartupple:
        if i == 0:
            if cigar_i[0] == 4:
                trim_start = cigar_i[1]
        if i == final_tupple:
            if cigar_i[0] == 4:
                trim_end = cigar_i[1]
        i += 1

    if trim_end == 0:
        return seq[trim_start:]
    return seq[trim_start:-trim_end]


def write_sequence(name, seq, out_file, replace_illegal_char=True, prefix='', suffix=''):
    if replace_illegal_char:
        name = name.replace('*', '__')  # add a double because i want to fine the 02 at the beginning easir
        name = name.replace(':', '_')
    out_file.write(
        '>{0}{1}{2}\n'.format(prefix, name, suffix))  # . We do mostly non human primates I want to differentiate.
    out_file.write('{0}\n'.format(seq))


def fasta_to_df(fasta_path=None, header_name='allele', sequence_name='SEQUENCE', as_df=False):
    fasta_sequences = SeqIO.parse(open(fasta_path), 'fasta')
    fasta_dict = {}
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        fasta_dict[name] = sequence
    fasta_sequences.close()
    if as_df:
        return pd.DataFrame(fasta_dict.items(), columns=[header_name, sequence_name])
    return fasta_dict


def single_file_per_fasta(fasta_path, single_fasta_dir, name_num_json_path):
    i = 0
    ipd_num = {}
    fasta_sequences = SeqIO.parse(open(fasta_path), 'fasta')
    for fasta in fasta_sequences:
        i += 1
        name, sequence = fasta.id, str(fasta.seq)
        if len(sequence) > 1:
            truncated_ref = os.path.join(single_fasta_dir, '{0}.fasta'.format(i))
            ipd_num[name] = '{0}'.format(i)
            with open(truncated_ref, 'w') as out_file:
                out_file.write('>{0}\n'.format(name))
                out_file.write('{0}\n'.format(sequence))
    with open(name_num_json_path, 'w') as convert_file:
        convert_file.write(json.dumps(ipd_num))
    return


def like_join(x, df_ref, column_i='SEQUENCE'):
    name_list = []
    for idx, row in df_ref.iterrows():
        if x in row['allele']:
            name_list.append(row[column_i])
    return name_list


def rev_comp_like_join(inner_seq, df_ref):
    name_list = []
    for idx, row in df_ref.iterrows():
        if (inner_seq in row['SEQUENCE']) or (rev_comp_st(inner_seq) in row['SEQUENCE']):
            name_list.append(row['allele'])
    return name_list


#####################################
#   Generate legal name for fasta   #
#####################################
# name_lengths = 0
# name_long = ''
print('generate legal names for the fasta files (under 40 characters, no special characters')
fasta_dict = {}
ipd_name_lookup = {}
for suffix in ['gen', 'exon']:
    fasta_sequences = SeqIO.parse(open(raw_fasta[suffix]), 'fasta')
    for fasta in fasta_sequences:
        name = fasta.name
        #         if len(name) > name_lengths:
        #             name_lengths = len(name)
        #             name_long = name
        if len(name) > 1:
            fasta_dict[name.split('|')[0]] = str(fasta.seq)
            name_short = '{0}-{1}'.format(name.split('|')[0], suffix)
            ipd_name_lookup[name_short] = '{0}-{1}'.format(name, suffix)
    name_list = list(fasta_dict.keys())
    name_list.sort()
    with open(os.path.join(config_dir, '{0}.fasta'.format(suffix)), 'w') as out_file:

        for name in name_list:
            write_sequence(name=name,
                           seq=fasta_dict[name],
                           out_file=out_file,
                           replace_illegal_char=True,
                           prefix='',
                           suffix='-{0}'.format(suffix))
# concatenate the exon2 and genomic files
with open(ipd_all_filepath, 'w') as out_file:
    for suffix in ['gen', 'exon']:
        fasta_sequences = SeqIO.parse(open(os.path.join(config_dir, '{0}.fasta'.format(suffix))), 'fasta')
        for fasta in fasta_sequences:
            name = fasta.name
            write_sequence(name=name,
                           seq=str(fasta.seq),
                           out_file=out_file,
                           replace_illegal_char=True,
                           prefix='',
                           suffix='')
# create legal names for miseq
with open(miseq_legacy_short_path, 'w') as out_file:
    fasta_sequences = SeqIO.parse(open(miseq_legacy_db_path), 'fasta')
    for fasta in fasta_sequences:
        name = fasta.description
        name = name.split('|')[0]
        if len(name) > 1:
            write_sequence(name=name,
                           seq=str(fasta.seq),
                           out_file=out_file,
                           replace_illegal_char=True,
                           prefix='')
# map the miseq data base to the ipd database.  This is to help trim what is missing.
os.system(f'! minimap2 -N 1 --secondary=no -a {miseq_legacy_short_path} {ipd_all_filepath} > {bam_inpath}')

# create the bam files to remove the unmapped files.
# This also creates lists of interest that are used to make a miseq database fasta for all inclusive purposes.
bf_in = pysam.AlignmentFile(bam_inpath, "r")
bf_out = pysam.AlignmentFile(bam_outpath, "w", template=bf_in)
bf_u_out = pysam.AlignmentFile(bam_u_outpath, "w", template=bf_in)
i = 0
read_list_mapped = []
read_list_unmapped = []
cigar_list_mapped = []
# sequence_list_map = []
trimmed_sequence = []
read_dict = {}

for r in bf_in.fetch(until_eof=True):
    i = i + 1
    mapped_read = str(r).split('\t')
    # if mapped_read[0].startswith('Mamu-A2__'):
    if int(mapped_read[1]) == 4:
        bf_u_out.write(r)
        read_list_unmapped.append(mapped_read[0])
    elif int(mapped_read[1]) == 0:
        read_list_mapped.append(mapped_read[0])
        cigar_list_mapped.append(r.cigartuples)
        # sequence_list_map.append(r.query_sequence)
        trimmed_sequence.append(r.query_alignment_sequence)

        bf_out.write(r)
bf_u_out.close()
bf_out.close()
bf_in.close()

# trimmed_sequence = [trim_read_for_diag(x,y) for x, y in zip(sequence_list_map,cigar_list_mapped)]
read_list_mapped, trimmed_sequence = zip(*sorted(zip(read_list_mapped, trimmed_sequence)))
print(trimmed_sequence[:3])
print(read_list_mapped[:3])

# sort by name
read_list_mapped, trimmed_sequence = zip(*sorted(zip(read_list_mapped, trimmed_sequence)))
with open(miseq_ungrouped_all_filepath, 'w') as out_file:
    for read_name, seq in zip(read_list_mapped, trimmed_sequence):
        write_sequence(name=read_name,
                       seq=seq,
                       out_file=out_file,
                       replace_illegal_char=True,
                       prefix='')

####################################
# deduplicate miseq database fasta #
####################################
print('deduplicate miseq database fasta')
df_diag = fasta_to_df(fasta_path=miseq_ungrouped_all_filepath, header_name='allele', sequence_name='SEQUENCE',
                      as_df=True)
df_diag_dedup = df_diag.drop_duplicates(subset=['SEQUENCE'], keep='first', ignore_index=False)
df_diag_dedup.rename(columns={'allele': 'allele_first'}, inplace=True)
df_diag_merge = df_diag.merge(df_diag_dedup, on=['SEQUENCE'], how='inner')
df_diag_merge['allele'] = [x[5:] for x in df_diag_merge['allele']]
df_diag_merge = (df_diag_merge.groupby(['SEQUENCE', 'allele_first'])
                 .agg({'allele': lambda x: ','.join(x)})
                 .reset_index())
df_diag_merge.sort_values('allele_first', inplace=True)
df_diag_merge = df_diag_merge.reset_index()
df_diag_merge['allele_group_name'] = [x.split('-')[1] for x in df_diag_merge['allele_first']]
df_diag_merge['allele_group_name'] = ['_'.join(x.split('_')[:3]) for x in df_diag_merge['allele_group_name']]
df_diag_merge['allele_group_name'] = [x.replace('__', '_') for x in df_diag_merge['allele_group_name']]
df_diag_merge['group_number'] = df_diag_merge.groupby(['allele_group_name'])['allele_group_name'].cumcount().add(1)
df_diag_merge['header'] = ['Mamu-{0}g{1}|{2}'.format(x, y, z) if z.count(',') > 0 else j for x, y, z, j in
                           zip(df_diag_merge['allele_group_name'],
                               df_diag_merge['group_number'],
                               df_diag_merge['allele'],
                               df_diag_merge['allele_first'])]
# add the previous vetted miseq file (this contains novel and links to the haplotype file(s))
df_diag_2 = fasta_to_df(fasta_path=miseq_legacy_db_path, header_name='header', sequence_name='SEQUENCE', as_df=True)
df_diag = pd.concat([df_diag_merge, df_diag_2], ignore_index=True)
# df_diag['allele'] = [x.split('|')[0] for x in df_diag['allele']]
df_diag_dedup = df_diag.drop_duplicates(subset=['SEQUENCE'], keep='first', ignore_index=False)

df_diag_fasta = df_diag_dedup[['header', 'SEQUENCE']]
df_diag_fasta.sort_values('header', inplace=True)

with open(miseq_long_all_filepath, 'w') as out_file:
    for indx, row in df_diag_fasta.iterrows():
        if len(row['SEQUENCE']) > 130:
            write_sequence(name=row['header'],
                           seq=row['SEQUENCE'],
                           out_file=out_file,
                           replace_illegal_char=True,
                           prefix='')

df_diag_fasta['header'] = ['{0}-miseq'.format(x.split('|')[0]) for x in df_diag_fasta['header']]
with open(miseq_all_filepath, 'w') as out_file:
    for indx, row in df_diag_fasta.iterrows():
        if len(row['SEQUENCE']) > 130:
            write_sequence(name=row['header'].split('|')[0],
                           seq=row['SEQUENCE'],
                           out_file=out_file,
                           replace_illegal_char=True,
                           prefix='')
        else:
            print("Seq shorter than min length:", len(row['SEQUENCE']), row['header'].split('|')[0])

#############################################
# Concatenate the ipd and miseq fasta files #
#############################################
print('Concatenate the ipd and miseq fasta files')
pathlist = [ipd_all_filepath, miseq_all_filepath]
with open(all_filepath, 'w') as out_file:
    print(all_filepath)
    for path_i in pathlist:
        fasta_sequences = SeqIO.parse(open(path_i), 'fasta')
        for fasta in fasta_sequences:
            name = fasta.description
            if len(name) > 1:
                write_sequence(name=name,
                               seq=str(fasta.seq),
                               out_file=out_file,
                               replace_illegal_char=True,
                               prefix='')

############################################################################
# Create IPD to Miseq dict & missing list of IPD sequences w/out Miseq match #
############################################################################
print('Create IPD to Miseq dict')
diag_fasta_dict = fasta_to_df(fasta_path=miseq_all_filepath)
ipd_fasta_dict = fasta_to_df(fasta_path=ipd_all_filepath)

ipd_to_diag_lookup = {}
for ipd_name, ipd_seq in ipd_fasta_dict.items():
    match_list = []
    for diag_name, diag_seq in diag_fasta_dict.items():
        if (diag_seq in ipd_seq) or (rev_comp_st(diag_seq) in ipd_seq):
            match_list.append(diag_name)
    ipd_to_diag_lookup[ipd_name] = match_list

# out put to gson
with open(ipd_to_miseq_lookup_path, 'w') as convert_file:
    convert_file.write(json.dumps(ipd_to_diag_lookup))
# get the missing list
missing_ipd_list = []
for ipd_seq, diag_seq_list in ipd_to_diag_lookup.items():
    if len(diag_seq_list) < 1:
        missing_ipd_list.append(ipd_seq)
missing_ipd_list.sort()
# output do data frame
df_missing_ipd = pd.DataFrame({'allele': missing_ipd_list})
df_missing_ipd.to_csv(missing_ipd_alleles_path, header=False, index=False)

############################################################################
# Create Miseq to IPD dict & missing list of Miseq sequences w/out IPD match #
############################################################################
print('Create Miseq to IPD dict')
diag_fasta_dict = fasta_to_df(fasta_path=miseq_all_filepath)
ipd_fasta_dict = fasta_to_df(fasta_path=ipd_all_filepath)

diag_to_ipd_lookup = {}
for diag_name, diag_seq in diag_fasta_dict.items():
    match_list = []
    for ipd_name, ipd_seq in ipd_fasta_dict.items():
        if (diag_seq in ipd_seq) or (rev_comp_st(diag_seq) in ipd_seq):
            match_list.append(ipd_name)
    diag_to_ipd_lookup[diag_name] = match_list

# out put to gson
with open(miseq_to_ipd_lookup_path, 'w') as convert_file:
    convert_file.write(json.dumps(diag_to_ipd_lookup))
# get the missing list
missing_diag_list = []
for diag_seq, ipd_seq_list in diag_to_ipd_lookup.items():
    if len(ipd_seq_list) < 1:
        missing_diag_list.append(diag_seq)
missing_diag_list.sort()
# output do data frame
df_missing_diag = pd.DataFrame({'allele': missing_diag_list})
df_missing_diag.to_csv(missing_miseq_alleles_path, header=False, index=False)

############################################################################
# Create IPD to Exon dict & missing list of IPD sequences w/out Exon match #
############################################################################
print('Create IPD to Exon dict')
# diag_fasta_dict = fasta_to_df(fasta_path=miseq_all_filepath)
ipd_fasta_dict = fasta_to_df(fasta_path=ipd_all_filepath)

gen_fasta_dict = {}
exon_fasta_dict = {}

for k, v in ipd_fasta_dict.items():
    if k.endswith('-gen'):
        gen_fasta_dict[k] = v
    if k.endswith('-exon'):
        exon_fasta_dict[k] = v
ipd_to_diag_lookup = {}
for ipd_name, ipd_seq in gen_fasta_dict.items():
    match_list = []
    for diag_name, diag_seq in exon_fasta_dict.items():
        if (diag_seq in ipd_seq) or (rev_comp_st(diag_seq) in ipd_seq):
            match_list.append(diag_name)
    ipd_to_diag_lookup[ipd_name] = match_list

# out put to gson
# ipd_to_exon_lookup_path
# missing_ipd_alleles_path
with open(ipd_to_exon_lookup_path, 'w') as convert_file:
    convert_file.write(json.dumps(ipd_to_diag_lookup))
# get the missing list
missing_ipd_list = []
for ipd_seq, diag_seq_list in ipd_to_diag_lookup.items():
    if len(diag_seq_list) < 1:
        missing_ipd_list.append(ipd_seq)
missing_ipd_list.sort()
# output do data frame
df_missing_ipd = pd.DataFrame({'allele': missing_ipd_list})
df_missing_ipd.to_csv(missing_ipd_exom_alleles_path, header=False, index=False)

# out put to gson
# ipd_to_exon_lookup_path
# missing_ipd_alleles_path
with open(ipd_to_exon_lookup_path, 'w') as convert_file:
    convert_file.write(json.dumps(ipd_to_diag_lookup))
# get the missing list
missing_ipd_list = []
for ipd_seq, diag_seq_list in ipd_to_diag_lookup.items():
    if len(diag_seq_list) < 1:
        missing_ipd_list.append(ipd_seq)
missing_ipd_list.sort()
# output do data frame
df_missing_ipd = pd.DataFrame({'allele': missing_ipd_list})
df_missing_ipd.to_csv(missing_ipd_exom_alleles_path, header=False, index=False)

############################################################################
# Create Exon to IPD dict & missing list of EXON sequences w/out IPD match #
############################################################################
print('Create Exon to IPD dict')
ipd_fasta_dict = fasta_to_df(fasta_path=ipd_all_filepath)

gen_fasta_dict = {}
exon_fasta_dict = {}

for k, v in ipd_fasta_dict.items():
    if k.endswith('-gen'):
        gen_fasta_dict[k] = v
    if k.endswith('-exon'):
        exon_fasta_dict[k] = v

diag_to_ipd_lookup = {}
for diag_name, diag_seq in exon_fasta_dict.items():
    match_list = []
    for ipd_name, ipd_seq in gen_fasta_dict.items():
        if (diag_seq in ipd_seq) or (rev_comp_st(diag_seq) in ipd_seq):
            match_list.append(ipd_name)
    diag_to_ipd_lookup[diag_name] = match_list

# out put to gson
with open(exon_to_ipd_lookup_path, 'w') as convert_file:
    convert_file.write(json.dumps(diag_to_ipd_lookup))
# get the missing list
missing_diag_list = []
for diag_seq, ipd_seq_list in diag_to_ipd_lookup.items():
    if len(ipd_seq_list) < 1:
        missing_diag_list.append(diag_seq)
missing_diag_list.sort()
# output do data frame
df_missing_diag = pd.DataFrame({'allele': missing_diag_list})
df_missing_diag.to_csv(missing_exp_ipd_alleles_path, header=False, index=False)

# ipd_to_exon_lookup_path
# exon_to_ipd_lookup_path
# missing_ipd_exom_alleles_path
# missing_exp_ipd_alleles_path
###############################
#  Create haplotype Library   #
# - Based on the legacy file  #
###############################
########################################
# open the lookup .json for haplotypes #
########################################
print('Create haplotype Library')
with open(haplotype_json_path) as f_in:
    haplo_dict = json.load(f_in)

# Convert dictionary to a dataframe
df = pd.DataFrame(haplo_dict[species])
df = df.reset_index()
df.rename(columns={'index': 'HAPLOTYPE_CALL'}, inplace=True)
haplotype_cols = list(df.columns)
haplotype_cols = [x for x in haplotype_cols if x.startswith('MHC')]
# stack the data frame so all the columns line up better
df_melt = pd.melt(df, id_vars=['HAPLOTYPE_CALL', 'PREFIX'], value_vars=haplotype_cols,
                  var_name='TYPE', value_name='allele')
# drop the na entries that are created from the dict to df conversion
df_melt.dropna(inplace=True)
df_melt['CALL_COUNT'] = [len(x) for x in df_melt['allele']]
# explode the lists so each list item gets a row
df_melt = df_melt.explode('allele')
# Convert the miseq fasta to a data frame
df_diag_haplo_ref = fasta_to_df(fasta_path=miseq_legacy_short_path, as_df=True)
# merge with the data frame by name with the haplotype calling dictionary
# it is unclear if it needs the "mamu"
# if species=='Mamu':
df_melt['haplo_allele'] = [like_join(x, df_diag_haplo_ref, 'allele') for x in df_melt['allele']]
df_haplo_call = df_melt.explode('haplo_allele')
df_haplo_call = df_haplo_call[~df_haplo_call['haplo_allele'].isnull()]

##########################################################################
# map the Miseq and IPD sequences to the renamed haplotoype fasta w/BBMAP #
##########################################################################
print('map the Miseq and IPD sequences to the renamed haplotoype ')
# - First convert the Miseq to shorter headers as the headernames are >256 characters, too long for pysam
#    - Next semiperfect the Miseq to diagnostic_haplotype file
# - Finally perfectmode map the ipd to diagnostic_haplotype fasta
os.system(f'java -ea -Xmx{ram}m -Xms{ram}m -cp {cp_path} align2.BBMap build=1 \
in={miseq_all_filepath} \
ref={miseq_legacy_short_path} \
outm={miseq_haplo_sam_path} \
threads={threads} nodisk=t semiperfectmode=t ambiguous=best overwrite=True')

os.system(f'java -ea -Xmx{ram}m -Xms{ram}m -cp {cp_path} align2.BBMap build=1 \
ref={ipd_all_filepath} \
in={miseq_legacy_short_path} \
outm={ipd_haplo_sam_path} \
threads={threads} \
nodisk=t ambiguous=all maxsites=1000 perfectmode=t overwrite=t')


ipd_diag_sequence = []
haplo_ref_sequence = []
read_dict = {}
bf_in_diag = pysam.AlignmentFile(miseq_haplo_sam_path, "r", header=True)
# bf_header = pysam.AlignmentHeader(bam_inpath_diag, "r")
for r in bf_in_diag.fetch(until_eof=True):
    mapped_read = str(r).split('\t')

    if int(mapped_read[1]) in [0,256]:
        ipd_diag_sequence.append(mapped_read[0])
        haplo_ref_sequence.append(bf_in_diag.get_reference_name(int(mapped_read[2])))

bf_in_diag.close()
bf_in_ipd = pysam.AlignmentFile(ipd_haplo_sam_path, "r", header=True)
for r in bf_in_ipd.fetch(until_eof=True):
    mapped_read = str(r).split('\t')
    if int(mapped_read[1]) in [0,256]:
        ipd_diag_sequence.append(bf_in_ipd.get_reference_name(int(mapped_read[2])))
        haplo_ref_sequence.append(mapped_read[0])
bf_in_ipd.close()
df_haplo_ipd_diag = pd.DataFrame({'haplo_allele':haplo_ref_sequence, 'allele_ipd':ipd_diag_sequence})
# merge the ipd_diag witi the haplotype lookup
df_haplo_call_merge = df_haplo_call.merge(df_haplo_ipd_diag, on=['haplo_allele'], how='inner')

df_haplo_call_merge[['HAPLOTYPE_CALL','PREFIX','TYPE','allele','CALL_COUNT','allele_ipd']]
df_haplo_call_merge.to_csv(ipd_miseq_allele_haplotype_path, index=False)
print(ipd_miseq_allele_haplotype_path)



# os.system(f'java -ea -Xmx{ram}m -Xms{ram}m -cp {cp_path} align2.BBMap build=1 \
# ref={ipd_all_filepath} \
# in={miseq_legacy_short_path} \
# outm={ipd_haplo_sam_path} \
# threads={threads} \
# nodisk=t ambiguous=all maxsites=1000 perfectmode=t overwrite=t')
