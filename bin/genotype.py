#!/usr/bin/env python3
from argparse import ArgumentParser, ArgumentTypeError
import os
import pysam
import pandas as pd
import time
import math


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
    parser.add_argument('--bam_dir',
                        type=str,
                        help='directory your input bam files reside',
                        required=True)
    parser.add_argument('--project_name',
                        type=str,
                        help='Project name for file prefix',
                        required=True)
    parser.add_argument('--config_dir',
                        type=str,
                        help='where the config files are stored',
                        required=True)

    parser.add_argument('--bait_fasta',
                        type=str,
                        help='Where you genomic-exon2-miseq combined fasta exists',
                        default=None,
                        required=False)
    parser.add_argument('--ipd_ref_matrix_dir',
                        type=str,
                        help='directory where the ipd reference matrix files exist',
                        default=None,
                        required=False)

    parser.add_argument('--unpaired_edge_threshold',
                        type=str,
                        help='unpaired_edge_threshold how far (bp) from the edge of reference sequence \
                        to allow the paired read to be unmapped (unmated) and still included',
                        default=500,
                        required=False)
    parser.add_argument('--depth_threshold',
                        type=str,
                        help='depth_threshold how minimum depth of coverage for all \
                        (except the near the edge by defied as edge distance threshold)',
                        default=3,
                        required=False)
    parser.add_argument('--edge_distance_threshold',
                        type=str,
                        help='how far from the edge of the ipd reference sequence to allow for zero depth of coverage \
                             depth of coverage does not apply close to the edge, \
                             short sequence miseq db sequences do not apply as it is set to zero.',
                        default=50,
                        required=False)
    parser.add_argument('--low_ram_mode',
                        type=str2bool,
                        help='low ram mode is slower if doing several as it will re read the ipd_ref_matrix from disk \
                             for every sample',
                        default=False,
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
bam_dir = args.bam_dir
project_name = args.project_name
config_dir = args.config_dir

unpaired_edge_threshold = args.unpaired_edge_threshold
depth_threshold = args.depth_threshold
edge_distance_threshold = args.edge_distance_threshold
low_ram_mode = args.low_ram_mode

bait_fasta = args.bait_fasta
ipd_ref_matrix_dir = args.ipd_ref_matrix_dir
if bait_fasta is None:
    bait_fasta = os.path.join(config_dir, 'bait.fasta')
if ipd_ref_matrix_dir is None:
    ipd_ref_matrix_dir = os.path.join(config_dir, 'ipd_ref_matrix')

os.makedirs(out_dir, exist_ok=True)

k = 0


def depth_screening(df, depth_threshold, edge_distance_threshold, ipd_length_dict):
    df_length = pd.DataFrame(list(ipd_length_dict.items()), columns=['ALLELE', 'MAX_POSITION'])
    df_length['MAX_POSITION'] = df_length['MAX_POSITION'] - 1
    df = df.merge(df_length, on=['ALLELE'], how='inner')
    # df['MAX_POSITION'] = df.groupby(['ALLELE'])['POSITION'].transform('max')
    df["DISTANCE_FROM_EDGE"] = [y - x if y / 2 < x else x for x, y in zip(df['POSITION'], df['MAX_POSITION'])]
    df['DEPTH_THRESHOLD'] = [0 if (x <= depth_threshold) & (y > edge_distance_threshold) else x for
                             x, y in zip(df['DEPTH'], df['DISTANCE_FROM_EDGE'])]
    df['INNER_THRESHOLD'] = [0 if (x <= depth_threshold) & (y > edge_distance_threshold) else 1 for
                             x, y in zip(df['DEPTH'], df['DISTANCE_FROM_EDGE'])]
    df['DEPTH_THRESHOLD'] = [depth_threshold if (x <= depth_threshold) & (y <= edge_distance_threshold) else x for
                             x, y in zip(df['DEPTH_THRESHOLD'], df['DISTANCE_FROM_EDGE'])]
    df['EDGE_THRESHOLD'] = [edge_distance_threshold if (x <= depth_threshold) else y for
                            x, y in zip(df['DEPTH'], df['DISTANCE_FROM_EDGE'])]
    df['SIDE'] = ['END' if y / 2 < x else 'START' for x, y in zip(df['POSITION'], df['MAX_POSITION'])]
    df['EDGE_THRESHOLD'] = df.groupby(['ALLELE', 'DISTANCE_FROM_EDGE', 'SIDE'])['EDGE_THRESHOLD'].transform('max')

    # Find the minimum depth threshold.  If it is zero it means there is at least 1 position
    # that the threshold was not met
    df['MIN_INNER_THRESHOLD'] = df.groupby(['ALLELE'])['INNER_THRESHOLD'].transform('min')
    df = df[df['MIN_INNER_THRESHOLD'] > 0]
    df_summary = df[['ALLELE', 'DEPTH_THRESHOLD', 'EDGE_THRESHOLD', 'SIDE']].groupby(['ALLELE', 'SIDE'])[
        'DEPTH_THRESHOLD', 'EDGE_THRESHOLD'].min().reset_index()
    # Filter out any alleles that the depth threshold was not met.
    df_summary_pass = df_summary[df_summary['DEPTH_THRESHOLD'] > 0]

    # Find a unique list of alleles, as the allele name is repeated for each positions
    # df_summary_pass['SAMPLE'] = num_i
    # df_summary_pass_all = pd.concat([df_summary_pass,df_summary_pass_all], ignore_index=True)
    allele_list = list(df_summary_pass['ALLELE'].unique())
    return allele_list


def filter_count_map_bam(df_bam_m, allele_list):
    df_bam_m = df_bam_m[df_bam_m['ALLELE'].isin(allele_list)]
    df_bam_m['NAME_COUNTS'] = 1 - ((df_bam_m['PAIRED'] - 1) / 2)
    df_bam_m['num_read_maps'] = df_bam_m.groupby(['NAME'])['NAME_COUNTS'].transform('sum')
    df_bam_m['unique_maps_per_allele'] = df_bam_m.groupby(['ALLELE'])['num_read_maps'].transform('min')
    df_bam_m.sort_values(by=['ALLELE', 'START', 'END'], inplace=True)
    return df_bam_m


def process_df_bam(df_bam, ipd_length_dict, sample_i, unpaired_edge_threshold):
    df_bam['SAMPLE_NUM'] = sample_i
    df_bam["DISTANCE_TO_END"] = [ipd_length_dict[y] - x for x, y in zip(df_bam['END'], df_bam['ALLELE'])]
    df_bam["INCLUDE_UNPAIRED"] = [
        1 if ((x < unpaired_edge_threshold) and z) or ((y < unpaired_edge_threshold) and not (z)) else 0 for x, y, z in
        zip(df_bam["START"], df_bam["DISTANCE_TO_END"], df_bam['REVERSED_MAPPING'])]
    df_bam = df_bam[(df_bam['PAIRED'] + df_bam["INCLUDE_UNPAIRED"]) > 1]
    df_bam['NAME_COUNTS'] = 1 - ((df_bam['PAIRED'] - 1) / 2)
    # df_bam_m['num_read_maps'] = df_bam_m.groupby(['NAME','ALLELE'])['NAME_COUNTS'].transform('sum')
    df_bam['num_read_maps'] = df_bam.groupby(['NAME'])['NAME_COUNTS'].transform('sum')
    df_bam['unique_maps_per_allele'] = df_bam.groupby(['ALLELE'])['num_read_maps'].transform('min')
    return df_bam


def fasta_length_dict(ipd_fasta_path):
    ipd_length_dict = {}
    genome_fasta_open = pysam.Fastafile(ipd_fasta_path)
    for ref in genome_fasta_open.references:
        ipd_length_dict[ref] = genome_fasta_open.get_reference_length(ref)
    genome_fasta_open.close()
    return ipd_length_dict


def get_gap_allele_list(df_bam_m, max_gap=70):
    max_gap = max_gap + 0.5

    df_bam_m_dedup = df_bam_m[
        ['NAME', 'REVERSED_MAPPING', 'START', 'END', 'PAIRED', 'ALLELE', 'SAMPLE_NUM']].drop_duplicates()
    df_bam_m_dedup['START_ALLELE_MIN'] = df_bam_m_dedup.groupby(['ALLELE', 'SAMPLE_NUM'])['START'].transform('min')
    df_bam_m_dedup['END_ALLELE_MAX'] = df_bam_m_dedup.groupby(['ALLELE', 'SAMPLE_NUM'])['END'].transform('max')

    df_bam_paired = df_bam_m_dedup[df_bam_m_dedup['PAIRED'] == 2]
    # df_bam_unpaired = df_bam_m_dedup[df_bam_m_dedup['PAIRED'] == 1]

    df_bam_paired['START_MIN'] = df_bam_paired.groupby(['ALLELE', 'SAMPLE_NUM', 'NAME'])['START'].transform('min')
    df_bam_paired['END_MAX'] = df_bam_paired.groupby(['ALLELE', 'SAMPLE_NUM', 'NAME'])['END'].transform('max')
    df_bam_paired['INSERT_LEN'] = df_bam_paired['END_MAX'] - df_bam_paired['START_MIN']
    df_bam_overlap = df_bam_paired[df_bam_paired['INSERT_LEN'] < 302]
    df_bam_no_overlap = df_bam_paired[df_bam_paired['INSERT_LEN'] > 301]
    df_bam_overlap['START'] = df_bam_overlap['END_MAX'] - 151
    df_bam_overlap['END'] = df_bam_overlap['START_MIN'] + 151
    df_bam_overlap['START'] = [y if x < y else x for x, y in
                               zip(df_bam_overlap['START_MIN'], df_bam_overlap['START_ALLELE_MIN'])]
    df_bam_overlap['END'] = [y if x > y else x for x, y in
                             zip(df_bam_overlap['END_MAX'], df_bam_overlap['END_ALLELE_MAX'])]

    df_bam_overlap = df_bam_overlap[
        ['NAME', 'REVERSED_MAPPING', 'START', 'END', 'PAIRED', 'ALLELE', 'SAMPLE_NUM']].drop_duplicates()

    df_bam_m_dedup = pd.concat([df_bam_m_dedup, df_bam_overlap], ignore_index=True)

    df_start = df_bam_m_dedup[['ALLELE', 'SAMPLE_NUM', 'START']].drop_duplicates()

    df_start.rename(columns={'START': 'POSITION'}, inplace=True)
    df_start.sort_values(['ALLELE', 'SAMPLE_NUM', 'POSITION'], inplace=True)
    df_start['GAP'] = df_start.groupby(['ALLELE', 'SAMPLE_NUM'])['POSITION'].diff()
    df_end = df_bam_m_dedup[['ALLELE', 'SAMPLE_NUM', 'END']].drop_duplicates()
    df_end.rename(columns={'END': 'POSITION'}, inplace=True)
    df_end.sort_values(['ALLELE', 'SAMPLE_NUM', 'POSITION'], inplace=True)
    df_end = df_end[['ALLELE', 'SAMPLE_NUM', 'POSITION']].drop_duplicates()
    df_end['GAP'] = df_end.groupby(['ALLELE', 'SAMPLE_NUM'])['POSITION'].diff()
    df_gap = pd.concat([df_start, df_end], ignore_index=True)
    df_gap.fillna(0, inplace=True)
    df_gap['MAX_GAP'] = df_gap.groupby(['ALLELE', 'SAMPLE_NUM'])['GAP'].transform('max')
    df_gap.sort_values(['ALLELE', 'SAMPLE_NUM', 'POSITION'], inplace=True)
    df_gap_summary = df_gap[df_gap['GAP'] == df_gap['MAX_GAP']]
    df_gap_summary = df_gap_summary.groupby(['ALLELE', 'SAMPLE_NUM', 'MAX_GAP'])['POSITION'].min().reset_index()
    df_gap_summary.rename(columns={'POSITION': 'GAP_POSITION'}, inplace=True)
    df_gap_filtered = df_gap[df_gap['MAX_GAP'] < max_gap + 0.5]

    allele_list = list(df_gap_filtered['ALLELE'].unique())
    #     df_bam_m = df_bam_m[df_bam_m['ALLELE'].isin(allele_list) & df_bam_m['ALLELE_2'].isin(allele_list)]
    return allele_list, df_gap, df_gap_summary


def make_range_from_fasta_lengths_dict(ipd_length_dict, allele_list):
    df_depth_range = pd.DataFrame()
    for allele_i in allele_list:
        df_depth_range_i = pd.DataFrame({'ALLELE': allele_i, 'POSITION': list(range(0, ipd_length_dict[allele_i]))})
        df_depth_range = pd.concat([df_depth_range, df_depth_range_i], ignore_index=True)
    return df_depth_range


def get_depth_miseq(df_bam_m, ipd_length_dict, sample_i):
    df_bam_m_dedup = df_bam_m[['SAMPLE_NUM', 'ALLELE', 'NAME', 'START', 'END']].drop_duplicates()

    df_bam_m_dedup['POSITION'] = [list(range(x, y)) for x, y in zip(df_bam_m_dedup['START'], df_bam_m_dedup['END'])]
    df_depth = df_bam_m_dedup.explode('POSITION')

    df_depth = df_depth.groupby(['SAMPLE_NUM', 'ALLELE', 'POSITION']).agg({'START': 'count'}).rename(
        columns={'START': 'DEPTH'}).reset_index()
    #     df_depth['DEPTH_ADJ'] = df_depth['DEPTH'] * df_depth['DEPTH'] / df_depth['num_read_maps']
    #     df_depth['DEPTH_ADJ'] = df_depth['DEPTH_ADJ'].round(2)
    #     df_depth['DEPTH_RATIO'] = df_depth['num_read_maps']  / df_depth['DEPTH']
    #     df_depth['DEPTH_RATIO'] = df_depth['DEPTH_RATIO'].round(2)
    df_depth['HAPLOTYPE_GROUP'] = [x.split('__')[0] for x in df_depth['ALLELE']]
    allele_list = list(df_depth['ALLELE'].unique())
    df_depth_range = make_range_from_fasta_lengths_dict(ipd_length_dict, allele_list)
    df_depth = df_depth.merge(df_depth_range, on=['ALLELE', 'POSITION'], how='outer')
    df_depth['SAMPLE_NUM'] = sample_i
    df_depth['HAPLOTYPE_GROUP'] = [x.split('__')[0] for x in df_depth['ALLELE']]
    df_depth.fillna(value=0, inplace=True)
    return df_depth


def get_depth(df_bam_m, ipd_length_dict, fillzero=True):
    df_bam_m_dedup = df_bam_m[
        ['SAMPLE_NUM', 'ALLELE', 'NAME', 'START', 'END', 'unique_maps_per_allele', 'num_read_maps']].drop_duplicates()
    df_bam_m_dedup['POSITION'] = [list(range(x, y)) for x, y in zip(df_bam_m_dedup['START'], df_bam_m_dedup['END'])]
    df_depth = df_bam_m_dedup.explode('POSITION')
    df_depth = df_depth.groupby(['SAMPLE_NUM', 'ALLELE', 'POSITION', 'unique_maps_per_allele']).agg(
        {'START': 'count', 'num_read_maps': 'sum'}).rename(columns={'START': 'DEPTH'}).reset_index()
    df_depth['DEPTH_ADJ'] = df_depth['DEPTH'] * df_depth['DEPTH'] / df_depth['num_read_maps']
    df_depth['DEPTH_ADJ'] = df_depth['DEPTH_ADJ'].round(2)
    df_depth['DEPTH_RATIO'] = df_depth['num_read_maps'] / df_depth['DEPTH']
    df_depth['DEPTH_RATIO'] = df_depth['DEPTH_RATIO'].round(2)
    df_depth['HAPLOTYPE_GROUP'] = [x.split('__')[0] for x in df_depth['ALLELE']]
    allele_list = list(df_depth['ALLELE'].unique())
    sample_i = list(df_depth['SAMPLE_NUM'].unique())[0]
    if fillzero:
        df_depth_range = make_range_from_fasta_lengths_dict(ipd_length_dict, allele_list)
        df_depth = df_depth.merge(df_depth_range, on=['ALLELE', 'POSITION'], how='outer')
    df_depth['SAMPLE_NUM'] = sample_i
    df_depth['HAPLOTYPE_GROUP'] = [x.split('_')[0] for x in df_depth['ALLELE']]
    df_depth.fillna(value=0, inplace=True)
    df_depth['unique_maps_per_allele'] = df_depth.groupby(['ALLELE'])['unique_maps_per_allele'].transform('max')
    return df_depth


def get_normalized_median_by_allele_sample(df_depth, median_groups):
    df_summary = df_depth.groupby(['ALLELE', 'SAMPLE_NUM', 'HAPLOTYPE_GROUP']).agg(
        {'unique_maps_per_allele': 'median', 'DEPTH_ADJ': 'median'}).reset_index()
    if len(median_groups) > 0:
        df_summary_median = df_summary[df_summary['HAPLOTYPE_GROUP'].isin(median_groups)]
    else:
        df_summary_median = df_summary
    median_depth = df_summary_median['DEPTH_ADJ'].median()
    df_summary['DEPTH_NORM'] = df_summary['DEPTH_ADJ'] / median_depth
    df_summary['DEPTH_NORM'] = df_summary['DEPTH_NORM'].round(3)
    return df_summary


def bam_to_df(bam_path):
    samfile = pysam.AlignmentFile(bam_path, 'rb')
    allele_list = []
    name_list = []
    # paired_list = []
    start_list = []
    end_list = []
    reverse_list = []
    for segment in samfile:
        if segment.is_unmapped:
            continue
        allele_list.append(segment.reference_name)
        name_list.append(segment.query_name)
        start_list.append(segment.reference_start)
        end_list.append(segment.reference_end)
        reverse_list.append(segment.is_reverse)
    #         if not segment.is_paired or segment.mate_is_unmapped or segment.is_duplicate:
    #             paired_list.append(0)
    #         else:
    #             paired_list.append(1)
    segment_count = len(name_list)
    df_bam_in = pd.DataFrame({'NAME': name_list, 'ALLELE_2': allele_list,
                              'START_2': start_list, 'END_2': end_list, 'REVERSED_MAPPING': reverse_list})
    df_bam_in.sort_values(by=['ALLELE_2', 'START_2', 'END_2'], inplace=True)
    return df_bam_in, segment_count


def diff_calc(w, x, y, z):
    if x != w:
        return 1
    if y < z:
        return 0
    return 1


def quick_filter_depth(df_bam, edge_distance_threshold):
    # quick filter by distance to the edge
    # Since nearly every allele has something in common, we have 50% viable alleles.
    # Expanding to get depth of coverage on 2000-4000 alleles is very time consuming
    # Instead we can quickly filter out 80-90% of bad alleles (down to 200-400 viable alleles)
    # and subsequent steps are faster
    # This was not an issue using already filtered for depth of coverage files, (which was a very slow process)
    df_bam['min_dist_start'] = df_bam.groupby(['ALLELE', 'SAMPLE_NUM'])['START'].transform('min')
    df_bam['min_dist_end'] = df_bam.groupby(['ALLELE', 'SAMPLE_NUM'])['DISTANCE_TO_END'].transform('min')
    df_bam = df_bam[df_bam['min_dist_start'] <= edge_distance_threshold]
    df_bam = df_bam[df_bam['min_dist_end'] <= edge_distance_threshold]
    print('Finished IPD quick filter')
    # look for gaps in depth of coverage
    # this does not look for the ending or start distance as the previous step already addressed it
    # Quick way to find a differential by shifting the row up and making sure the alleles are still the same.
    df_bam.sort_values(by=['ALLELE', 'START', 'END'], inplace=True)
    allele_list = list(df_bam['ALLELE'])
    start_list = list(df_bam['START'])
    start_list.pop(0)
    start_list.append(0)
    allele_list.pop(0)
    allele_list.append(allele_list[-1])
    df_bam['DIFF'] = [diff_calc(w, x, y, z) for w, x, y, z in
                      zip(df_bam['ALLELE'], allele_list, df_bam['END'], start_list)]

    df_diff = df_bam.groupby(['ALLELE'])['DIFF'].min().reset_index()

    allele_list = list(df_diff[df_diff['DIFF'] > 0]['ALLELE'])
    return allele_list, df_bam, df_diff


# 0 qname
# 1 flag
# 2 ref name
# 3 ref start
# 4 mapping quality
# 5 cigar string
# 6 ref id of mate
# 7 next ref start
# 8 inferred insert size overlap (reverse is negative, (END2- START1) non paired is 0)
# 9 query sequence
# 10 query sequence score


def bam_to_df_2(bam_path):
    sam_file = pysam.AlignmentFile(bam_path, 'rb')
    name_list = []
    segment_list = []
    reverse_list = []
    mapping_quality_list = []
    reference_start_list = []
    query_sequence_list = []
    query_qualities_list = []
    reference_end_list = []
    query_flag_list = []
    i = 0
    for segment in sam_file:
        name_list.append(segment.query_name)
        if not isinstance(segment.reference_end, int):
            reverse_list.append(None)
        else:
            reverse_list.append(segment.is_reverse)
        mapping_quality_list.append(segment.mapping_quality)
        reference_start_list.append(segment.reference_start)
        reference_end_list.append(segment.reference_end)
        query_sequence_list.append(segment.query_sequence)
        query_qualities_list.append(segment.to_string().split('\t')[10])
        query_flag_list.append(segment.to_string().split('\t')[1])
    df = pd.DataFrame({'NAME': name_list,
                       # 'SEGMENT':segment_list,
                       'REVERSED_MAPPING': reverse_list,
                       'MAPPING_QUALITY': mapping_quality_list,
                       'REFERENCE_START': reference_start_list,
                       'REFERENCE_END': reference_end_list,
                       'QUERY_SEQUENCE': query_sequence_list,
                       'QUERY_QUALITIES': query_qualities_list,
                       'QUERY_FLAG': query_flag_list
                       })
    return df


def calc_flag(read_mapped, mate_unmapped=True, reverse=True, first_in_pair=1, read_paired=1):
    # assumes all paired reads
    # read map is 4, mate unmapped is 8 (4+4)
    # reverse is 16/32 (4+16 + 64)
    # 2 is read mapped in proper pair
    flag = read_paired
    if read_mapped and not mate_unmapped:
        flag += 2
    if not read_mapped:
        flag += 4
    if mate_unmapped:
        flag += 8
    if reverse:
        flag += 16
    if not mate_unmapped and not reverse:
        flag += 32
    if first_in_pair == 1:
        flag += 64
    if first_in_pair > 1:
        flag += 128
    return flag


def cigar_string(start, end):
    seg_len = end - start
    if seg_len < 151:
        split_len = 151 - seg_len
        if start > 0:
            return '{0}={1}S'.format(int(seg_len), int(split_len))
        return '{0}S{1}='.format(int(split_len), int(seg_len))
    return '151='


# 1	1	Read paired
# 2	2	Read mapped in proper pair
# 3	4	Read unmapped
# 4	8	Mate unmapped
# 5	16	Read reverse strand
# 6	32	Mate reverse strand
# 7	64	First in pair
# 8	128	Second in pair
# 9	256	Not primary alignment
# 10	512	Read fails platform/vendor quality checks
# 11	1024	Read is PCR or optical duplicate
# 12	2048	Supplementary alignment


def makebam(df_bam_in, df_mapped, ipd_length_dict, out_sam):
    df_read_names = df_mapped[['NAME']].drop_duplicates()
    df_mapped_dedup = df_mapped[['NAME', 'REVERSED_MAPPING', 'START', 'END', 'PAIRED', 'ALLELE']].drop_duplicates()
    df_mapped_bam = df_bam_in.merge(df_read_names, on=['NAME'], how='inner')

    df_mapped_bam = df_mapped_bam.merge(
        df_mapped_dedup[['NAME', 'REVERSED_MAPPING', 'START', 'END', 'PAIRED', 'ALLELE']],
        on=['NAME', 'REVERSED_MAPPING'],
        how='left')

    df_unmapped = df_mapped_bam[df_mapped_bam['PAIRED'].isna()]
    # Get mapped as they are treated differently
    df_unpaired = df_mapped_bam[df_mapped_bam['PAIRED'] == 1]
    # Get mapped as they are treated differently
    df_paired = df_mapped_bam[df_mapped_bam['PAIRED'] == 2]

    df_unmapped['CIGAR'] = '*'
    df_unmapped['REVERSED_MAPPING'] = False
    df_unmapped['FIRST_IN_PAIR'] = 2
    df_unmapped['READ_MAPPED'] = False
    df_unmapped['MATE_UNMAPPED'] = False
    df_unmapped['MAPPING_QUALITY'] = 0
    df_unmapped['FLAG'] = [calc_flag(read_mapped=w,
                                     mate_unmapped=x,
                                     reverse=y,
                                     first_in_pair=z,
                                     read_paired=1) for w, x, y, z in zip(df_unmapped['READ_MAPPED'],
                                                                          df_unmapped['MATE_UNMAPPED'],
                                                                          df_unmapped['REVERSED_MAPPING'],
                                                                          df_unmapped['FIRST_IN_PAIR'])]

    df_unmapped.drop(columns=['ALLELE', 'START', 'END'], inplace=True)
    df_unpaired_m = df_unpaired[['NAME', 'ALLELE', 'START', 'END']]
    df_unmapped = df_unmapped.merge(df_unpaired_m, on=['NAME'], how='inner')
    df_unmapped['INFERRED_OVERLAP'] = 0

    # df_ipd_unpaired['REVERSED_MAPPING'] = False
    df_unpaired['CIGAR'] = [cigar_string(start, end) for start, end in zip(df_unpaired['START'],
                                                                           df_unpaired['END'])]
    df_unpaired['FIRST_IN_PAIR'] = 1
    df_unpaired['READ_MAPPED'] = True
    df_unpaired['MATE_UNMAPPED'] = True
    df_unpaired['INFERRED_OVERLAP'] = 0
    df_unpaired['FLAG'] = [calc_flag(read_mapped=w,
                                     mate_unmapped=x,
                                     reverse=y,
                                     first_in_pair=z,
                                     read_paired=1) for w, x, y, z in zip(df_unpaired['READ_MAPPED'],
                                                                          df_unpaired['MATE_UNMAPPED'],
                                                                          df_unpaired['REVERSED_MAPPING'],
                                                                          df_unpaired['FIRST_IN_PAIR'])]

    df_paired['CIGAR'] = [cigar_string(start, end) for start, end in zip(df_paired['START'],
                                                                         df_paired['END'])]
    df_paired['FIRST_IN_PAIR'] = df_paired.groupby(['ALLELE', 'NAME'])['START'].rank(method='first')
    df_paired['READ_MAPPED'] = True
    df_paired['MATE_UNMAPPED'] = False
    df_paired['START_1'] = df_paired.groupby(['ALLELE', 'NAME'])['START'].transform('min')
    df_paired['END_1'] = df_paired.groupby(['ALLELE', 'NAME'])['END'].transform('min')
    df_paired['INFERRED_OVERLAP'] = df_paired['END_1'] - df_paired['START_1']
    df_paired['FLAG'] = [calc_flag(read_mapped=w,
                                   mate_unmapped=x,
                                   reverse=y,
                                   first_in_pair=z,
                                   read_paired=1) for w, x, y, z in zip(df_paired['READ_MAPPED'],
                                                                        df_paired['MATE_UNMAPPED'],
                                                                        df_paired['REVERSED_MAPPING'],
                                                                        df_paired['FIRST_IN_PAIR'])]

    df_all_pair = pd.concat([df_paired, df_unpaired, df_unmapped], ignore_index=True)
    df_all_pair.sort_values(by=['ALLELE', 'END', 'START', 'NAME', 'FIRST_IN_PAIR'], inplace=True)

    allele_list = list(set(df_all_pair['ALLELE']))
    allele_list.sort()
    with open(out_sam, mode='w') as fout:
        fout.write('@HD\tVN:1.4\tSO:unsorted\n')
        for allele_i in allele_list:
            header_str = '@SQ\tSN:{0}\tLN:{1}\n'.format(allele_i, int(ipd_length_dict[allele_i]))
            fout.write(header_str)
        for idx, row in df_all_pair.iterrows():
            sam_line_list = list()
            sam_line_list.append(row['NAME'])
            sam_line_list.append(format(int(row['FLAG'])))
            sam_line_list.append(row['ALLELE'])
            sam_line_list.append(format(int(row['START'] + 1)))
            sam_line_list.append(format(int(row['MAPPING_QUALITY'])))
            sam_line_list.append(row['CIGAR'])
            sam_line_list.append(row['ALLELE'])
            sam_line_list.append(format(int(row['START'] + 1)))
            sam_line_list.append(format(int(row['INFERRED_OVERLAP'])))
            sam_line_list.append(row['QUERY_SEQUENCE'])
            sam_line_list.append(row['QUERY_QUALITIES'])
            sam_line = '{0}\n'.format('\t'.join(sam_line_list))
            fout.write(sam_line)

    # 0 qname
    # 1 flag
    # 2 ref name
    # 3 ref start
    # 4 mapping quality
    # 5 cigar string
    # 6 ref id of mate
    # 7 next ref start
    # 8 inferred insert size overlap (reverse is negative, (END2- START1) non paired is 0)
    # 9 query sequence
    # 10 query sequence score

    out_bam = '{0}bam'.format(out_sam[:-3])
    create_bam_cmd = 'samtools view -o {0}.bam {1}'.format(out_bam, out_sam)
    os.system(create_bam_cmd)
    os.remove(out_sam)
    sort_bam_cmd = 'samtools sort -o {0} {0}.bam'.format(out_bam)
    os.system(sort_bam_cmd)
    os.remove('{0}.bam'.format(out_bam))
    print(out_bam)
    return out_bam


def annotate_depth(x, y, z):
    y = str(int(y))
    y = y.zfill(2)
    ambig = '00'
    if z < 0.75:
        ambig = '02'
    return '{0}.{1}{2}'.format(int(x), y, ambig)


def annotate_depth_miseq(x, y):
    ambig = '00'
    if y > 1:
        ambig = '02'
    return '{0}.99{1}'.format(int(x), ambig)


def short_sequence_analysis(df_ss2, sample, suffix):
    df_ss = df_ss2.copy()
    allele_s_list = list(df_ss['ALLELE'].unique())
    df_ss['SAMPLE_NUM'] = sample

    df_ss = filter_count_map_bam(df_ss, allele_s_list)
    print(df_ss)

    df_depth_ss = get_depth(df_ss,
                            ipd_length_dict,
                            fillzero=True)
    print('Finished {0} depth 1 of 3: {1}'.format(suffix, sample))
    allele_s_list = depth_screening(df_depth_ss,
                                    depth_threshold=2,
                                    edge_distance_threshold=-1,
                                    ipd_length_dict=ipd_length_dict)
    df_ss = filter_count_map_bam(df_ss, allele_s_list)
    df_depth_ss = get_depth(df_ss, ipd_length_dict, fillzero=True)
    print('Finished {0} depth 2 of 3: {1}'.format(suffix, sample))
    allele_s_list = depth_screening(df_depth_ss,
                                    depth_threshold=2,
                                    edge_distance_threshold=-1,
                                    ipd_length_dict=ipd_length_dict)
    df_ss = filter_count_map_bam(df_ss, allele_s_list)

    df_ss_2 = df_ss2[df_ss2['ALLELE'].isin(allele_s_list)][
        ['NAME', 'ALLELE_2', 'START_2', 'END_2', 'REVERSED_MAPPING', 'ALLELE']]
    df_ss_2.rename(columns={'ALLELE': 'ALLELE_3'}, inplace=True)

    df_ss_m = df_ss.merge(df_ss_2, on=['NAME', 'ALLELE_2', 'START_2', 'END_2', 'REVERSED_MAPPING'])
    df_ss_m.drop(columns=['ALLELE_2', 'START_2', 'END_2'], inplace=True)
    df_ss_m.rename(columns={'ALLELE_3': 'ALLELE_2'}, inplace=True)
    df_ss_m = filter_count_map_bam(df_ss_m, allele_s_list)
    os.makedirs(os.path.join(out_dir, 'expanded_maps_{0}'.format(suffix)), exist_ok=True)
    df_ss_m.to_csv(os.path.join(out_dir,
                                'expanded_maps_{0}'.format(suffix),
                                '{0}_expanded_maps_{1}.csv'.format(sample, suffix)),
                   index=False)
    df_depth_ss = get_depth(df_ss, ipd_length_dict, fillzero=True)
    os.makedirs(os.path.join(out_dir, 'depth_{0}'.format(suffix)), exist_ok=True)
    df_depth_ss.to_csv(os.path.join(out_dir,
                                    'depth_{0}'.format(suffix),
                                    '{0}_depth_{1}.csv'.format(sample, suffix)),
                       index=False)
    print('Finished {0}  depth 3 of 3: {1}'.format(suffix, sample))
    df_median_ss = get_normalized_median_by_allele_sample(df_depth_ss, median_groups=[])
    os.makedirs(os.path.join(out_dir, 'normalized_median_{0}'.format(suffix)), exist_ok=True)
    df_median_ss.to_csv(os.path.join(out_dir,
                                     'normalized_median_{0}'.format(suffix),
                                     '{0}_norm_median_{1}.csv'.format(sample, suffix)),
                        index=False)
    df_median_ss['SAMPLE_MEDIAN'] = df_median_ss.groupby('SAMPLE_NUM')['DEPTH_ADJ'].transform('median')
    df_median_ss['DEPTH_NORM'] = df_median_ss['DEPTH_ADJ'] / df_median_ss['SAMPLE_MEDIAN']
    df_median_ss['DEPTH_NORM'] = df_median_ss['DEPTH_NORM'].round(3)
    df_median_ss.drop(columns='SAMPLE_MEDIAN', inplace=True)
    df_median_ss['DEPTH_ADJ'] = [annotate_depth_miseq(x, y) for x, y in zip(df_median_ss['DEPTH_ADJ'],
                                                                    df_median_ss[
                                                                        'unique_maps_per_allele'])]
    df_median_ss['DEPTH_ADJ'] = df_median_ss['DEPTH_ADJ'].astype(float)
    print('--- {0} complete {1} ---'.format(suffix, sample))
    return df_median_ss, df_depth_ss, df_ss_m


ipd_matrix_df_list = []
# def pair_reads_from_bam(out_dir, num_i, ipd_length_dict, unpaired_edge_threshold=500):


pd.options.mode.chained_assignment = None
t1 = time.time()
t2 = time.time()
# bam_name_path = os.path.join(out_dir, '{0}.merged.bam'.format(num_i))

ipd_length_dict = fasta_length_dict(bait_fasta)

ipd_matrix_list = os.listdir(ipd_ref_matrix_dir)
ipd_matrix_list = [os.path.join(ipd_ref_matrix_dir, x) for x in ipd_matrix_list if
                   x.endswith('_matrix.csv.gz') and not x.startswith('._')]
i = 1
file_count = len(ipd_matrix_list)

bam_filelist = os.listdir(bam_dir)
bam_filelist = [os.path.join(bam_dir, x) for x in bam_filelist if x.endswith('.bam') and not x.startswith('._')]

df_median_exon = pd.DataFrame()
df_median_miseq = pd.DataFrame()
df_norm_median = pd.DataFrame()
df_read_ct = pd.DataFrame()
df_gap_summary = pd.DataFrame()
# ipd_matrix_df_list = []
bam_count = len(bam_filelist)
bam_iter = 1
for bam_file_i in bam_filelist:
    print('--- {0} of {1} ---'.format(bam_iter, bam_count))
    bam_iter += 1
    sample_i = os.path.basename(bam_file_i).split('.')[0]
    print(bam_file_i)
    df_bam_in, segment_count = bam_to_df(bam_file_i)
    df_all = pd.DataFrame()

    if low_ram_mode:
        for ipd_matrix_i in ipd_matrix_list:
            print('low ram mode', i, file_count)
            i += 1
            df_mat = pd.read_csv(ipd_matrix_i)
            df_all_i = df_bam_in.merge(df_mat, on=['ALLELE_2', 'START_2', 'END_2'], how='inner')
            df_all = pd.concat([df_all, df_all_i], ignore_index=True)
    else:
        if len(ipd_matrix_df_list) < 1:
            for ipd_matrix_i in ipd_matrix_list:
                print('read_files high ram mode', i, file_count)
                i += 1
                df_mat = pd.read_csv(ipd_matrix_i)
                ipd_matrix_df_list.append(df_mat)
        i = 1
        for df_mat in ipd_matrix_df_list:
            print('merge files high ram mode', i, file_count)
            i += 1
            df_all_i = df_bam_in.merge(df_mat, on=['ALLELE_2', 'START_2', 'END_2'], how='inner')
            df_all = pd.concat([df_all, df_all_i], ignore_index=True)

    df_all['PAIRED'] = df_all.groupby(['NAME', 'ALLELE'])['REVERSED_MAPPING'].transform('count')
    df_all[['START', 'END', 'START_2', 'END_2']] = df_all[['START', 'END', 'START_2', 'END_2']].astype(int)
    df_miseq = df_all[df_all['ALLELE'].str.endswith('-miseq')]
    df_ipd = df_all[df_all['ALLELE'].str.endswith('-gen')]
    df_exon = df_all[df_all['ALLELE'].str.endswith('-exon')]
    df_bam = process_df_bam(df_bam=df_ipd,
                            ipd_length_dict=ipd_length_dict,
                            sample_i=sample_i,
                            unpaired_edge_threshold=unpaired_edge_threshold)

    print('Finished IPD pair_reads_from_bam: {0}'.format(sample_i))

    allele_list, df_bam, df_diff = quick_filter_depth(df_bam, edge_distance_threshold=edge_distance_threshold)
    df_bam = filter_count_map_bam(df_bam, allele_list)

    # dropping columns does not reduce RAM as the object already reserved the space
    # Saving over the same TYPE Does.
    # changing Type also does not save ram
    # but if you do a filter step and garbage collect to a new object it will temporarily use 1x-2x ram, then use less.
    #     allele_list = depth_screening(df_depth, depth_threshold, edge_distance_threshold, ipd_length_dict)
    # print(len(allele_list))
    print('Finished quick screening depth 1 of 4: {0}'.format(sample_i))
    # df_bam = df_bam[df_bam['ALLELE'].isin(allele_list)]
    # df_depth.to_csv(os.path.join(out_dir, '{0}_m_detph.csv'.format(num_i)), index=False)
    df_depth = get_depth(df_bam, ipd_length_dict, fillzero=True)
    # df_depth.to_csv(os.path.join(out_dir, '{0}_m_quick_depth.csv'.format(sample_i)), index=False)
    print('Finished genomic depth with zero fill and depth of coverage2 of 4: {0}'.format(sample_i))
    allele_list = depth_screening(df_depth, depth_threshold, edge_distance_threshold, ipd_length_dict)
    # print(len(allele_list))
    df_bam = filter_count_map_bam(df_bam, allele_list)

    # Remove Unpaired reads that are paired with another viable allele
    df_bam['MAX_PAIRS'] = df_bam.groupby(['NAME'])['PAIRED'].transform('max')
    df_bam = df_bam[df_bam['MAX_PAIRS'] == df_bam['PAIRED']]

    df_depth = get_depth(df_bam, ipd_length_dict, fillzero=True)
    print('Finished genomic depth 3 of 4: {0}'.format(sample_i))
    allele_list = depth_screening(df_depth, depth_threshold, edge_distance_threshold, ipd_length_dict)
    df_bam = filter_count_map_bam(df_bam, allele_list)
    # print(len(allele_list))
    # join matching alleles
    allele_list_, df_gap, df_gap_summary_i = get_gap_allele_list(df_bam)
    df_gap_summary = pd.concat([df_gap_summary, df_gap], ignore_index=True)
    df_ipd_2 = df_ipd[df_ipd['ALLELE'].isin(allele_list)][
        ['NAME', 'ALLELE_2', 'START_2', 'END_2', 'REVERSED_MAPPING', 'ALLELE']]
    df_ipd_2.rename(columns={'ALLELE': 'ALLELE_3'}, inplace=True)
    df_bam_m = df_bam.merge(df_ipd_2, on=['NAME', 'ALLELE_2', 'START_2', 'END_2', 'REVERSED_MAPPING'])
    df_bam_m.drop(columns=['ALLELE_2', 'START_2', 'END_2'], inplace=True)
    df_bam_m.rename(columns={'ALLELE_3': 'ALLELE_2'}, inplace=True)
    os.makedirs(os.path.join(out_dir, 'expanded_maps_ipd'), exist_ok=True)
    df_bam_m.to_csv(os.path.join(out_dir, 'expanded_maps_ipd', '{0}_expanded_maps_ipd.csv'.format(sample_i)),
                    index=False)
    df_depth = get_depth(df_bam_m, ipd_length_dict, fillzero=True)
    print('Finished genomic depth 4 of 4: {0}'.format(sample_i))
    os.makedirs(os.path.join(out_dir, 'depth_gen'), exist_ok=True)
    df_depth.to_csv(os.path.join(out_dir, 'depth_gen', '{0}_depth_gen.csv'.format(sample_i)), index=False)

    df_norm_median_i = get_normalized_median_by_allele_sample(df_depth, median_groups=['Mamu-A1', 'Mamu-B'])
    df_norm_median_i = df_norm_median_i.merge(df_gap_summary_i, on=['ALLELE',
                                                                    'SAMPLE_NUM'], how='inner')
    os.makedirs(os.path.join(out_dir, 'normalized_median_gen'), exist_ok=True)
    df_norm_median_i.to_csv(os.path.join(out_dir, 'normalized_median_gen', '{0}_norm_median.csv'.format(sample_i)),
                            index=False)
    # print(allele_list)
    print('--- Genomic complete {0} ---'.format(sample_i))
    df_median_miseq_i, df_depth_miseq, df_miseq_m = short_sequence_analysis(df_ss2=df_miseq,
                                                                            sample=sample_i,
                                                                            suffix='miseq')

    df_median_exon_i, df_depth_exon, df_exon_m = short_sequence_analysis(df_ss2=df_exon,
                                                                         sample=sample_i,
                                                                         suffix='exon')

    df_read_ct_i = pd.DataFrame({'sample_read_ct': [segment_count], 'gs_id': [sample_i]})
    df_read_ct = pd.concat([df_read_ct, df_read_ct_i], ignore_index=True)

    df_norm_median = pd.concat([df_norm_median, df_norm_median_i], ignore_index=True)
    df_median_miseq = pd.concat([df_median_miseq, df_median_miseq_i], ignore_index=True)
    df_median_exon = pd.concat([df_median_exon, df_median_exon_i], ignore_index=True)

    print('--- start make filtered and expanded bam files ---')
    df_bam_in = bam_to_df_2(bam_file_i)
    os.makedirs(os.path.join(out_dir, 'miseq_bam'), exist_ok=True)
    os.makedirs(os.path.join(out_dir, 'gen_bam'), exist_ok=True)
    os.makedirs(os.path.join(out_dir, 'exon_bam'), exist_ok=True)
    # miseq db
    makebam(df_bam_in=df_bam_in,
            df_mapped=df_miseq_m,
            ipd_length_dict=ipd_length_dict,
            out_sam=os.path.join(out_dir, 'miseq_bam', '{0}_miseq.sam'.format(sample_i)))
    # exon db
    makebam(df_bam_in=df_bam_in,
            df_mapped=df_exon_m,
            ipd_length_dict=ipd_length_dict,
            out_sam=os.path.join(out_dir, 'exon_bam', '{0}_exon.sam'.format(sample_i)))
    # gen db
    makebam(df_bam_in=df_bam_in,
            df_mapped=df_bam_m,
            ipd_length_dict=ipd_length_dict,
            out_sam=os.path.join(out_dir, 'gen_bam', '{0}_gen.sam'.format(sample_i)))
    print('--- SAMPLE complete {0} ---'.format(sample_i))
    print(time.time() - t2)
    t2 = time.time()

df_norm_median['DEPTH_ADJ'] = [annotate_depth(x, y, z) for x, y, z in zip(df_norm_median['DEPTH_ADJ'],
                                                                          df_norm_median['unique_maps_per_allele'],
                                                                          df_norm_median['DEPTH_NORM'])]

df_norm_median['DEPTH_ADJ'] = df_norm_median['DEPTH_ADJ'].astype(float)
df_norm_median.to_csv(os.path.join(out_dir, '{0}_norm_median_gen.csv'.format(project_name)), index=False)
df_median_miseq.to_csv(os.path.join(out_dir, '{0}_norm_median_miseq.csv'.format(project_name)), index=False)
df_median_exon.to_csv(os.path.join(out_dir, '{0}_norm_median_exon.csv'.format(project_name)), index=False)
df_read_ct.to_csv(os.path.join(out_dir, '{0}_read_ct.csv'.format(project_name)), index=False)
df_gap_summary.to_csv(os.path.join(out_dir, '{0}_gaps.csv'.format(project_name)), index=False)
df_norm_median['DB'] = 'gen'
df_median_miseq['DB'] = 'miseq'
df_median_exon['DB'] = 'exon'
df_norm_median_all = pd.concat([df_norm_median, df_median_miseq, df_median_exon], ignore_index=True)
df_norm_median_all.to_csv(os.path.join(out_dir, '{0}_norm_median_all.csv'.format(project_name)), index=False)
print(out_dir)
print(time.time() - t1)
print('---Pipeline Complete--')
