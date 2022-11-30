#!/usr/bin/env python
from argparse import ArgumentParser, ArgumentTypeError
import pandas as pd
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
    parser.add_argument('--depth_input_path',
                        type=str,
                        help='relative path to depth of coverage file (zipped)',
                        required=True)
    parser.add_argument('--merged_bam_path',
                        type=str,
                        help='relative path to the merged bam file',
                        required=True)
    parser.add_argument('--filtered_allele_list_outpath',
                        type=str,
                        help='relative path to the filtered_allele_list',
                        required=True)
    parser.add_argument('--filtered_merged_bam_outpath',
                        type=str,
                        help='relative path to the filtered_merged_bam_',
                        required=True)
    parser.add_argument('--depth_threshold',
                        type=int,
                        default=10,
                        help='depth of coverage requirement for the center of the read',
                        required=False)
    parser.add_argument('--edge_distance_threshold',
                        type=int,
                        help='distance from edge to consider depth of coverage in positions',
                        default=0,
                        required=False)
    parser.add_argument('--maximum_start_position_gap',
                        type=int,
                        help='how far apart the start position can be between reads for the chimera criteria',
                        default=30,
                        required=False)
    parser.add_argument('--minimum_bridge_read_length',
                        type=int,
                        help='How mapped (untrimmed) portion of hte read must be to be considered',
                        default=70,
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
depth_input_path = args.depth_input_path
merged_bam_path = args.merged_bam_path
filtered_allele_list_outpath = args.filtered_allele_list_outpath
filtered_merged_bam_outpath = args.filtered_merged_bam_outpath
edge_distance_threshold = args.edge_distance_threshold
depth_threshold = args.depth_threshold
maximum_start_position_gap = args.maximum_start_position_gap
minimum_bridge_read_length = args.minimum_bridge_read_length

# open the depth file
##########################################
# Start filter for the depth of coverage #
##########################################
# ignoring the edges up to edge_distance_threshold positions from the edge
df = pd.read_csv(depth_input_path, header=None, sep='\t', names=['allele', 'position', 'depth'])
# figure out the distance from the edge and if the threshold for each position is met
# Find the max position per allele, as each allele can be a different length.
df['MAX_POSITION'] = df.groupby(['allele'])['position'].transform('max')
df["DISTANCE_FROM_EDGE"] = [y - x if y / 2 < x else x - 1 for x, y in zip(df['position'], df['MAX_POSITION'])]
df['DEPTH_THRESHOLD'] = [0 if (x <= depth_threshold) & (y > edge_distance_threshold) else 1 for
                         x, y in zip(df['depth'], df['DISTANCE_FROM_EDGE'])]

# Find the minimum depth threshold.  If it is zero it means there is at least 1 position that the threshold was not met
df_summary = df[['allele', 'DEPTH_THRESHOLD']].groupby(['allele'])['DEPTH_THRESHOLD'].min().reset_index()
# Filter out any alleles that the depth threshold was not met.
df_summary_pass = df_summary[df_summary['DEPTH_THRESHOLD'] > 0]
# Find a unique list of alleles, as the allele name is repeated for each positions
allele_list = list(df_summary_pass['allele'].unique())
# df_summary_median = df[['allele', 'depth']].groupby(['allele'])['depth'].median().reset_index()
# allele_median_depth_list = []
# for allele_i in allele_list:
#   df_summary_median_i = df_summary_median[df_summary_median['allele'] == allele_i]
#   median_depth_list_i = list(df_summary_median_i['allele'])
#   median_depth_i = median_depth_list_i[0]
#   allele_median_depth_list.append(median_depth_i)

#############################################
# Start the Computational Chimera Filtering #
#############################################
# Start an empty list to appned to of passing alleles
passing_alleles = []
# open the bam file merged
infile = pysam.AlignmentFile(merged_bam_path, "rb")
bam_header = infile.header.copy().to_dict()

outfile_groups = {}
i = 1

ref_dict = {}
ref_name_old = 'ZXAZXV'
position_old = 0
compare_start_old = 0
for segment in infile:
    if i % 100000 == 0:
        print(i)
    i += 1
    # Get the start and end of the ampping.
    compare_start = segment.reference_start
    compare_end = segment.reference_end
    if (compare_end is None) or (compare_start is None):
        continue
    # get the reference name (allele of interest)
    ref_name = segment.reference_name
    # check if the allele is on the passing allele list based on depth, of not skip it.
    if ref_name in allele_list:
        # This leverages a sorted file.
        # If they are the same, then we can just use the privious start position to compare
        # if they are not the same, it is a new allele, and the previous start position does not exist/not relavent
        if ref_name == ref_name_old:
            # This compares the current read to the previous read start position, if it exceeds the threshold
            # This allele will be "failed"
            # because we do not know how many more reads will be of the same allele, it continues to loop the allele.
            if compare_start - compare_start_old > maximum_start_position_gap:
                # saves the failed allele to a dictionary
                ref_dict[ref_name] = False

        else:
            # assume the reference allele it passes until we have information to fail it.
            ref_dict[ref_name] = True
        # make sure the length of mapped portion of the read criteria is passing
        if compare_end - compare_start > minimum_bridge_read_length:
            compare_start_old = compare_start
        # save the current name to the old name (previous loop) as the loop is finished.
        ref_name_old = ref_name
# close the bam file we are reading.
infile.close()

# loop through the pass/fail dictionary and make a list of passing reference alleles
for key, value in ref_dict.items():
    if value:
        passing_alleles.append(key)
# Create a csv of the passing alleles.
df_passing = pd.DataFrame({'allele': passing_alleles})


df_summary_median = df[['allele', 'depth']].groupby(['allele'])['depth'].median().reset_index()
df_passing = df_passing.merge(df_summary_median, on=['allele'], how='inner')
df_passing.to_csv(filtered_allele_list_outpath,  index=False, sep='\t')
# df_depth.to_csv(filtered_allele_list_outpath, header=False, index=False)
# make a bam file of the passing allele
infile = pysam.AlignmentFile(merged_bam_path, "rb")
bam_header = infile.header.copy().to_dict()

outfile = pysam.AlignmentFile(filtered_merged_bam_outpath, "wb", header=bam_header)
i = 1
# Since we already know which ones pass, we do not need to worry about the start posistion comparisons.
# read through the bam file again
for segment in infile:
    if i % 100000 == 0:
        print(i)
    i += 1
    ref_name = segment.reference_name
    # if the reference allele is passing, add the data to the bam file.
    if ref_name in passing_alleles:
        outfile.write(segment)
# close the read and write bam files
infile.close()
outfile.close()
