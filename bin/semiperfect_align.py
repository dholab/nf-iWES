#!/usr/bin/env python
from argparse import ArgumentParser, ArgumentTypeError
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
    parser.add_argument('--fastq_dir',
                        type=str,
                        help='directory to your input fastq files',
                        required=True)
    parser.add_argument('--bam_dir',
                        type=str,
                        help='directory your output bam files will reside',
                        required=True)
    parser.add_argument('--cp_dir',
                        type=str,
                        help='directory your bbmap executables are',
                        required=True)

    parser.add_argument('--config_dir',
                        type=str,
                        help='where the config files are stored',
                        required=True)

    parser.add_argument('--bait_fasta',
                        type=str,
                        help='Where you ipd-diag combined fasta exists',
                        default=None,
                        required=False)
    parser.add_argument('--threads',
                        type=int,
                        help='Where you ipd-diag combined fasta exists',
                        default=1,
                        required=False)
    parser.add_argument('--ram',
                        type=int,
                        help='Where you ipd-diag combined fasta exists',
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
cp_dir = args.cp_dir
fastq_dir = args.fastq_dir
bam_dir = args.bam_dir
config_dir = args.config_dir
bait_fasta = args.bait_fasta
threads = args.threads
ram = args.ram

if bait_fasta is None:
    bait_fasta = os.path.join(config_dir, 'bait.fasta')

fastq_filelist = os.listdir(fastq_dir)
fastq_filelist = [os.path.join(fastq_dir, x) for x in fastq_filelist if
                  x.endswith('R1_001.fastq.gz') and not x.startswith('._')]
fastq2_filelist = [x.replace('_R1_', '_R2_') for x in fastq_filelist]
sample_list = [os.path.basename(x).split('_')[0] for x in fastq_filelist]
# os.makedirs(bam_dir, exist_ok=True)
for in_1, in_2, sample_i in zip(fastq_filelist, fastq2_filelist, sample_list):
    print(in_1)
    if os.path.exists(os.path.join(bam_dir, '{0}.bam'.format(sample_i))):
        print(os.path.join(bam_dir, '{0}.bam'.format(sample_i)))
        continue
    if not os.path.exists(in_1) and not os.path.exists(in_2):
        print(in_1)
        continue
    os.system('java -ea -Xmx{0}m -Xms{0}m -cp {1} align2.BBMap build=1 \
    in={2} \
    in2={3} \
    ref={4} \
    outm={5}/{6}.bam \
    semiperfectmode=t \
    threads={7} \
    nodisk=f'.format(int(ram), cp_dir, in_1, in_2, bait_fasta, bam_dir, sample_i, int(threads)))