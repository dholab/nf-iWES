#!/usr/bin/env bash

usage () { echo "Usage : $0 \
-d --data_dir FASTQ_PATH of <FASTQ_PATH>"; }
echo $1

echo $@

for arg in "$@"; do
  shift
  case "$arg" in
    "--config_dir") set -- "$@" "-c" ;;
    "--cp_dir") set -- "$@" "-e" ;;
    "--fastq_dir") set -- "$@" "-f" ;;
    "--threads") set -- "$@" "-t" ;;
    "--ram") set -- "$@" "-r" ;;
    "--bam_dir") set -- "$@" "-b" ;;
    "--out_dir") set -- "$@" "-o" ;;
    "--project_name") set -- "$@" "-n" ;;
    "--animal_lookup_path") set -- "$@" "-a" ;;
    *)        set -- "$@" "$arg"
  esac
done


while getopts c:e:f:t:r:b:o:n:a: opt ; do
   case $opt in
      c) config_dir=$OPTARG ;;
	    e) cp_dir=$OPTARG ;;
	    f) fastq_dir=$OPTARG ;;
	    t) threads=$OPTARG ;;
      r) ram=$OPTARG ;;
	    b) bam_dir=$OPTARG ;;
	    o) out_dir=$OPTARG ;;
	    n) project_name=$OPTARG ;;
      a) animal_lookup_path=$OPTARG ;;
      *) usage; exit 1;;
   esac
done


config_dir

# Launch the ipd ref matrix creator
python3 /Users/dabaker3/github/iwes_genotyping_v2/create_ipd_ref_matrix.py \
--config_dir=${config_dir}
# Align files
python3 /Users/dabaker3/github/iwes_genotyping_v2/semiperfect_align.py \
--cp_dir=${cp_dir} \
--fastq_dir=${fastq_dir} \
--bam_dir=${bam_dir} \
--config_dir=${config_dir} \
--threads=${threads} \
--ram=${ram}


# Expand the alignment and create depth of coverage plots
python3 /Users/dabaker3/github/iwes_genotyping_v2/main.py \
--project_name=${project_name} \
--out_dir=${out_dir} \
--config_dir=${config_dir} \
--bam_dir=${bam_dir}


# create a pivot table
python3 /Users/dabaker3/github/iwes_genotyping_v2/create_pivot_table.py \
--project_name=${project_name} \
--out_dir=${out_dir} \
--config_dir=${config_dir} \
--animal_lookup_path=${animal_lookup_path}