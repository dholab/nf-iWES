#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {
	
	
	// input channels
	ch_reads = Channel
		.fromFilePairs ( '${params.data_dir}/*_R{1,2}_001.fastq.gz', flat: true )
	
	ch_ref_alleles = Channel
		.fromPath ( '${params.avrl_amplicon_path}/*.fasta' )
	
	
	// Workflow steps 
	BAIT_MHC ( 
		ch_reads
	)
	
	EXHAUSTIVE_MAPPING ( 
		BAIT_MHC.out.baited_fastq,
		ch_ref_alleles
	)
	
	SORT_SAM ( 
		EXHAUSTIVE_MAPPING.out
	)
	
	INDEX_FASTA ( 
		ch_ref_alleles
	)
	
	CONVERT_BAM ( 
		SORT_SAM.out,
		ch_ref_alleles
	)
	
	MERGE_BAM ( 
		CONVERT_BAM.out.bams.collect()
	)
	
	SAVE_REF_WITH_BAM ( 
		CONVERT_BAM.out.bams.collect()
	)
	
	COMPUTE_DEPTH ( 
		MERGE_BAM.out
	)
	
	COMPRESS_DEPTH ( 
		COMPUTE_DEPTH.out
	)
	
	FILTER_DEPTH_CHIMERAS ( 
		COMPRESS_DEPTH.out,
		MERGE_BAM.out
	)
	
	CREATE_ALLELE_LIST_FULL_LENGTH ( )
	
	GENOTYPE ( )
	
}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config
params.baiting_output = params.results + "/01-bait-mhc"
params.exhaustive_mapping_output = params.results + "/02-exhaustive-map"
params.converted = params.results + "/03-converted_to_bam"
params.merged_bams = params.results + "/04-merged"
params.depth_stats = params.results + "/05-depth"
// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

process BAIT_MHC {
	
	// baiting extracts MHC reads from iWES FASTQ
	// this makes it faster to do exhaustive searching for MHC-mapped reads that map perfectly
	// it also simplifies troubleshooting and optimizing genotyping parameters
	// download FASTQ from SRA with fasterq-dump and send to stdout, then
	// map FASTQ to CY0333 gDNA MHC reads
	// save mapped reads as FASTQ file
	
	tag "${accession}"
	
	publishDir params.baiting_output, mode: 'copy', pattern: '*.fastq.gz'
	publishDir params.baiting_output, mode: 'copy', pattern: '*read_count.txt'
	
	// memory
	// cpus 
	
	input:
	tuple val(accession), path(reads1), path(reads2)
	
	output:
	path "*.fastq.gz", emit: baited_fastq
	path("*.read_count.txt")
	
	when:
	params.use_avrl_amplicons == true
	
	script:
	"""
	#!/usr/bin/env python3
	
	import subprocess
	
	main_cmd = 'java -ea -Xmx{0}m -Xms{0}m -cp {1} align2.BBMap build=1'.format(${task.memory}, params.bbmap_cp)
	
	required_arg_list = ['in={0}'.format(${reads1}),
						 'in2={0}'.format(${reads2}),
						 'ref={0}'.format(${params.bait_fasta}),
						 'outm={0}'.format(${accession}.fastq.gz)]
	default_arg_dict = {'semiperfectmode': 't',
						'threads': ${task.cpus},
						'nodisk':'t'}
	
	run_cmd = create_cmd_req_list_optional_dict(main_cmd=main_cmd,
		required_arg_list=required_arg_list,
		default_arg_dict=default_arg_dict,
		optional_arg_dict=params.optional_arg_dict,
		delimiter='=')
	subprocess.call(run_cmd,shell=True)
	shell("zcat ${accession}.fastq.gz | echo $((`wc -l`/4)) > ${accession}.read_count.txt")
	shell("touch ${accession}.read_count.txt")
	"""
}

process EXHAUSTIVE_MAPPING {
	
	// exhautively map baited reads by mapping one at a time
	// to reference sequences that have exons 2-4
	// using bbmap
	// suppress output stderr and stdout because it consumes a lot of unnecessary space
	// 1 thread
	
	tag "${accession}-${allele}"
	// publishDir params.exhaustive_mapping_output

	// memory 1.GB
	cpus 1
	
	input:
	each path(baited_fastq)
	each path(ref_allele)
	
	output:
	path "*.sam"
	
	script:
	allele = ref_allele.getSimpleName()
	accession = baited_fastq.getSimpleName()
	"""
	#!/usr/bin/env python3
	
	main_cmd = 'java -ea -Xmx{0}m -Xms{0}m -cp {1} align2.BBMap build=1'.format(${task.memory}, ${params.bbmap_cp})
	
	required_arg_list = ['in={0}'.format(${baited_fastq}),
						 'ref={0}'.format(${ref_allele}),
						 'outm={0}'.format(${accession}-${allele}.sam)]
	default_arg_dict = {'semiperfectmode': 't',
						'threads': str(threads),
						'int': 't',
						'nodisk': 't',
						'>/dev/null': None,
						'2>&1': None}
	
	run_cmd = create_cmd_req_list_optional_dict(main_cmd=main_cmd,
		required_arg_list=required_arg_list,
		default_arg_dict=default_arg_dict,
		delimiter='=')
	subprocess.call(run_cmd,shell=True)
	"""
}

process SORT_SAM {
	
	// sort SAM files following bbmap
	// make list of output files to use with samtools merge
	
	tag "${mapping}"
	publishDir params.exhaustive_mapping_output, mode: 'symlink'
	
	cpus 1
	
	input:
	path sam
	
	output:
	path "*.sam"
	
	script:
	mapping = sam.getSimpleName()
	"""
	samtools sort ${sam} -o ${mapping}.sorted.sam
	"""
}

process INDEX_FASTA {
	
	// create index for each individual FASTA file
	// sometimes fails for no reason, so retry three times
	
	tag "${allele}"
	publishDir params.avrl_amplicon_path, mode: 'copy'
	
	cpus 1
	errorStrategy 'retry'
	maxRetries 3
	
	input:
	path fasta
	
	output:
	path "*.fasta.fai"
	
	script:
	allele = fasta.getSimpleName()
	"""
	samtools faidx ${fasta}
	touch ${allele}.fasta.fai
	"""
}

process CONVERT_BAM {
	
	// convert SAM file to BAM
	
	tag "${mapping}"
	publishDir params.converted, mode: 'copy'
	
	input:
	path sam
	each path(ref_allele)
	
	output:
	path "*.bam", emit: bams
	path "*.bam.bai", emit: bam_indices
	
	when:
	ref_allele.getSimpleName() == sam.getSimpleName().split("-")[1]
	
	script:
	mapping = sam.getSimpleName()
	allele = sam.getSimpleName().split("-")[1]
	"""
	samtools view -b -h -T {input.indv_fasta} -o ${mapping}.bam ${sam} \
	&& samtools index ${mapping}.bam
	"""
}

process MERGE_BAM {
	
	// merge sorted SAM files into single SAM file
	// Some nodes/systems have file limits set.  
	// By looping over 200 files at a time we should not meet that file open limit.
	// The default file limit to be around 1024, but there is a lot opened in the background 
	// (and possibly on the shared underlying node) that can cause this to be 
	// reached as low as 500 files.
	
	publishDir params.merged_bams, mode: 'copy'
	
	input:
	path bam_list
	
	output:
	path "*.merged.bam"
	
	script:
	"""
	#!/usr/bin/env python3
	
	import os
	merged_bam = os.path.join(${params.results},'04_merged_bam_files.txt')
	split_files = os.path.join(${params.results},'04_split_files.txt')
	converted_dir = os.path.join(${params.results},'04-converted')
	filelist = os.path.join(${params.results},'04_filelist.txt')
	bam_split_dir = os.path.join(${params.results},'04_bam_split')
	bam_segment = os.path.join(bam_split_dir,'bam_segment')
	shell("find {converted_dir}  -name '*.bam' > {filelist}")
	shell("mkdir -p {bam_split_dir}")
	shell("split -l 200 {filelist} {bam_segment}")
	shell("find {bam_split_dir} -name 'bam_segment*' > {split_files}")
	shell("for f in `cat {split_files}`; do \
	samtools merge ${{f}}.bam -b ${{f}} && samtools index ${{f}}.bam; \
	done")
	shell("find {bam_split_dir} -name '*.bam' > {merged_bam}")
	shell("samtools merge {output.out_bam} -b {merged_bam} && samtools index {output.out_bam}")
	shell("rm -rf {bam_split_dir}")
	"""
}

process SAVE_REF_WITH_BAM {
	
	// save the reference sequence with the BAM file for visual inspection
	
	input:
	path bam_list
	
	script:
	"""
	cp ${params.full_ref_fasta} ${params.converted}
	"""
}

process COMPUTE_DEPTH {
	
	// use samtools depth to compute depth at each position of SAM file
	// by default several types of reads are filtered from depth calculations
	// add back with -g [UNMAP,SECONDARY,QCFAIL,DUP]
	
	tag "${accession}"
	
	cpus 1
	
	input:
	path merged_bam
	
	output:
	tuple path("*.txt"), val(accession)
	
	script:
	accession = merged_bam.getSimpleName()
	"""
	samtools depth -a ${merged_bam} -g UNMAP,SECONDARY,QCFAIL,DUP -o ${accession}.depth.txt
	"""
}

process COMPRESS_DEPTH {
	
	// ZIP compress depth file
	
	tag "${accession}"
	publishDir params.depth_stats, mode: 'copy'
	
	cpus 1
	
	input:
	tuple path(txt), val(accession)
	
	output:
	path "*.txt.gz"
	
	script:
	"""
	gzip -c ${txt} > ${accession}.txt.gz
	"""
}

process FILTER_DEPTH_CHIMERAS {
	
	// 1.) filters the depth of coverage (based on a depth and position from the edge)
	// 2.) next step filters out chimeras:
	// 2a.)It ensures there is at least one mapping read a set gap (default 30 positions) apart across the allele
	// 3.) Next a list of passing alleles is created as a single column csv file
	// 4.) Next a bam files with mapped reads to the the "passing" alleles 
	
	tag "${accession}"
	publishDir params.depth_stats, mode: 'copy', pattern: '*.allele_list.tsv'
	publishDir params.merged_bams, mode: 'copy', pattern: '*.filtered.merged.bam'
	
	cpus 1
	
	input:
	each path(depth_stats)
	each path(merged_bam)
	
	output:
	path "*.allele_list.tsv"
	tuple path("*.filtered.merged.bam"), val(accession)
	
	when:
	depth_stats.getSimpleName() == merged_bam.getSimpleName()
	
	script:
	accession = merged_bam.getSimpleName()
	"""
	#!/usr/bin/env python3
	
	import subprocess
	
	main_cmd = 'filter_depth_chimera.py'
	
	required_arg_list = ['--depth_input_path={0}'.format(${depth_stats}),
						 '--merged_bam_path={0}'.format(${merged_bam}),
						 '--filtered_allele_list_outpath={0}'.format("${accession}.allele_list.tsv"),
						 '--filtered_merged_bam_outpath={0}'.format("${accession}.filtered.merged.bam")]
	default_arg_dict = {'--edge_distance_threshold': '0',
						'--depth_threshold': '10',
						"--maximum_start_position_gap": '45',
						'--minimum_bridge_read_length': '70'}
	
	run_cmd = create_cmd_req_list_optional_dict(main_cmd=main_cmd,
		required_arg_list=required_arg_list,
		default_arg_dict=default_arg_dict,
		optional_arg_dict={},
		delimiter='=')
	subprocess.call(run_cmd,shell=True)
	"""
}

process CREATE_ALLELE_LIST_FULL_LENGTH {
	
	// exhautively map baited reads by mapping one at a time
	// to reference sequences that have exons 2-4
	// using bbmap
	// suppress output stderr and stdout because it consumes a lot of unnecessary space
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	cpus 1
	
	input:
	
	
	output:
	
	
	when:
	params.use_avrl_amplicons == true
	
	script:
	"""
	#!/usr/bin/env python3
	
	import pandas as pd
	import json
	import os
	import shutil
	from pathlib import Path


	# with open("${params.ipd_avrl_dict") as f_in:
	with open("${params.ipd_avrl_dict") as f_in:
		ipd_diag_dict = json.load(f_in)
	with open("${params.ipd_avrl_dict") as f_in:
		ipd_num_dict = json.load(f_in)

	sample = ${accession}
	# open the list of missing alleles do to no corresponding diag region to the fl.
	df_missing = pd.read_csv("${params.missing_alleles}",sep='\t',header=None,names=['allele'])
	missing_allele_list = list(df_missing['allele'].unique())
	# ipd_num_dict
	included = False
	# open the list of alleles that past the depth/computational chimera filter
	df = pd.read_csv("${accession}.allele_list.tsv"",sep='\t')
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
	# write teh two lists of alleles and corresponding number list for the fasta files.
	with open("${accession}.allele_list_fl.txt",'w') as f:
		f.write('\n'.join(ipd_allele_list))
	ipd_num_list = [ipd_num_dict[x] for x in ipd_allele_list]
	with open("${accession}.allele_list_fl_num.txt",'w') as f:
		f.write('\n'.join(ipd_num_list))

	shell("touch ${accession}_finished.txt")
	"""
}

process GENOTYPE {
	
	// identify SAM mappings where there is complete read support at 
	// least specified depth
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	memory 1.GB
	cpus 1
	time '10minutes'
	
	input:
	
	
	output:
	path "*.genotypes.csv"
	
	when:
	params.use_avrl_amplicons == false
	
	script:
	"""
	#!/usr/bin/env python3
	
	import os
	import shutil
	from pathlib import Path
	
	genotype = "${accession}.genotypes.csv"
	allele_diag_list = "${accesion}.avrl_allele_list.tsv"
	allele_fl_list = "${accesion}.allele_list.tsv"
	with open(${params.ipd_avrl_dict}) as f_in:
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
	df['accession'] = ${accession}
	if len(missing_diag_list) > 0 :
		df_diag_missing = pd.DataFrame({'allele':missing_diag_list})
		df_diag_missing = df_diag.merge(df_diag_missing, on=['allele'], how='inner')
		df_diag_missing.rename(columns={'depth':'read_ct'}, inplace=True)
		df_diag_missing['read_ct'] = df_diag_missing['read_ct'].round(decimals=0)
		df_diag_missing['read_ct'] = df_diag_missing['read_ct'] - .01
		df_diag_missing['accession'] = ${accession}
		df = pd.concat([df, df_diag_missing], ignore_index=True)
	df.to_csv(genotype, index=False)
	
	"""
}


// --------------------------------------------------------------- //
