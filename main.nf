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
		CONVERT_BAM.out.accessions
			.unique()
			.collectFile(name: 'accessions.txt', newLine: true)
	)
	
	SAVE_REF_WITH_BAM ( 
		CONVERT_BAM.out.bams.collect()
	)
	
	COMPUTE_DEPTH ( 
		MERGE_BAM.out.merged_bam
	)
	
	COMPRESS_DEPTH ( 
		COMPUTE_DEPTH.out
	)
	
	FILTER_DEPTH_CHIMERAS ( 
		COMPRESS_DEPTH.out,
		MERGE_BAM.out
	)
	
	CREATE_ALLELE_LIST_FULL_LENGTH ( 
		FILTER_DEPTH_CHIMERAS.out.allele_list
	)
	
	GENOTYPE ( 
		FILTER_DEPTH_CHIMERAS.out.allele_list
	)
	
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
params.allele_lists = params.results + "/06-allele_lists"
params.genotypes = params.results + "/07-full_length_genotypes"
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
	path "*.read_count.txt" 
	
	when:
	params.use_avrl_amplicons == true
	
	script:
	"""
	java -ea -Xmx${task.memory}m -Xms${task.memory}m -cp ${params.bbmap_cp} align2.BBMap build=1 \
	in=${reads1} in2=${reads2} ref=${params.bait_fasta} outm=${accession}.fastq.gz \
	semiperfectmode=t threads=${task.cpus} nodisk=t
	
	zcat ${accession}.fastq.gz | echo $((`wc -l`/4)) > ${accession}.read_count.txt
	
	touch ${accession}.read_count.txt
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
	java -ea -Xmx${task.memory}m -Xms${task.memory}m -cp ${params.bbmap_cp} align2.BBMap build=1 \
	in=${baited_fastq} ref=${ref_allele} outm=${accession}-${allele}.sam \
	semiperfectmode=t threads=${task.memory} int=t nodisk=t >/dev/null=None 2>&1=None
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
	val accession, emit: accessions
	
	when:
	ref_allele.getSimpleName() == sam.getSimpleName().split("-")[1]
	
	script:
	mapping = sam.getSimpleName()
	accession = sam.getSimpleName().split("-")[0]
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
	path accession_list
	
	output:
	tuple path("*.merged.bam"), emit: merged_bam
	path "*.bam.bai", emit: merged_index
	
	shell:
	"""
	for i in `cat !{accession_list}`;
	do
		find . -name '${i}*.bam' > ${i}_bam_list.txt
		mkdir -p ${i}_bam_split/
		split -l 200 ${i}_bam_list.txt ${i}_bam_split/bam_segment
		find ${i}_bam_split/ -name 'bam_segment*' > ${i}_split_files.txt
		for f in `cat ${i}_split_files.txt`;
		do
			samtools merge ${f}.bam -b ${f} && samtools index ${f}.bam
		done
		
		find ${i}_bam_split/ -name '*.bam' > ${i}_merged_bam_files.txt
		samtools merge ${i}.merged.bam -b ${i}_merged_bam_files.txt && samtools index ${i}.merged.bam
		rm -rf ${i}_bam_split/
	done
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
	
	input:
	each path(depth_stats)
	each path(merged_bam)
	
	output:
	path "*.allele_list.tsv", emit: allele_list
	tuple path("*.filtered.merged.bam"), val(accession), emit: filtered_bam
	
	when:
	depth_stats.getSimpleName() == merged_bam.getSimpleName()
	
	script:
	accession = merged_bam.getSimpleName()
	"""
	filter_depth_chimeras.py \
	--depth_input_path=${depth_stats} \
	--merged_bam_path=${merged_bam} \
	--filtered_allele_list_outpath="${accession}.allele_list.tsv" \
	--filtered_merged_bam_outpath=="${accession}.filtered.merged.bam" \
	--edge_distance_threshold=${params.edge_distance_threshold} \
	--depth_threshold=${params.depth_threshold} \
	--maximum_start_position_gap=${params.maximum_start_position_gap} \
	--minimum_bridge_read_length=${params.minimum_bridge_read_length} \
	"""
}

process CREATE_ALLELE_LIST_FULL_LENGTH {
	
	// exhautively map baited reads by mapping one at a time
	// to reference sequences that have exons 2-4
	// using bbmap
	// suppress output stderr and stdout because it consumes a lot of unnecessary space
	
	tag "${tag}"
	publishDir params.allele_lists, mode: 'copy'
	
	cpus 1
	
	input:
	path allele_list
	
	output:
	path "*.allele_list_fl.txt"
	path "*.allele_list_fl_num.txt"
	
	when:
	params.use_avrl_amplicons == true
	
	script:
	"""
	create_allele_list_full_length.py \
	${params.ipd_avrl_dict} \
	${ipd_num_lookup} \
	${allele_list} \
	${params.missing_alleles}
	"""
}

process GENOTYPE {
	
	// identify SAM mappings where there is complete read support at 
	// least specified depth
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	cpus 1
	
	input:
	path allele_list
	
	output:
	path "*.genotypes.csv"
	
	when:
	params.use_avrl_amplicons == false
	
	script:
	"""
	genotype.py ${accession} ${params.ipd_avrl_dict}
	"""
}


// --------------------------------------------------------------- //
