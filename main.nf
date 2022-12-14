#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {

	// Disclaimer
	println "NOTE:"
	println "THIS WORKFLOW CURRENTLY ONLY SUPPORTS RHESUS MACAQUE DATA"
	
	// input channels
	ch_reads = Channel
		.fromFilePairs ( '${params.fastq_dir}/*_R{1,2}_001.fastq.gz', flat: true )
	
	// ch_samples = Channel
	// 	.fromPath( params.samplesheet )
	// 	.splitCsv( )
	// 	.map { row -> tuple( row.accession, row.animal, File(row.reads1), File(row.reads2) )  }
	
	ch_alignments = Channel
		.fromPath ( "${params.bam_dir}/*.bam" )
	
	// Workflow steps 
	CREATE_REF_MATRIX ( )
	
	SEMIPERFECT_ALIGN ( 
		CREATE_REF_MATRIX.out.cue,
		ch_reads
	)
	
	GENOTYPE_MAMU ( 
		CREATE_REF_MATRIX.out.cue,
		SEMIPERFECT_ALIGN.out
			.mix ( ch_alignments )
			.collect()
	)

	// GENOTYPE_MAFA ( 
	// 	SEMIPERFECT_ALIGN.out.collect()
	// )
	
	CREATE_PIVOT_TABLE ( 
		GENOTYPE.out
	)
	
	
}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config
params.config_dir = params.resources + "/" + params.animal
params.bait_fasta = params.config_dir + "/" + "bait.fasta"
params.run_animal_lookup = params.config_dir + "/" + "baylor_33_mamu_lookup.csv"
params.haplotype_lookup = params.config_dir + "/" + "haplotype_lookup.csv"
params.ipd_avrl_dict = params.config_dir + "/" + "ipd_to_diag_lookup.json"
params.ipd_num_lookup = params.config_dir + "/" + "ipd_num_lookup.json"

if( params.bam_dir.isEmpty() ){
	params.bam_dir = params.results + "/" + "01-" + params.run_name + "-alignments"
}
params.genotypes = params.results + "/" + "02-" + params.run_name + "-genotypes"
params.pivot_tables = params.results + "/" + "03-" + params.run_name + "-pivot_table"

// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

process CREATE_REF_MATRIX {
	
	// This process does something described here
	
	output:
	val "cue", emit: cue
	
	script:
	ref_matrix_dir = file('${params.config_dir}/${params.animal)_ref_matrix')
	if( !ref_matrix_dir.exists() )
		ref_matrix_dir.mkdir()
	"""
	create_ref_fasta_lookups.py \
	--config_dir=${params.config_dir} \
	--miseq_legacy_db_path= \
	--gen_db_path= \
	--exon_db_path= \
	--haplotype_json_path= \
	--species=${params.animal} \
	--cp_path=/miniconda2/bin/bbmap/current \
	--threads=${task.cpus} \
	--ram=${task.memory}
	"""
}

process SEMIPERFECT_ALIGN {
	
	// This process does something described here
	
	tag "${accession}"
	publishDir params.bam_dir, mode: 'copy'
	
	memory params.ram
	cpus 4
	
	input:
	each val(cue)
	tuple val(accession), path(reads1), path(reads2)
	
	output:
	path "*.bam"
	
	script:
	"""
	semiperfect_align.py \
	--cp_dir=/miniconda2/bin/bbmap/current \
	--fastq_dir=. \
	--bam_dir=. \
	--bait_fasta=${params.bait_fasta} \
	--config_dir=${params.config_dir} \
	--threads=${task.cpus} \
	--ram=${task.memory}
	"""
}

process GENOTYPE_MAMU {
	
	// This process does something described here
	
	tag "${params.animal}"
	publishDir params.genotypes, mode: 'copy'
	
	input:
	each val(cue)
	path bam_list
	
	output:
	tuple path("*.csv"), val(params.animal)

	when:
	params.animal.toLowerCase() == "mamu"
	
	script:
	"""
	genotype.py \
	--out_dir=. \
	--bam_dir=. \
	--project_name=${params.run_name} \
	--config_dir=${params.config_dir} \
	--bait_fasta=${params.bait_fasta} \
	--ipd_ref_matrix_dir='${params.config_dir}/${params.animal}_ref_matrix' \
	--unpaired_edge_threshold=${params.unpaired_edge_threshold} \
	--depth_threshold=${params.depth_threshold} \
	--edge_distance_threshold=${params.edge_distance_threshold} \
	--low_ram_mode=${params.genotype_in_low_ram_mode}
	"""
}

process GENOTYPE_MAFA {
	
	// This process does something described here
	
	tag "${params.animal}"
	publishDir params.genotypes, mode: 'copy'
	
	input:
	each val(cue)
	path bam_list
	
	output:
	tuple path("*.csv"), val(params.animal)

	when:
	params.animal.toLowerCase() == "mafa"
	
	script:
	"""
	genotype.py \
	--out_dir=. \
	--bam_dir=. \
	--project_name=${params.run_name} \
	--config_dir=${params.config_dir} \
	--bait_fasta=${params.bait_fasta} \
	--ipd_ref_matrix_dir='${params.config_dir}/${params.animal}_ref_matrix' \
	--unpaired_edge_threshold=${params.unpaired_edge_threshold} \
	--depth_threshold=${params.depth_threshold} \
	--edge_distance_threshold=${params.edge_distance_threshold} \
	--low_ram_mode=${params.genotype_in_low_ram_mode}
	"""
}

process CREATE_PIVOT_TABLE {
	
	// This process does something described here
	
	tag "${params.animal}"
	publishDir params.results, mode: 'copy'
	
	input:
	tuple path(csv_files), val(animal)
	
	output:
	path ".xlsx"
	
	script:
	"""
	create_pivot_table.py \
	--out_dir=./ \
	--project_name=${params.run_name} \
	--config_dir=${params.config_dir} \
	--animal_lookup_path=${params.run_animal_lookup} \
	--bait_fasta=${params.bait_fasta} \
	--haplotype_lookup=${params.haplotype_lookup} \
	--diag_to_ipd_json=${params.ipd_avrl_dict} \
	--exon_to_ipd_json=
	"""
}

// --------------------------------------------------------------- //
