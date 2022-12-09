#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {
	
	
	// input channels
	ch_reads = Channel
		.fromFilePairs ( '${params.data_dir}/*_R{1,2}_001.fastq.gz', flat: true )
	
	
	// Workflow steps 
	SEMIPERFECT_ALIGN ( )
	
	GENOTYPE ( )
	
	CREATE_PIVOT_TABLE ( )
	
	
}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config
if( params.bam_dir.isEmpty() ){
	params.bam_dir = params.results + "/" + "01-" + params.run_name + "-alignments"
}
params.genotypes = params.results + "/" + "02-" + params.run_name + "-genotypes"
params.pivot_tables = params.results + "/" + "02-" + params.run_name + "-pivot_table"

// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

process SEMIPERFECT_ALIGN {
	
	// This process does something described here
	
	// tag "${tag}"
	publishDir params.bam_dir, mode: 'copy'
	
	memory params.ram
	cpus 4
	
	input:
	tuple path(reads1), path(reads2)
	
	output:
	path "*.bam"
	
	when:
	params.bam_dir.isEmpty()
	
	script:
	"""
	semiperfect_align.py \
	--cp_dir=/Users/dabaker3/anaconda3/bin/bbmap/current \
	--fastq_dir=. \
	--bam_dir=. \
	--config_dir=${params.resources} \
	--threads=${task.cpus} \
	--ram=${task.memory}
	"""
}

process GENOTYPE {
	
	// This process does something described here
	
	// tag "${tag}"
	publishDir params.genotypes, mode: 'copy'
	
	output:
	path "*"
	
	when:
	!params.bam_dir.isEmpty()
	
	script:
	"""
	genotype.py \
	--project_name=${params.run_name} \
	--out_dir=. \
	--config_dir=${params.resources} \
	--bam_dir=.
	"""
}

process CREATE_PIVOT_TABLE {
	
	// This process does something described here
	
	// tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	memory 1.GB
	cpus 1
	time '10minutes'
	
	input:
	path ""
	
	output:
	
	
	when:
	
	
	script:
	"""
	create_pivot_table.py \
	--project_name=${params.run_name} \
	--out_dir=./ \
	--config_dir=/Users/dabaker3/github/iwes_genotyping_v2/config/Mamu \
	--animal_lookup_path=/Users/dabaker3/github/iwes_genotyping_v2/config/Mamu/baylor_32_animal_lookup.csv

	"""
}

// --------------------------------------------------------------- //
