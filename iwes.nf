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
	CREATE_REF_FASTAS ( )
	
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
if( params.animal.toLowerCase() == "mamu" ){
	params.legacy_fasta = params.mamu_legacy_fasta
} else {
	params.legacy_fasta = params.mafa_legacy_fasta
}
params.ipd_fasta = params.config_dir + "/" + "ipd_raw.fasta"
params.exon2_fasta = params.config_dir + "/" + "exon2_raw.fasta"
params.haplotype_lookup = params.config_dir + "/" + "haplotype_lookup.json"
params.haplotype_lookup_csv = params.config_dir + "/" + "haplotype_lookup.csv"
params.run_animal_lookup = params.config_dir + "/" + "baylor_33_mamu_lookup.csv"
params.miseq_to_ipd_lookup = params.config_dir + "/" + "miseq_to_ipd_lookup.json"
params.exon_to_ipd_json = params.config_dir + "/" + "exon_to_gen_lookup.json"

if( params.bam_dir.isEmpty() ){
	params.bam_dir = params.results + "/" + "01-" + params.run_name + "-alignments"
}
params.genotypes = params.results + "/" + "02-" + params.run_name + "-genotypes"
params.pivot_tables = params.results + "/" + "03-" + params.run_name + "-pivot_table"

// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

process CREATE_REF_FASTAS {
	
	// Here we produce folders of animal-specific reference "keys" or "lookups." These 
	// files tell the code downstream which of our reference sequences draw from the same 
	// alleles. This is a necessity with our current approach to genotyping, which uses all 
	// available gDNA sequences from the Immuno Polymorphism Database (IPD) plus exon 2 
	// "fall-backs." Some alleles in IPD are cDNA-only, with no introns. For those sequences, 
	// we fall back to either 1) an exon 2-only sequence or 2) a legacy ~150-bp genotyping 
	// amplicon from the MHC genotyping workflow at the AIDS Vaccine Research Laboratory (AVRL). 
	// In this step, the workflow uses mapping to identify which exon-2-only sequences, which 
	// AVRL amplicon sequences, and which IPD alleles correspond with each other. That way, 
	// the workflow is able to prioritize the best match to the longest available reference 
	// sequence downstream.
	
	output:
	val "cue", emit: cue
	
	script:
	ref_matrix_dir = file('${params.config_dir}/${params.animal)_ref_matrix')
	if( !ref_matrix_dir.exists() )
		ref_matrix_dir.mkdir()
	"""
	create_ref_fasta_lookups.py \
	--config_dir=${params.config_dir} \
	--miseq_legacy_db_path=${params.legacy_fasta} \
	--gen_db_path=${params.ipd_fasta} \
	--exon_db_path=${params.exon2_fasta} \
	--haplotype_json_path=${params.haplotype_lookup} \
	--species=${params.animal} \
	--cp_path=/miniconda2/bin/bbmap/current \
	--threads=${task.cpus} \
	--ram=${task.memory}
	"""
}

process SEMIPERFECT_ALIGN {
	
	// Here the workflow maps reads from each animal to all possible reference allele sequences. 
	// As described above, these reference allele sequences include full-length IPD gDNA sequences, 
	// exon-2-only sequences, and legacy AVRL amplicon sequences. Importantly, this step uses BBMap 
	// in "semiperfect mode," which allows that a read must match perfectly to the reference allele, 
	// but allows some overhang, e.g., the read hangs over the 5' or the 3' end of the reference 
	// sequences. This allows us to maintain stringent read-mapping while retaining as much data as 
	// possible.
	
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
	
	// In this step, the workflow identifies the best-possible allele sequence matches and uses 
	// pysam to remove a number of issues that may appear in the data. These include what we call 
	// "computational chimeras;" these are cases where a read from one allele maps perfectly to a 
	// different allele. With short-read data like that used here, this happens often in the MHC 
	// region, where many alleles are extremely similar to one another, save for a few key differences 
	// separated by hundreds of bases. To account for this, our script requires that all reads must 
	// start within 50 bases of one another. It also requires that paired reads must be within 500 
	// bases either end of the reference, as otherwise, it's likely that the two reads were sequenced 
	// from different alleles. In that case, we keep one of the two paired alleles. Finally, we 
	// require that each allele has a certain baseline read depth-of-coverage with the parameter 
	// "depth_threshold," though we make this less stringent near the edges of the reference alleles 
	// with the parameter "edge_distance_threshold." All of these parameters can be controlled in the 
	// file nextflow.config. We recomment that users do not modify the sources code in this .nf script 
	// or in all the bin/ scripts. This process runs only for Rhesus macaques.
	
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
	
	// This process does all the same thomgs as GENOTYPE_MAMU, except for Cynomolgus Macaque reads and
	// reference allele sequences.
	
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
	
	// Now that mapping, genotyping, and filtering is complete, this step simply formats all the data
	// in a convenient, Pivot table excel file. This file indicates cases where reads mapped ambiguously
	// between two alleles, when alleles have untrustworthy depth-of-coverage, when reads mapped to a
	// shorter reference allele instead of a full-length gDNA with introns, etc.
	
	tag "${params.animal}"
	publishDir params.pivot_tables, mode: 'copy'
	
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
	--haplotype_lookup=${params.haplotype_lookup_csv} \
	--diag_to_ipd_json=${params.miseq_to_ipd_lookup} \
	--exon_to_ipd_json=${params.exon_to_ipd_json}
	"""
}

// --------------------------------------------------------------- //
