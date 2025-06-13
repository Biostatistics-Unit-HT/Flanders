include { INPUT_COLUMNS_VALIDATION } from "./modules/local/input_columns_validation"
include { RUN_MUNGING } from "./workflows/munging"
include { RUN_FINEMAPPING } from "./workflows/finemap"
include { RUN_COLOCALIZATION } from "./workflows/coloc"
include { PROCESS_BFILE } from "./modules/local/process_bfile"
include { GET_STUDY_FROM_ANNDATA } from "./modules/local/get_study_from_anndata"
include { completionSummary } from "./modules/local/pipeline_utils"
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

workflow {
	// validate and log pipeline parameters
	validateParameters()
	log.info paramsSummaryLog(workflow)

	// Define asset files
  	chain_file = file("${projectDir}/assets/hg19ToHg38.over.chain")
  	outdir_abspath = file(params.outdir).toAbsolutePath().toString()
	
	// In case we are running a test profile, we need to set the base_dir to the projectDir
	base_dir = params.is_test_profile ? "${projectDir}" : "${launchDir}"
	
	// Initialize empty channels for finemapping results
	credible_sets_from_finemapping = Channel.empty()
	credible_sets_from_input = Channel.empty()
	
	// Initialize empty channel for data from previous studies or coloc exclusions
	previous_h5ad_studies = Channel.empty()
	excluded_studies_from_input = Channel.empty()

	// --- COLOCALIZATION ---
	if (params.summarystats_input) {
		sumstats_input_file = file(params.summarystats_input, checkIfExists:true)

		// Use nf-schema to read and validate the sample sheet
		samplesheetToList(params.summarystats_input, 'assets/summarystats_input_schema.json')

		// Validate input file
		INPUT_COLUMNS_VALIDATION(sumstats_input_file, base_dir)
		
		// Collect and process distinct bim datasets
		INPUT_COLUMNS_VALIDATION.out.table_out
			.splitCsv(header:true, sep:"\t")
			.map{ row -> 
				def bfile_dataset = params.is_test_profile ? file("${projectDir}/${row.bfile}.{bed,bim,fam}") : file("${row.bfile}.{bed,bim,fam}")
				tuple(
					row.process_bfile,
					row.bfile,
					"${row.grch_bfile ? row.grch_bfile : row.grch}",
					"${params.run_liftover ? "T" : "F"}",
					bfile_dataset
				)
			}
			.unique()
			.branch { process_bfile_flag, bfile_id, grch_bfile, run_liftover, bfile_dataset ->
				need_processing: process_bfile_flag in ["T", "t", "TRUE", "true", "True"] || (grch_bfile == "37" && run_liftover == "T") || process_bfile_flag == null
				processed: true 
			}
			.set { bfile_datasets }

		PROCESS_BFILE(bfile_datasets.need_processing, chain_file)

		processed_bfile_datasets = bfile_datasets.processed
			.map { process_bfile_flag, bfile_id, grch_bfile, run_liftover, bfile_dataset -> 
				tuple(bfile_id, bfile_dataset)
			}
			.mix(PROCESS_BFILE.out.processed_dataset)

		// Generate a channel with finemapping configuration
		finemapping_config = INPUT_COLUMNS_VALIDATION.out.table_out
		.splitCsv(header:true, sep:"\t")
		.map{ row -> 
			tuple(
				row.bfile,
				[
				"study_id": row.study_id
				],
				[  
				"skip_dentist": params.skip_dentist,
				"maf": row.maf,
				"hole": row.hole,
				"cs_thresh": row.cs_thresh
				]
			)
		}
		.combine(processed_bfile_datasets, by: 0)
		.map { bfile_id, study_id, finemap_config, bfile_dataset ->
			tuple(study_id, finemap_config, bfile_dataset)
		}

		// Define input channel for munging of GWAS sum stats
		sumstas_input_ch = INPUT_COLUMNS_VALIDATION.out.table_out
		.splitCsv(header:true, sep:"\t")
		.map { row -> 
			def bfile_string = params.is_test_profile ? "${projectDir}/${row.bfile}" : "${row.bfile}"
			def gwas_file = params.is_test_profile ? file("${projectDir}/${row.input}", checkIfExists:true) : file("${row.input}", checkIfExists:true)
			def sdY_string = file(row.sdY).exists() ? file(row.sdY).name : row.sdY
			def sdY_file = file(row.sdY).exists() ? file(row.sdY) : file('NO_SDY_FILE')
			tuple(
				row.bfile,
				[
				"study_id": row.study_id
				],
				[
				"is_molQTL": row.is_molQTL,
				"run_liftover": params.run_liftover ? "T" : "F",
				"key": row.key,
				"chr_lab": row.chr_lab,
				"pos_lab": row.pos_lab,
				"rsid_lab": row.rsid_lab,
				"a1_lab": row.a1_lab,
				"a0_lab": row.a0_lab,
				"freq_lab": row.freq_lab,
				"n_lab": row.n_lab,
				"effect_lab": row.effect_lab,
				"se_lab": row.se_lab,
				"pvalue_lab": row.pvalue_lab,
				"type": row.type,
				"sdY": sdY_string,
				"s": row.s,
				"grch": row.grch,
				"maf": row.maf,
				"p_thresh1": row.p_thresh1,
				"p_thresh2": row.p_thresh2,
				"hole": row.hole
				],
				gwas_file,
				sdY_file
			)
		}
		.combine(processed_bfile_datasets, by: 0)
		.map { bfile_id, study_id, munging_config, gwas_file, sdY_file, bfile_dataset ->
			tuple(study_id, munging_config, gwas_file, sdY_file, bfile_dataset)
		}

		RUN_MUNGING(sumstas_input_ch, chain_file)
		RUN_FINEMAPPING(
			finemapping_config, 
			RUN_MUNGING.out.finemapped_loci, 
			RUN_MUNGING.out.munged_stats, 
			outdir_abspath
			)
		credible_sets_from_finemapping = RUN_FINEMAPPING.out.finemap_anndata
	}


	// --- COLOCALIZATION ---
	if (params.coloc_h5ad_input) {
		credible_sets_from_input = Channel.fromPath(params.coloc_h5ad_input, checkIfExists:true)
		if (params.coloc_skip_previous_studies) {
			GET_STUDY_FROM_ANNDATA(credible_sets_from_input)
			previous_h5ad_studies = GET_STUDY_FROM_ANNDATA.out.previous_study_ids
		}	
		full_credible_sets = credible_sets_from_input
			.mix(credible_sets_from_finemapping)
			.collect()
	} else {
		full_credible_sets = credible_sets_from_finemapping.collect()
		previous_h5ad_studies = Channel.value(file('NO_FILE'))
	}

	exclude_studies_file = params.coloc_exclude_studies_table ? 
        Channel.fromPath(params.coloc_exclude_studies_table, checkIfExists: true) : 
        Channel.value(file('NO_FILE'))

	if (params.run_colocalization || params.coloc_h5ad_input) {		
		RUN_COLOCALIZATION( full_credible_sets, previous_h5ad_studies, exclude_studies_file )
	}

	// At the end store params in yml and input files
	file("${params.outdir}/pipeline_inputs").mkdirs()
	Channel
		.fromList(params.entrySet())
		.map { entry -> "${entry.key}: ${entry.value}" }
		.collectFile(name: 'params.yml', storeDir: "${params.outdir}/pipeline_inputs", newLine: true)
		
		if (params.summarystats_input) {
			file(params.summarystats_input).copyTo("${params.outdir}/pipeline_inputs/summarystats_input.tsv")
		}
		if (params.coloc_h5ad_input) {
			file(params.coloc_h5ad_input).collectFile(
				newLine: false,
				name: "input_h5ad_inputs.txt",
				storeDir: "${params.outdir}/pipeline_inputs"
			)
		}

	workflow.onComplete {
		// At the end store log status
		completionSummary()
	}
}