include { INPUT_COLUMNS_VALIDATION } from "./modules/local/input_columns_validation"
include { RUN_MUNGING } from "./workflows/munging"
include { RUN_FINEMAPPING } from "./workflows/finemap"
include { RUN_COLOCALIZATION } from "./workflows/coloc"
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

	if (params.summarystats_input) {
		sumstats_input_file = file(params.summarystats_input, checkIfExists:true)

		// Use nf-schema to read and validate the sample sheet
		samplesheetToList(params.summarystats_input, 'assets/summarystats_input_schema.json')

		// Validate input file
		INPUT_COLUMNS_VALIDATION(sumstats_input_file, base_dir)

		finemapping_config = Channel  
		.of(sumstats_input_file)
		.splitCsv(header:true, sep:"\t")
		.map{ row -> 
			def bfile_string = params.is_test_profile ? "${projectDir}/${row.bfile}" : "${row.bfile}"
			tuple(
				[
				"study_id": row.study_id
				],
				[  
				"p_thresh3": row.p_thresh3,
				"p_thresh4": row.p_thresh4,
				"bfile": bfile_string,
				"skip_dentist": params.skip_dentist,
				"maf": row.maf,
				"hole": row.hole,
				"cs_thresh": row.cs_thresh
				]
			)
		}

		// Define input channel for munging of GWAS sum stats
		sumstas_input_ch = INPUT_COLUMNS_VALIDATION.out.table_out
		.splitCsv(header:true, sep:"\t")
		.map { row -> 
			def bfile_string = params.is_test_profile ? "${projectDir}/${row.bfile}" : "${row.bfile}"
			def gwas_file = params.is_test_profile ? file("${projectDir}/${row.input}", checkIfExists:true) : file("${row.input}", checkIfExists:true)
			tuple(
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
				"sdY": row.sdY,
				"s": row.s,
				"grch": row.grch,
				"bfile": bfile_string,
				"maf": row.maf,
				"p_thresh1": row.p_thresh1,
				"p_thresh2": row.p_thresh2,
				"hole": row.hole
				],
				gwas_file
			)
		}

		RUN_MUNGING(sumstas_input_ch, chain_file)
		RUN_FINEMAPPING(
			finemapping_config, 
			RUN_MUNGING.out.finemapped_loci, 
			RUN_MUNGING.out.munged_stats, 
			outdir_abspath
			)
		credible_sets_from_finemapping = RUN_FINEMAPPING.out.coloc_master.collectFile(
			newLine: false, 
			name: "all_credible_sets_from_finemapping.tsv", 
			storeDir: "${params.outdir}/results/finemap",
			keepHeader: true,
			skip: 1
		)
	}

	if (params.coloc_input) {
		credible_sets_from_input = Channel.fromPath(params.coloc_input, checkIfExists:true)

		// Use nf-schema to read and validate the sample sheet
		samplesheetToList(params.coloc_input, 'assets/coloc_input_schema.json')
	}

	if (params.run_colocalization || params.coloc_input) {
		file("${params.outdir}/results/coloc_info_tables").mkdirs()
		
		full_credible_sets = credible_sets_from_input
			.mix(credible_sets_from_finemapping)
			.collectFile(
				newLine: false, 
				name: "ALL_COMBINED_coloc_info_master_table.tsv", 
				storeDir: "${params.outdir}/results/coloc_info_tables",
				keepHeader: true,
				skip: 1
			)
		
		colocalization_input = full_credible_sets
			.splitCsv(header:true, sep:"\t")
			.map{ row ->
				[
				"chr_cs": row.chr
				]
			}
			.unique()
			.combine(full_credible_sets)

		RUN_COLOCALIZATION( colocalization_input )
	}
}
