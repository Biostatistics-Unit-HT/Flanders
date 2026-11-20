include { INPUT_COLUMNS_VALIDATION } from "./modules/local/input_columns_validation"
include { RUN_MUNGING } from "./workflows/munging"
include { RUN_FINEMAPPING } from "./workflows/finemap"
include { RUN_COLOCALIZATION } from "./workflows/coloc"
include { RUN_JOINT_FINEMAPPING } from "./workflows/joint_finemapping"
include { PROCESS_BFILE } from "./modules/local/process_bfile"
include { GET_STUDY_FROM_ANNDATA } from "./modules/local/get_study_from_anndata"
include { completionSummary } from "./modules/local/pipeline_utils"
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { parseIntervals } from './lib/utils'

// Internal default: DENTIST execution is currently unsupported.
params.skip_dentist = true

workflow {

	// Joint finemapping requires performing colocalization.
    if (params.run_joint_finemapping && !params.run_colocalization) {
        error "Joint finemapping requires colocalization. Set --run_colocalization true or --run_joint_finemapping false."
    }

	// validate and log pipeline parameters
	validateParameters()
	log.info paramsSummaryLog(workflow)

	// additional validation for chromosome intervals
	def parsed_intervals = parseIntervals(params.chromosomes)
	log.info "Parsed chromosome intervals: ${parsed_intervals}"

	// Define asset files
  	chain_file = file("${projectDir}/assets/hg19ToHg38.over.chain")
  	outdir_abspath = file(params.outdir).toAbsolutePath().toString()
	
	// In case we are running a test profile, we need to set the base_dir to the projectDir
	base_dir = params.is_test_profile ? "${projectDir}" : "${launchDir}"
	
	// Initialize empty channels for finemapping results
	credible_sets_from_finemapping = Channel.empty()
	credible_sets_from_input = Channel.empty()

	// Output a WARNING when coloc guide and exclusion filters are both provided
	// We are going to ignore everything and just apply the coloc guide table provided
	if (params.coloc_guide_table && (params.coloc_exclude_studies_table || params.coloc_skip_previous_studies)) {
		log.warn """
		=== ⚠️  WARNING  ⚠️ ===
		You provided an external coloc guide table (${params.coloc_guide_table}), but also a coloc_exclude_studies_table or coloc_skip_previous_studies is true.
		Only the external guide table will be used.
		"""
	} 

	// --- INPUT VALIDATION, MUNGING and FINEMAPPING ---
	if (params.summarystats_input) {
		sumstats_input_file = file(params.summarystats_input, checkIfExists:true)
	
		// samplesheetToList returns positional rows. This order must match the property order in assets/summarystats_input_schema.json.
		def schema_columns = [
			'input',
			'study_id',
			'chr_lab',
			'pos_lab',
			'rsid_lab',
			'a1_lab',
			'a0_lab',
			'freq_lab',
			'n_lab',
			'effect_lab',
			'se_lab',
			'pvalue_lab',
			'is_molQTL',
			'key',
			'type',
			'sdY',
			's',
			'grch',
			'bfile',
			'maf',
			'p_thresh1',
			'p_thresh2',
			'cs_thresh',
			'hole',
			'grch_bfile',
			'process_bfile'
		]

		// Resolve relative paths against the launch directory, or the pipeline directory for the test profile.
		def resolve_input_path = { value ->
			def path_string = value.toString()
			file(path_string.startsWith('/') ? path_string : "${base_dir}/${path_string}")
		}

		resolved_input_rows_ch = Channel
			.of(params.summarystats_input)
			// Read and validates against the schema, emits each parsed row as a list of values.
			.flatMap { samplesheetToList(it, "${projectDir}/assets/summarystats_input_schema.json") }
			// Convert each positional samplesheet row into a named map.
			.map { values ->
				def row = [:]
				schema_columns.eachWithIndex { column, index ->
					row[column] = values[index]
				}
				row
			}
			.collect()
			// If liftover is disabled, check that all rows use the same grch and that any supplied grch_bfile matches its row’s grch.
			// The .flatMap then emits the rows individually again (after being gathered into one list by .collect() ).
			.flatMap { rows ->
				if (!params.run_liftover) {
					def grch_builds = rows
						.collect { it.grch.toString().trim().toDouble() }
						.unique()
					
					if (grch_builds.size() > 1) {
						throw new IllegalArgumentException(
							"Multiple GRCh builds detected in grch column but run_liftover is false: ${grch_builds}"
						)
					}

                if (rows.any { row ->
                    def build = row.grch_bfile
                    build != null &&
                        build.toString().trim() &&
                        build.toString() != '[]' &&
                        build.toString().trim().toDouble() != row.grch.toString().trim().toDouble()
                	}) {
						throw new IllegalArgumentException(
							"grch_bfile and grch columns do not match while run_liftover is false."
						)
					}
				}
				// If liftover is enabled, previous checks are skipped
				rows
			}
			// Resolve GWAS path and the PLINK prefix
			.map { row ->
				// These Path objects are reused by the validation and downstream channels.
				def gwas_file = resolve_input_path(row.input)
				def bfile_prefix = resolve_input_path(row.bfile)
				def plink1_files = ['.bed', '.bim', '.fam'].collect { ext -> file("${bfile_prefix}${ext}") }
				def plink2_files = ['.pgen', '.pvar', '.psam'].collect { ext -> file("${bfile_prefix}${ext}") }
				def has_plink1 = plink1_files.every { it.exists() }
				def has_plink2 = plink2_files.every { it.exists() }
				if (!has_plink1 && !has_plink2) {
					throw new IllegalArgumentException("Neither PLINK1 nor PLINK2 files found for prefix: ${bfile_prefix}")
				}
				// Prefer PLINK2 when both complete file sets are available.
				def bfile_files = has_plink2 ? plink2_files : plink1_files
				def bfile_file_type = has_plink2 ? 'plink2' : 'plink1'

				// sdY may be a number/NA or a path; keep a placeholder path so the validation process always receives a path input.
				def sdY_candidate = resolve_input_path(row.sdY)
				def sdY_is_file = sdY_candidate.exists()
				def sdY_file = sdY_is_file 
					? sdY_candidate 
					: file("${projectDir}/assets/NO_SDY_FILE")
				
				// Each emitted item contains the original row map, plus resolved paths
				[
					row: row,
					gwas_file: gwas_file,
//					bfile_prefix: bfile_prefix,
				 	bfile_files: bfile_files,
					bfile_file_type: bfile_file_type,
				 	sdY_file: sdY_file,
					sdY_is_file: sdY_is_file
				]
			}

		// Collect all rows and their files for one single validation task. The metadata is
		// written to a task-local manifest by the process, avoiding a large CLI argument.
		validation_input_ch = resolved_input_rows_ch
			.collect()
			.map { items ->
				def manifest_records = items.collect { validation_entry ->
					def metadata = new LinkedHashMap(validation_entry.row)
					metadata.remove('input')
					[metadata: metadata, sdY_is_file: validation_entry.sdY_is_file]
				}
				tuple(
					manifest_records,
					items.collect { it.gwas_file },
					items.collect { it.sdY_file }
				)
			}

		// Validate input file
		INPUT_COLUMNS_VALIDATION(validation_input_ch)

		validated_input_table = INPUT_COLUMNS_VALIDATION.out.validated_rows.collect().map { validated_rows -> sumstats_input_file }

		downstream_input_ch = resolved_input_rows_ch
			.combine(validated_input_table)
			.map { item, input_table -> item }


		// Build processing and downstream channels from the same resolved file objects.
		bfile_datasets = downstream_input_ch
			.map { item ->
				def row = item.row
				tuple(
					row.process_bfile,
					row.bfile,
					"${row.grch_bfile ?: row.grch}",
					"${params.run_liftover ? 'T' : 'F'}",
					item.bfile_file_type, item.bfile_files
				)
			}
			.unique()
			.branch { process_bfile_flag, bfile_id, grch_bfile, run_liftover, file_type, bfile_dataset ->
				need_processing: process_bfile_flag in [true, "T", "t", "TRUE", "true", "True"] || (grch_bfile.toString() == "37" && run_liftover == "T") || process_bfile_flag == null
				processed: true
			}


		PROCESS_BFILE(bfile_datasets.need_processing, chain_file)

		processed_bfile_datasets = bfile_datasets.processed
			.map { process_bfile_flag, bfile_id, grch_bfile, run_liftover, file_type, bfile_dataset -> 
				tuple(bfile_id, file_type, bfile_dataset)
			}
			.mix(PROCESS_BFILE.out.processed_dataset)

		// Generate a channel with finemapping configuration
		finemapping_config = downstream_input_ch
			.map { item ->
				def row = item.row
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
			.map { bfile_id, study_id, finemap_config, file_type, bfile_dataset ->
				tuple(study_id, finemap_config, file_type, bfile_dataset)
			}

		// Define input channel for munging of GWAS sum stats
		sumstas_input_ch = downstream_input_ch
			.map { item ->
				def row = item.row
				def sdY_string = item.sdY_is_file ? item.sdY_file.name : row.sdY
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
				"accepted_chr": params.chromosomes,
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
				item.gwas_file,
				item.sdY_file
			)
		}
		.combine(processed_bfile_datasets, by: 0)
		.map { bfile_id, study_id, munging_config, gwas_file, sdY_file, file_type, bfile_dataset ->
			tuple(study_id, munging_config, gwas_file, sdY_file, file_type, bfile_dataset)
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
	full_credible_sets = credible_sets_from_finemapping.collect()
	previous_h5ad_studies = Channel.value(file('NO_PREVIOUS_H5AD_STUDIES'))
	if (params.coloc_h5ad_input) {
		credible_sets_from_input = Channel.fromPath(params.coloc_h5ad_input, checkIfExists:true)
		if (params.coloc_skip_previous_studies) {
			GET_STUDY_FROM_ANNDATA(credible_sets_from_input)
			previous_h5ad_studies = GET_STUDY_FROM_ANNDATA.out.previous_study_ids
		}	
		full_credible_sets = full_credible_sets
			.mix(credible_sets_from_input)
			.collect()
	}

	exclude_studies_file = params.coloc_exclude_studies_table ? 
        Channel.fromPath(params.coloc_exclude_studies_table, checkIfExists: true) : 
        Channel.value(file('NO_EXCLUDE_STUDIES'))

	if (params.run_colocalization || params.coloc_h5ad_input) {		
		RUN_COLOCALIZATION( full_credible_sets, previous_h5ad_studies, exclude_studies_file )
	}


	// --- MODULES and JOINT FINEMAPPING ---
	if (params.run_joint_finemapping) {		
		RUN_JOINT_FINEMAPPING( RUN_COLOCALIZATION.out.coloc_results_all, RUN_COLOCALIZATION.out.merged_h5ad)
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
			Channel.fromPath(params.coloc_h5ad_input)
			.map { "$it" }
			.collectFile(
				newLine: true,
				name: "input_h5ad_inputs.txt",
				storeDir: "${params.outdir}/pipeline_inputs"
			)
		}

	workflow.onComplete {
		// At the end store log status
		completionSummary()
	}
}
