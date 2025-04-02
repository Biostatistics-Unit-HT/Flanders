#!/usr/bin/env nextflow

nextflow.enable.dsl=2 // specify the Domain Specific Language version to be used in the Nextflow script

// Set parameter for pipeline version - only for documenting sake
//params.pipe_vers="1.0"
// Source all processes
include { MUNG_AND_LOCUS_BREAKER  }  from "./modules/local/mung_and_locus_breaker"
include { SUSIE_FINEMAPPING       }  from "./modules/local/susie_finemapping"
include { COJO_AND_FINEMAPPING    }  from "./modules/local/cojo_and_finemapping"
include { APPEND_TO_MASTER_COLOC  }  from "./modules/local/append_to_master_coloc"
include { APPEND_TO_IND_SNPS_TAB  }  from "./modules/local/append_to_ind_snps_tab"
include { INPUT_COLUMNS_VALIDATION } from "./modules/local/input_columns_validation"
include { samplesheetToList } from 'plugin/nf-schema'


workflow {

  // Define asset files
  chain_file = file("${projectDir}/assets/hg19ToHg38.over.chain")

  // Ensure the folder to store not finemapped loci exists
  file("${params.outdir}/results/not_finemapped_loci").mkdirs()

  // Use nf-schema to read and validate the sample sheet
  samplesheetToList(params.inputFileList, params.schema)

  // Validate input file
  // In case we are running a test profile, we need to set the base_dir to the projectDir
  base_dir = params.is_test_profile ? "${projectDir}" : "${launchDir}"
  INPUT_COLUMNS_VALIDATION(file(params.inputFileList), base_dir)

  // Define input channel for munging of GWAS sum stats
  gwas_input = INPUT_COLUMNS_VALIDATION.out.table_out
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
  // Run MUNG_AND_LOCUS_BREAKER process on gwas_input channel
  MUNG_AND_LOCUS_BREAKER(gwas_input, chain_file)

// Output channel of LOCUS_BREAKER *** process one locus at a time ***
  MUNG_AND_LOCUS_BREAKER.out.loci_table
    .splitCsv(header:true, sep:"\t")
    .map{ row -> tuple(
      [
        "study_id": row.study_id
      ],
      [
        "chr":row.chr,
        "start":row.start,
        "end":row.end,
        "locus_size":row.locus_size,
        "phenotype_id":row.phenotype_id
      ]
      )
    }
    .branch { study_meta, locus_meta ->
      small: locus_meta.locus_size.toInteger() < params.large_locus_size
      large: locus_meta.locus_size.toInteger() >= params.large_locus_size
    }
    .set { loci_for_finemapping }

// Large loci are collected in a file and published
  loci_for_finemapping.large
    .map{ study, locus ->
      [study.study_id, locus.chr, locus.start, locus.end, locus.locus_size, locus.phenotype_id].join('\t')
    }
    .collectFile(
      newLine: true, 
      name: "large_loci.tsv", storeDir: "${params.outdir}/results/not_finemapped_loci",
      seed: "study_id\tchr\tstart\tend\tlocus_size\tphenotype_id\n")

// Define metadata channel for COJO/SUSIE FINEMAPPING
// Only small loci will be used for finemapping
  finemapping_input = Channel  
    .of(file(params.inputFileList))
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
        "skip_dentist": row.skip_dentist,
        "maf": row.maf,
        "hole": row.hole,
        "cs_thresh": row.cs_thresh
      ]
    )
  }
  .combine(loci_for_finemapping.small, by:0)
  .combine(MUNG_AND_LOCUS_BREAKER.out.dataset_munged_aligned, by:0)
//  .map{study_id, meta_finemapping, meta_loci, gwas_final, gwas_final_index -> tuple(study_id, meta_finemapping+meta_loci, gwas_final, gwas_final_index)}


// Run SUSIE_FINEMAPPING process on finemapping_input channel
  SUSIE_FINEMAPPING(finemapping_input, params.outdir)


// Run COJO on failed SUSIE loci (only for specific errors!! Stored in the R script of susie)
  COJO_AND_FINEMAPPING(SUSIE_FINEMAPPING.out.failed_susie_loci, params.outdir)
   
// Append all to independent SNPs table /// What if the channel is empty?? Can you check and behave accordingly?
  append_ind_snps = COJO_AND_FINEMAPPING.out.ind_snps_table
    .mix(SUSIE_FINEMAPPING.out.ind_snps_table)
    .groupTuple()
    .map{ it.flatten().collect() }
    .map{ tuple( it[0], it[1..-1])}

  APPEND_TO_IND_SNPS_TAB(append_ind_snps)


// Append all to coloc_info_master_table
  append_input_coloc = COJO_AND_FINEMAPPING.out.cojo_info_coloc_table
    .mix(SUSIE_FINEMAPPING.out.susie_info_coloc_table)
    .groupTuple()
    .map{ it.flatten().collect() }
    .map{ tuple( it[0], it[1..-1])}

  APPEND_TO_MASTER_COLOC(append_input_coloc)
}