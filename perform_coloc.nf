#!/usr/bin/env nextflow

nextflow.enable.dsl=2 // specify the Domain Specific Language version to be used in the Nextflow script

// Source all processes
include {FIND_CS_OVERLAP_BY_CHR} from "$projectDir/modules/find_cs_overlap_by_chr"
include {COLOC} from "$projectDir/modules/coloc"
include {IDENTFY_REG_MODULES} from "$projectDir/modules/identify_reg_modules"

// Need to specify here the chunck size! Rather than hardcode it the channel

// Define the main workflow
workflow {

// Define a channel for each process

// Define input channel for identifying overlapping cs among submitted conditional datasets
  coloc_info_table = Channel.fromPath(params.inputFileList)

// Define input channel of unique chromosomes present in the list of finemapped conditional datasets
  all_cond_datasets_cs = Channel
    .of(file(params.inputFileList))
    .splitCsv(header:true, sep:"\t")
    .map{ row ->
      [
        "chr_cs": row.chr
      ]
    }
    .unique()
    .combine(coloc_info_table)

  // Run FIND_CS_OVERLAP process on all_cond_datasets_cs channel
  FIND_CS_OVERLAP_BY_CHR(all_cond_datasets_cs)


  // Define input channel for performing pair-wise colocalisation analysis (in batches)
  coloc_pairs_by_batches = FIND_CS_OVERLAP_BY_CHR.out.coloc_pairwise_guide_table
    .splitText(by:5000, keepHeader:true, file:true) //// can the batch size not be hardcoded??

  // Run COLOC process on coloc_pairs_by_batches channel
  COLOC(coloc_pairs_by_batches)

  // Define input channel for identifying regulatory modules - collect all tables (coloc performed in batches of n pairwise tests)
  coloc_results_all = COLOC.out.colocalization_table_all_by_chunk
    .collectFile(name: 'colocalization_table_all_merged.txt', keepHeader: true)
    .combine(coloc_info_table)
  
  // Run COLOC process on coloc_pairs_by_batches channel
  IDENTFY_REG_MODULES(coloc_results_all)

}