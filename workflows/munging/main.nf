include { MUNG_AND_LOCUS_BREAKER  }  from "../../modules/local/mung_and_locus_breaker"

workflow RUN_MUNGING {
	take:
    sumstats_input // input channel for munging of GWAS sum stats
		chain_file // file of hg19ToHg38 chain

	main:
  // Ensure the folder to store not finemapped loci exists
  file("${params.outdir}/results/not_finemapped_loci").mkdirs()
  
  // Run MUNG_AND_LOCUS_BREAKER process on gwas_input channel
  MUNG_AND_LOCUS_BREAKER(sumstats_input, chain_file)

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
      seed: "study_id\tchr\tstart\tend\tlocus_size\tphenotype_id")

  emit:
    finemapped_loci = loci_for_finemapping.small
    munged_stats = MUNG_AND_LOCUS_BREAKER.out.dataset_munged_aligned
}