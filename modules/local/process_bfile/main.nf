process PROCESS_BFILE {
    tag "${bfile_dataset[0].baseName}"
    label "process_small"

    input:
    tuple val(bfile_id), val(genome_build), val(run_liftover), path(bfile_dataset)
    path chain_file

    output:
    tuple val(bfile_id), path("${bfile_dataset[0].baseName}.GRCh38.alpha_sorted_alleles.{bed,bim,fam}"), emit: processed_dataset
    
    script:
	"""
    s01_alpha_sort_alleles_snpid.R \
        --bfile ${bfile_dataset[0].baseName} \
        --grch ${genome_build} \
        --run_liftover ${run_liftover}
	"""
}