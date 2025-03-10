Example of `nf-hcoloc_run.sbatch` to run the pipeline.
For the first nextflow run (`nextflow run nf-hcoloc/main.nf`), user must provide to `--inputFileList` a table listing all required metadata for each GWAS summary statistics you want to analysis (see example below).

In the second nextflow run (`nextflow run nf-hcoloc/perform_coloc.nf`), user must assemble and provide to `--inputFileList` a file listing all credible sets for the traits to be tested - this is achievable by catting *_coloc_info_master_table.tsv (found in results/coloc_info_tables) for all the traits of interest (an helper script to do so is coming). Other arguments are:
`--pph3` and `--pph4` are the significance threshold for posterior probability of hypothesis 3 (no colocalisation) and 4 (colocalisation) respectively.
`--coloc_id` is a label identifier for your colocalisation analysis (all output files will have this id as prefix).
`--results_path` is the path where the `results` folder, containing pre-processing outputs for ALL the traits you're planning to colocalise, lives.



```
#!/bin/bash
#SBATCH --job-name nf-hcoloc
#SBATCH --output nf-hcoloc_%A.log
#SBATCH --partition cpuq
#SBATCH --cpus-per-task 1
#SBATCH --mem 8G
#SBATCH --time 20-00:00:00

module load nextflow/23.10.0 #singularity/3.8.5
 
nextflow run nf-hcoloc/main.nf \
  --inputFileList hcoloc_submission_table_example.tsv -profile ht_cluster -resume

nextflow run nf-hcoloc/perform_coloc.nf \
  --inputFileList all_cond_datasets_cs_example.tsv \
  --pph3 0.75 \
  --pph4 0.75 \
  --coloc_id test \
  --results_path /path/to/results/folder/ \
  -profile ht_cluster -resume
```

<br>
<br>
<br>
<br>

Example of input table `hcoloc_submission_table_example.tsv`


| input | study_id | chr_lab | pos_lab | rsid_lab | a1_lab | a0_lab | freq_lab | n_lab | effect_lab | se_lab | pvalue_lab | type | sdY | s | grch | p_thresh1 | p_thresh2 | hole | bfile | p_thresh3 | p_thresh4 | maf | is_molQTL | key | cs_thresh | skip_dentist |
|-------|----------|---------|---------|----------|--------|--------|----------|-------|------------|--------|------------|------|-----|---|-----|-----------|-----------|------|-------|-----------|-----------|-----|-----------|-----|-----------|--------------|
| /group/pirastu/prj_014_huvec_coloc/input/diseases/GWAS_Cardioembolic_Stroke_Eur_Mishra_2022_Nature_hg38.tsv.gz | Cardioembolic_Stroke_Eur_Mishra_2022_Nature | CHROM | GENPOS | SNP | ALLELE1 | ALLELE0 | A1FREQ | N | BETA | SE | P | cc | NA | 0.1 | 38 | 5.00E-08 | 1.00E-05 | 250000 | /ssu/bsssu/ghrc38_reference/ukbb_all_chrs_grch38_maf0.01_30000_random_unrelated_white_british_alpha_sort_alleles | 1.00E-04 | 1.00E-06 | 1.00E-04 | FALSE | NA | 0.99 |    FALSE     |
| /group/pirastu/prj_014_huvec_coloc/input/ATAC/ATAC_cis_eQTLs_chr22.tsv.gz | ATAC_chr22 | CHROM | GENPOS | SNP | ALLELE1 | ALLELE0 | A1FREQ | N | BETA | SE | P | quant | NA | NA | 38 | 5.00E-08 | 1.00E-05 | 250000 | /ssu/bsssu/ghrc38_reference/ukbb_all_chrs_grch38_maf0.01_30000_random_unrelated_white_british_alpha_sort_alleles | 1.00E-04 | 1.00E-06 | 1.00E-04 | TRUE | trait | 0.99 |    FALSE     |
