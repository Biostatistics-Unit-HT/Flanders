# Flanders : Finemapping coLocalization AND plEiotRopy Solver	

Set of tools and pipeline to efficiently colocalise association signals across large sets of traits.
Colocalization is mainly composed of two steps, which are generally merged in a single function finemapping of the locus and colocalization itself. Given that the expensive part of the computation is the finemapping, to minimise its cost we have split the process in two parts:    

## Step 1: Munging, significant genomic region identification and Finemapping with Susie

1) Munging and harmonising of summary statistics
   GWAS summary statistics are lifted over to build 38 and snps IDs are converted to the FIGARO internal coding which looks like chromosome:position_bp:Effect_allele:Non_effect_allele. Liftover is strongly suggested although optional.
   Alleles are ordered according to alphabetical order so that the effect allele is always the one that comes first in alphabetical order. Effect sizes are flipped accordingly if needed.
   WARNING: This is different from the often used REF/ALT schema. This choice allows to be able to reconstruct the effect allele directly from the SNP id without any further information
   
2) Identification of associated genomic regions.
   The next step in the workflow is to identify the genomic regions (ie. chr1:1247282-1256282) which contain significant associated snps.
   To identify genomic regions, we employ an in-house developed algorithm called Locusbreaker, which identifies each association peak based on the distance between the end of a peak and the start of the next one.
   Briefly, given a defined p-value threshold higher than the genome-wide significant threshold (default: 1x10-6), it first selects all the snps below such threshold. This will create groups of snps close to each other if they belong to the 
   same association peak. If two consecutive SNPs are closer to each other than a set threshold (default: 250kb), then they are kept together in the same genomi region, while if two snps are further than the threshold, they define the 
   boundaries between peaks.  Genomc regions with at least a significant SNPs are retained and their boundries are enlarged by 100kb to include the surrounding data which is usefull for finemapping.

3) For each genomic region finemapping is performed using Susie-rss. This step also requires reference genotypes to compute LD between the SNPs. WARNING: Whenever possible, it is best to use in sample LD. This is particularly true when 
   finemapping molecular omic phenotypes where the explained variance can be very large.

4) 99% credible sets from the fine mapping step are QC'd based on [SODBO to describe], and lABFs for each SNP in each credible set are stored in an AnnData object (https://anndata.dynverse.org/index.html). In addition to the lABFs for each snp 
   we store the average lABF of all the SNPs which are not in the 99% cs to be used later in the iColoc step.
   AnnData allows to efficiently store the lABF from all credible sets in a single sparse matrix while allowing to have metadata for credible sets and snps in the obs and var matrices. This way of storing the data has several advantages: a) all 
   the information needed for colocalization are stored in a single object b) it allows to store and reuse the finemapping results c) the resulting file is many times smaller than the initial files, for example ........... 

## Step 2: Colocalization analysis

1) To maximise efficiency and reduce as much as possible, the first step is to compute all pairs of credible sets which share at least a SNP. Given the way colocalization is computed it is not possible for two credible sets to colocalize if they don't share at least a SNP.
2) Colocalization for each identified pair is performed in parallel and results are all stored in a single final table using iCOLOC. iCOLOC imputes all the values of lABF outside the credible set to their average. We have run extensive simulations under many different scenarios and this performs as well as standard coloc for discovery while it reduces false positives due to two causal SNPs being in LD.








# Running Flanders


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
