# Flanders : Finemapping coLocalization AND plEiotRopy Solver BETA

Set of tools and pipeline to efficiently colocalise association signals across large sets of traits.
Colocalization is mainly composed of two steps, which are generally merged in a single function finemapping of the locus and colocalization itself. Given that the expensive part of the computation is the finemapping, to minimise its cost we have split the process in two parts:    
</br>
## Step 1: Munging, significant genomic region identification and Finemapping with Susie

1) Munging and harmonising of summary statistics
   GWAS summary statistics are lifted over to build 38 and snps IDs are converted to the Flanders internal coding which looks like `chromosome:position_bp:Effect_allele:Non_effect_allele`. Liftover is strongly suggested although optional. Alleles are ordered according to alphabetical order so that the effect allele is always the one that comes first in alphabetical order. Effect sizes are flipped accordingly if needed.</br>
   __WARNING: This is different from the often used REF/ALT schema.__ This choice allows to be able to reconstruct the effect allele directly from the SNP id without any further information
   
3) Identification of associated genomic regions.
   The next step in the workflow is to identify the genomic regions (ie. chr1:1247282-1256282) which contain significant associated snps. To identify genomic regions, we employ an in-house developed algorithm called Locusbreaker, which identifies each association peak based on the distance between the end of a peak and the start of the next one.
   Briefly, given a defined p-value threshold higher than the genome-wide significant threshold (default: 1x10<sup>-6</sup>), it first selects all the snps below such threshold. This will create groups of snps close to each other if they belong to the 
   same association peak. If two consecutive SNPs are closer to each other than a set threshold (default: 250kb), then they are kept together in the same genomi region, while if two snps are further than the threshold, they define the 
   boundaries between peaks.  Genomc regions with at least a significant SNPs are retained and their boundries are enlarged by 100kb to include the surrounding data which is usefull for finemapping.

4) For each genomic region finemapping is performed using Susie-rss. This step also requires reference genotypes to compute LD between the SNPs.</br>__WARNING: Whenever possible, it is best to use in sample LD.__ This is particularly true when finemapping molecular omic phenotypes where the explained variance can be very large.

5) 99% credible sets from the fine mapping step are QC'd based on [SODBO to describe], and lABFs for each SNP in each credible set are stored in an AnnData object (https://anndata.dynverse.org/index.html). In addition to the lABFs for each snp 
   we store the average lABF of all the SNPs which are not in the 99% cs to be used later in the iColoc step.
   AnnData allows to efficiently store the lABF from all credible sets in a single sparse matrix while allowing to have metadata for credible sets and snps in the obs and var matrices. This way of storing the data has several advantages: a) all 
   the information needed for colocalization are stored in a single object b) it allows to store and reuse the finemapping results c) the resulting file is many times smaller than the initial files, for example ........... 
</br>

## Step 2: Colocalization analysis

1) To maximise efficiency and reduce as much as possible, the first step is to compute all pairs of credible sets which share at least a SNP. Given the way colocalization is computed it is not possible for two credible sets to colocalize if they don't share at least a SNP.
2) Colocalization for each identified pair is performed in parallel and results are all stored in a single final table using iCOLOC. iCOLOC imputes all the values of lABF outside the credible set to their average. We have run extensive simulations under many different scenarios and this performs as well as standard coloc for discovery while it reduces false positives due to two causal SNPs being in LD.
</br>
</br>

## Requirements
Before running the pipeline, ensure you have the following installed:

[Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) (version 24.04 or higher)</br>
[Docker](https://www.docker.com/) or [Conda](https://docs.conda.io/en/latest/) for environment management
</br>
</br>

## Installation
1. Clone the Repository:
```
git clone https://github.com/Biostatistics-Unit-HT/Flanders.git
```
2. Set Up the Environment:

   Using Conda:
```
conda env create -f pipeline_environment.yml
conda activate flanders_env
```
</br>
</br>

## Input Data
- GWAS summary statistics in tab or comma separated format, even gzipped
- LD panel in plink .bed/.bim/.fam format
- Input table in .tsv format

  

#### A parenthesis - the importance of in sample LD reference
Providing, if possible, in sample LD greatly improves the accuracy of fine-mapping
</br>

## Basic Usage
Run the pipeline using:

```
nextflow run Flanders/main.nf \
   -profile [docker|singularity|conda] \
   --summarystats_input /path/to/input_table.tsv \
   --run_liftover T \
   --run_colocalization T \
   -w ./work \
   -resume
```
| Parameter                     | Description                                                     |
|-------------------------------|-----------------------------------------------------------------|
| `--summarystats_input`        | Path to input table                                             |
| `--run_liftover`              | Whether to lift input data from hg37 to hg38                    |
| `--run_colocalization`        | Wheter to follow-up fine-mapping with colocalization            |




## Configuration
Check nextflow.config for all customizabile parameters

| Parameter                     | Description                                                               |
|-------------------------------|---------------------------------------------------------------------------|
| **Input data** |                                       |
| `--summarystats_input`        | Input table with summary statistics and parameters for mungin and finemap |
| `--coloc_input `              | Input file for coloc                                                      |
| `--coloc_id`                  | ID label for the coloc analysis output                                    |
| **Output settings** |                                       |
| `--outdir`                    | Directory where output will be generated                        |
| **Munging and finemapping settings** |                                       |


</br>

## Quick run with example dataset
```
nextflow run Flanders/main.nf -profile test,conda -w ./work
```
</br>

## Output
- Processed (and lifted) LD reference plink bfiles (.bed/.bim/.fam)
- Munged and index (and lifted) GWAS summary statistics
- Significant loci table
- Fine-mapping output per locus in multiple .rds files (lABFs, conditional beta and se, qc metrics and metadata for each cs)
- Fine-mapping output in a single AnnData file (lABFs, conditional beta and se, qc metrics and metadata for each cs). Check [here](https://github.com/Biostatistics-Unit-HT/flanders_r?tab=readme-ov-file#anndata-column-specifications) a detailed description of AnnData structure.
- Tables collecting pairwise colocalisation test results
</br>

## Performances
</br>

## Credits
Developed by Biostatistics and Genome Analysis Units at [Human Technopole](https://humantechnopole.it/en/)<br>
[Arianna Landini](mailto:arianna.landini@fht.org)<br>
[Sodbo Sharapov](mailto:sodbo.sharapov@fht.org)<br>
[Edoardo Giacopuzzi](mailto:edoardo.giacopuzzi@fht.org)<br>
[Bruno Ariano](mailto:bruno.ariano@fht.org)<br>
[Nicola Pirastu](mailto:nicola.pirastu@fht.org)<br>
