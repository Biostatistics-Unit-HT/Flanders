# Flanders : Finemapping coLocalization AND plEiotRopy Solver BETA

Set of tools and pipeline to efficiently colocalise association signals across large sets of traits.
Colocalization is mainly composed of two steps, which are generally merged in a single function finemapping of the locus and colocalization itself. Given that the expensive part of the computation is the finemapping, to minimise its cost we have split the process in two parts:    
</br>
## Step 1: Munging, significant genomic region identification and Finemapping with Susie

### 1) Munging and harmonising of summary statistics
   GWAS summary statistics are lifted over to build 38 and snps IDs are converted to the Flanders internal format: `chromosome:position_bp:Effect_allele:Non_effect_allele`. Alleles are ordered according to alphabetical order so that the effect allele is always the one that comes first in alphabetical order. For example rs1570295954 would be called chr1:146054827:A:G and all effect sizes would be flipped to match the A allele being the coded one.
   
   __WARNING: This is different from the commonly used REF/ALT schema.__
   
   This choice allows reconstruction of the effect and non-effect alleles directly from the SNP id. 
   Liftover is strongly suggested, although optional. </br>
      
### 2) Identification of associated genomic regions.
   The next step in the workflow is to identify the genomic regions (ie. chr1:1247282-1256282) which contain significant associated snps. To identify genomic regions, we employ an in-house developed algorithm called Locusbreaker, which identifies each association peak based on the distance between the end of a peak and the start of the next one.
   Briefly, given a defined p-value threshold higher than the genome-wide significant threshold (default: 1x10<sup>-6</sup>), it first selects all the SNPs below such threshold. This will create groups of SNPs close to each other if they belong to the 
   same association peak. If two consecutive SNPs are closer to each other than a set threshold (default: 250kb), then they are kept together in the same genomic region, while if two SNPs are further than the threshold, they define the 
   boundaries between peaks.  Genomic regions with at least a significant SNP are retained and their boundaries are enlarged by 100kb to include the surrounding data, which is useful for finemapping.

### 3) Finemapping of each locus
For each genomic region finemapping is performed using Susie-rss. This step also requires reference genotypes to compute LD between the SNPs.</br>

__WARNING: In-sample LD should always be preferred when possible.__ 

Out of sample LD is particularly problematic and will lead to a very large number of false positives cs when fine-mapping molecular omic phenotypes where the explained variance can be very large.

### 4) Credible sets are stored in csAnnData 
99% credible sets from the fine mapping step are QC'd based , and lABFs for each SNP in each credible set are stored in an [csAnnData](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/csAnnData-specifications) object. This is a specific implementation of general AnnData described [here](https://anndata.dynverse.org/index.html). In addition to the lABFs for each snp in the 99% credible set we store the average lABF of all the SNPs which are not in the 99% cs to be used later in the iColoc step.
Furthermore, for each SNP we store in the additional layers the beta and se extracted directly from the IBSS joint multivariable regression from Susie. 

csAnnData allows to efficiently store the lABF from all credible sets in a single sparse matrix while allowing metadata both for credible sets and snps in the obs and var matrices.

This format has several advantages: 
- efficient storing the information usefull for colocalization analyses in a very limited space (~1000x less storage compared to the original sumstats). 
- easy merging with other csAnnData objects containing finemapping results from other studies for cross study colocalization
- R and pythion compatibility
</br>

## Step 2: Colocalization analysis

1) To maximise efficiency and reduce as much as possible, the first step is to compute all pairs of credible sets which share at least a SNP. Given the way colocalization is computed it is not possible for two credible sets to colocalize if they don't share at least a SNP.
2) Colocalization for each identified pair is performed in parallel and results are all stored in a single final table using iCOLOC. iCOLOC imputes all the values of lABF outside the credible set to their average. We have run extensive simulations under many different scenarios and this performs as well as standard coloc for discovery while it reduces false positives due to two causal SNPs being in LD.
</br>
</br>

## Requirements
Before running the pipeline, ensure you have the following installed:

Nextflow (>= version 24.04) </br>
Docker / Singularity / Conda
</br>

## Installation
Clone the repository and navigate into the project directory:
```
git clone https://github.com/Biostatistics-Unit-HT/Flanders.git
```
</br>

## Input Data
### Fine-mapping
- [GWAS summary statistics (.csv/.csv.gz/.tsv/.tsv.gz)](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/Inputs#gwas-summary-statistics).
- [Genotypes for LD reference panel in Plink .bed/.bim/.fam format](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/Inputs#genotypes-for-ld-reference-panel). __Providing, if possible, in sample LD greatly improves the accuracy of fine-mapping__
- [Metadata and GWAS-specific parameters table](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/Inputs#metadata-and-gwas-specific-parameters-table).

### Colocalisation
- [csAnnData](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/csAnnData-specifications)
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
- Fine-mapping output in a single [csAnnData](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/csAnnData-specifications) file (lABFs, conditional beta and se, qc metrics and metadata for each cs)
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
