# Flanders : Finemapping coLocalization AND plEiotRopy Solver

Flanders is a modular pipeline and toolkit for scalable **fine-mapping** and **colocalization** of genetic association signals across large scale datasets and multiple traits.
Implemented using **Nextflow** and mostly **R**, it separates computationally intensive fine-mapping from downstream colocalization to optimize reusability and performance.

---

## ✅ Requirements
Before running the pipeline, ensure you have the following installed:

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) (v24.04+)
- For environment management, either:
  - [Docker](https://www.docker.com/)
  - [Conda](https://docs.conda.io/en/latest/)
</br>

## 📦 Installation
#### 1. Clone the Repository:
```
git clone https://github.com/Biostatistics-Unit-HT/Flanders.git
```
#### 2. If using Conda, Set up the environment
```bash
conda env create -f pipeline_environment.yml
conda activate flanders_env
```
</br>

## ▶️ Running the pipeline
### Example: Fine-Mapping + Colocalization
```bash
nextflow run Flanders/main.nf    -profile [docker|singularity|conda]    --summarystats_input /path/to/input_table.tsv    --run_colocalization true    -w ./work    -resume
```

### Example: Run Only Colocalization (with existing `.h5ad`)
```bash
nextflow run Flanders/main.nf    -profile [docker|singularity|conda]    --coloc_input /path/to/finemapping_output.h5ad    --run_colocalization true    --coloc_id my_coloc_run    -w ./work    -resume
```

### Quick run with example dataset
```bash
nextflow run Flanders/main.nf -profile test,conda -w ./work
```
</br>

## Pipeline overview
Flanders separates the fine-mapping and colocalization process into two distinct steps:
</br>
</br>
### Step 1: Fine-mapping

#### Inputs
|Input | File description |
|------|------------------|
| [GWAS summary statistics](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/Inputs#gwas-summary-statistics) | Tab/comma-separated files (.tsv/.csv, optionally gzipped) |
| LD reference panel | PLINK-format genotype reference panel (.bed/.bim/.fam) — preferably from the same sample population used in the GWAS |
| [Metadata and GWAS-specific parameters table](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/Inputs#metadata-and-gwas-specific-parameters-table) | A .tsv file listing GWAS summary statistics paths and optional trait-specific parameters for munging, liftover, and fine-mapping |
</br>

This step includes:

#### 1. Munging of GWAS summary statistics
  - Format harmonization (e.g. harmonized column names)
  - Imputation of missing information (e.g. missing allele frequency is calculated from the LD reference panel)
  - Optional liftover to GRCh38
  - Alphabetical ordering of alleles, ensuring the first one in alphabetical order is the effect allele (effect sizes and allele frequencies are flipped/inverted where needed)
  - Conversion of SNP IDs to Flanders internal coding of `"chr"CHR:POS:EA:NEA` (where EA is the first allele in alphabetical order)
    </br>⚠️ This differs from common REF/ALT conventions. This SNP ID format allows for robust variants matching between multiple GWAS summary statistics and LD reference panel. It further allows to reconstruct the effect allele directly from the SNP ID.
</br>

#### 2. Identification of significantly associated genomic regions
>Identifies genomic regions containing significant associated SNPs by employing `Locusbreaker`, an in-house developed algorithm which defines each association peak based on the   distance between the end of a peak and the start of the next one.
</br>`Locusbreaker` first selects all SNPs below a given a p-value threshold (default: 1x10<sup>-6</sup>), identifying groups of SNPs positionally close to each other. If two consecutive SNPs are closer to each other than a set distance threshold (default: 250kb), they are grouped into the same locus, while if they are further apart than the distance threshold, they are used to define the boundaries between peaks. Loci with at least a significant SNPs (default: 5x10<sup>-8</sup>) are retained and their boundaries are enlarged by 100kb to fully capture the shape of the association peak.
</br>


#### 3. Fine-Mapping with SuSiE-RSS
>For each genomic region, finemapping is performed using [SuSiE-RSS](https://stephenslab.github.io/susieR/reference/susie_rss.html) and LD calculated from input PLINK files
</br>⚠️ Whenever possible, it is best to use in sample LD, especially for molecular omic phenotypes where the explained variance can be very large.
</br>

#### 4. Saving of fine-mapping results to AnnData object
>Log approximate Bayes factors (lABFs) and metadata for the 99% credible sets are stored in an [AnnData object](https://anndata.dynverse.org/index.html) (.h5ad).
</br>
</br>

### Step 2: Colocalization analysis

#### Inputs
|Input | File description |
|------|------------------|
| Fine-mapping AnnData    |	An .h5ad file containing lABFs and metadata of credible sets (output from the fine-mapping step)  |


#### 1. Generation of colocalization guide table
>Pairs of credible sets sharing at least one SNP are identified and listed in the analysis guide table (it is not possible for credible sets to colocalize without sharing at least a SNP).

#### 2. iCOLOC Test
>Performs pair-wise colocalization for pair of credible sets listed in the guide table by employing `iCOLOC`, a framework extending traditional [colocalization analysis using Bayes Factors](https://chr1swallace.github.io/coloc/reference/coloc.abf.html) by imputing lABFs of SNPs outside of credible sets to the minimum lABF value in the locus.
This approach:
1. Significantly reduces storage requirements by saving in the AnnData object only exact lABF values of credible sets SNPs
2. Enhances colocalization accuracy compared to tradional coloc by reducing false positives due to two causal SNPs being in strong LD.
</br>
</br>


## Output
| Output Type                                     | Description                                                                                                      |
| ----------------------------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| `gwas_and_loci_tables/*_dataset_aligned.tsv.gz` | Harmonized (and optionally lifted) GWAS summary statistics                                                       |
| `gwas_and_loci_tables/*_loci.tsv`               | Boundaries of identified association regions and GWAS summary statistics for the sentinel SNP                    |
| `finemapping_exceptions/`                       | Multiple tables reporting information about loci that were not fine-mapped with the standard procedure or at all |
| `finemapping/*_susie_finemap.rds` _(optional)_  | Individual RDS files for each fine-mapped locus                                                                  |
| `anndata/*.h5ad`                                | AnnData object with lABFs, CS metadata and SNP annotations resulting from fine-mapping                           |
| `coloc/coloc_guide_table.csv`                   | Colocalization analysis guide table, listing all colocalization tests performed                                  |
| `coloc/*_colocalization.table.*.tsv`            | Colocalization analysis results (all, filterd by PPH4 threshold and filtered by PPH3 threshold)                  |
</br>
</br>


## Credits
Developed by the Biostatistics and Genome Analysis Units at [Human Technopole](https://humantechnopole.it/en/)<br>
[Arianna Landini](mailto:arianna.landini@fht.org)<br>
[Sodbo Sharapov](mailto:sodbo.sharapov@fht.org)<br>
[Edoardo Giacopuzzi](mailto:edoardo.giacopuzzi@fht.org)<br>
[Bruno Ariano](mailto:bruno.ariano@fht.org)<br>
[Nicola Pirastu](mailto:nicola.pirastu@fht.org)<br>
