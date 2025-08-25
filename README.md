# Flanders : Finemapping coLocalization AND plEiotRopy Solver

Flanders is a modular pipeline and toolkit for scalable **fine-mapping** and **colocalization** of genetic association signals across large scale datasets and multiple traits.
Implemented using **Nextflow** and mostly **R**, it separates computationally intensive fine-mapping from downstream colocalization to optimize reusability and performance.

---

## ✅ Requirements
Before running the pipeline, ensure you have the following installed:

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) (v24.04+)
- For environment management, one of:
  - [Docker](https://www.docker.com/)
  - [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html)
  - [Conda](https://docs.conda.io/en/latest/)
</br>


## ▶️ Running the pipeline
### Example: Fine-Mapping + Colocalization
```bash
nextflow run Biostatistics-Unit-HT/Flanders -r 1.0    -profile [docker|singularity|conda]    --summarystats_input /path/to/input_table.tsv    --run_colocalization true    --finemap_id my_finemap_run    --coloc_id my_coloc_run    -w ./work    -resume
```

### Example: Run Only Colocalization (with existing `.h5ad`)
```bash
nextflow run Biostatistics-Unit-HT/Flanders -r 1.0    -profile [docker|singularity|conda]    --coloc_input /path/to/finemapping_output.h5ad    --run_colocalization true    --coloc_id my_coloc_run    -w ./work    -resume
```

### Quick run with example dataset
```bash
nextflow run Biostatistics-Unit-HT/Flanders -r 1.0 -profile test,singularity -w ./work
```
</br>

## 🧠 Pipeline overview
Flanders separates the fine-mapping and colocalization process into two distinct steps:
</br>
</br>
### Step 1: Fine-mapping

#### Required Inputs
 Input | Description |
|-------|-------------|
| [GWAS summary statistics](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/Fine‐mapping-inputs#gwas-summary-statistics) | `.tsv`/`.csv` (optionally gzipped) |
| LD reference panel | PLINK-format reference panel (`.bed/.bim/.fam`) — preferably from the same sample population used in the GWAS |
| [Metadata and GWAS-specific parameters table](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/Fine‐mapping-inputs#metadata-and-gwas-specific-parameters-table) | `.tsv` file listing GWAS summary statistics paths and trait-specific parameters |
</br>

#### Steps

1. **Munging of GWAS summary statistics**
  - Format harmonization and imputation of missing information (e.g. missing allele frequency calculated from the LD reference panel)
  - Optional liftover to GRCh38
  - Optional restriction of analysis to enlisted chromosomes
    </br>⚠️ In addition to autosomes, chromosomes **X** and **Y** are also accepted.
  - Alphabetical ordering of alleles, ensuring the first one in alphabetical order is the effect allele (effect sizes and allele frequencies are flipped/inverted where needed)
  - Conversion of SNP IDs to Flanders internal coding of `"chr"CHR:POS:EA:NEA` ***where EA is the first allele in alphabetical order***
    </br>⚠️ This differs from common REF/ALT conventions and allows for robust variants matching between multiple GWAS summary statistics and LD reference panel.
</br>

#### 2. Identification of significantly associated genomic regions
  - Identifies genomic regions containing significant associated SNPs by employing `Locusbreaker`, an in-house developed algorithm which defines each association peak based on the   distance between the end of a peak and the start of the next one.

<details>
  
`Locusbreaker` first selects all SNPs below a given a p-value threshold (suggested value 1x10<sup>-6</sup>, customizable at the column `p_thresh2` of the [Metadata and GWAS-specific parameters table](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/Fine‐mapping-inputs#metadata-and-gwas-specific-parameters-table)), identifying groups of SNPs positionally close to each other.
</br>If two consecutive SNPs are closer to each other than a set distance threshold (suggested value 250kb, customizable at the column `hole` of the [Metadata and GWAS-specific parameters table](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/Fine‐mapping-inputs#metadata-and-gwas-specific-parameters-table)), they are grouped into the same locus, while if they are further apart than the distance threshold, they are used to define the boundaries between peaks.
</br>Loci with at least a significant SNPs (suggested value 5x10<sup>-8</sup>, customizable at the column `p_thresh1` of the [Metadata and GWAS-specific parameters table](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/Fine‐mapping-inputs#metadata-and-gwas-specific-parameters-table)) are retained and their boundaries are enlarged by 100kb to fully capture the shape of the association peak.
</details>
  
</br>


#### 3. Fine-Mapping with SuSiE-RSS
  - For each genomic region, finemapping is performed using [SuSiE-RSS](https://stephenslab.github.io/susieR/reference/susie_rss.html) and LD calculated from input PLINK files
</br>⚠️ Whenever possible, in sample LD is strongly recommended (especially for molecular omic phenotypes where the explained variance can be very large).
</br>⚠️ Be aware that only SNPs in common between the GWAS summary statistics and the LD reference panel are taken into account for fine-mapping, while all other SNPs are discarded.
</br>

#### 4. Saving fine-mapping results to AnnData object
  - Log approximate Bayes factors (lABFs) and metadata for the 99% credible sets are stored in an [AnnData object](https://anndata.dynverse.org/index.html) (.h5ad).
</br>
</br>

### Step 2: Colocalization analysis

#### Inputs
|Input | File description |
|------|------------------|
| [Fine-mapping AnnData](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/csAnnData-specifications)    |	An `.h5ad` file containing lABFs and metadata of credible sets (output from the fine-mapping step)  |
</br>

#### Steps

#### 1. Generation of colocalization guide table
- Lists all pairs of credible sets that share at least one SNP  (it is not possible for credible sets to colocalize without sharing at least a SNP).

#### 2. Colocalization with iCOLOC
  - Performs pair-wise colocalization for pair of credible sets listed in the guide table by employing `iCOLOC`, a framework extending traditional [colocalization analysis using Bayes Factors](https://chr1swallace.github.io/coloc/reference/coloc.abf.html) by imputing lABFs of SNPs outside of credible sets to the minimum lABF value in the locus.

<details>
  
iCOLOC approach allows to:
  1. Significantly reducing storage requirements by saving in the AnnData object only exact lABF values of credible sets SNPs
  2. Enhancing colocalization accuracy compared to tradional coloc by reducing false positives due to two causal SNPs being in strong LD.
</details>

</br>
</br>


## 📁 Output
| Output Type                                     | Description                                                                                                      |
| ----------------------------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| `gwas_and_loci_tables/*_dataset_aligned.tsv.gz` | Harmonized (and optionally lifted) GWAS summary statistics                                                       |
| `gwas_and_loci_tables/*_loci.tsv`               | Boundaries of identified association regions and GWAS summary statistics for the sentinel SNP                    |
| `finemapping_exceptions/`                       | [Multiple tables reporting information about loci that were not fine-mapped with the standard procedure or at all](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/Fine%E2%80%90mapping-exceptions) |
| `finemapping/*_susie_finemap.rds` _(optional)_  | Individual RDS files for each fine-mapped locus                                                                  |
| `anndata/*.h5ad`                                | AnnData object with lABFs, CS metadata and SNP annotations resulting from fine-mapping                           |
| `coloc/coloc_guide_table.csv`                   | Colocalization analysis guide table, listing all colocalization tests performed                                  |
| `coloc/*_colocalization.table.*.tsv`            | Colocalization analysis results (all, filterd by PPH4 threshold and filtered by PPH3 threshold)                  |
</br>
</br>


## 👩‍🔬 Credits
Developed by the Biostatistics and Genome Analysis Units at [Human Technopole](https://humantechnopole.it/en/)<br>
-  [Arianna Landini](mailto:arianna.landini@fht.org)<br>
-  [Sodbo Sharapov](mailto:sodbo.sharapov@fht.org)<br>
-  [Edoardo Giacopuzzi](mailto:edoardo.giacopuzzi@fht.org)<br>
-  [Bruno Ariano](mailto:bruno.ariano@fht.org)<br>
-  [Nicola Pirastu](mailto:nicola.pirastu@fht.org)<br>
