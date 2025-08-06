# Quick start

## Requirements

Before running the pipeline, ensure you have the following installed:

- Nextflow (>= version 24.04)
- Docker / Singularity / Conda

## Basic Usage

Run the pipeline directly from our repository using:

```bash
nextflow run Biostatistics-Unit-HT/Flanders [-r <version/commit>] \
   -profile [docker|singularity|conda] \
   --summarystats_input /path/to/input_table.tsv \
   --run_liftover T \
   --run_colocalization T \
   --outdir /path/to/output_directory \
   -w ./work \
```

| Parameter              | Description                                           |
| ---------------------- | ----------------------------------------------------- |
| `--summarystats_input` | Path to input table                                   |
| `--run_liftover`       | Whether to lift input data from hg37 to hg38          |
| `--run_colocalization` | Whether to follow-up fine-mapping with colocalization |
| `--outdir`             | Directory where output will be generated              |

## Essential configuration

| Parameter                        | Description                                                               |
| -------------------------------- | ------------------------------------------------------------------------- |
| **Input data**                   |                                                                           |
| `--summarystats_input`           | Input table with summary statistics and parameters for mungin and finemap |
| `--coloc_input`                  | Input file for coloc                                                      |
| `--coloc_id`                     | ID label for the coloc analysis output                                    |
| **Output settings**              |                                                                           |
| `--outdir`                       | Directory where output will be generated                                  |
| **Munging/finemapping settings** |                                                                           |


Check `nextflow.config` for a full list of customizabile parameters

## Input Data

### Fine-mapping

- GWAS summary statistics (.csv/.csv.gz/.tsv/.tsv.gz).
- Genotypes for LD reference panel in Plink .bed/.bim/.fam format. Providing, if possible, in sample LD greatly improves the accuracy of fine-mapping
- Metadata and GWAS-specific parameters table.

### Colocalisation

- csAnnData

## Test 

You can perform a test run to ensure everything is set up correctly. This will use a small test dataset and run the pipeline using Singularity.

```bash
nextflow run Biostatistics-Unit-HT/Flanders -profile test,singularity -w ./work
```