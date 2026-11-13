# Flanders: Finemapping coLocalization AND plEiotRopy Solver

Flanders is a Nextflow pipeline for fine-mapping and colocalization of genetic association signals across multiple studies and traits. It supports three stages:
1. **Fine-mapping**: harmonizes GWAS summary statistics and uses an LD reference panel to identify credible sets of candidate causal variants using [SuSiE](https://stephenslab.github.io/susieR/index.html).
2. **Colocalization**: tests pairs of credible sets for evidence of a shared causal variant.
3. **Module identification and joint fine-mapping**: groups colocalizing signals into modules and combines their evidence to refine candidate causal variants.
The stages can be run together, or existing fine-mapping results can be reused for downstream analyses.

## Requirements

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) v24.04+.
- One execution environment: [Docker](https://www.docker.com/), [Singularity](https://docs.sylabs.io/), or [Conda](https://docs.conda.io/en/latest/).


## Running the pipeline

Colocalization and joint fine-mapping are both enabled by default.

| Mode | Input | Stage settings |
|------|-------|----------------|
| Full workflow | `--summarystats_input` | Default settings |
| Fine-mapping only | `--summarystats_input` | `--run_colocalization false --run_joint_finemapping false` |
| Colocalization only | `--coloc_h5ad_input` | `--run_colocalization true --run_joint_finemapping false` |
| Colocalization and joint fine-mapping from existing results | `--coloc_h5ad_input` | Default settings |

For colocalization-only runs, omit `--summarystats_input`. For fine-mapping-only runs, omit `--coloc_h5ad_input`. Joint fine-mapping requires `--run_colocalization true`; the combination `--run_joint_finemapping true --run_colocalization false` is rejected.

### Full workflow

These examples use `-profile singularity`; replace it with `docker` or `conda` as appropriate. For local edits, replace `Biostatistics-Unit-HT/Flanders` with `/path/to/Flanders/main.nf`.


```bash
nextflow run Biostatistics-Unit-HT/Flanders \
    -profile singularity \
    --summarystats_input /path/to/input_table.tsv \
    --finemap_id my_finemap_run \
    --coloc_id my_coloc_run \
    --outdir flanders_output \
    -w ./work \
    -resume
```

### Fine-mapping only

```bash
nextflow run Biostatistics-Unit-HT/Flanders \
    -profile singularity \
    --summarystats_input /path/to/input_table.tsv \
    --finemap_id my_finemap_run \
    --run_colocalization false \
    --run_joint_finemapping false \
    --outdir flanders_output \
    -w ./work \
    -resume
```

### Colocalization only

```bash
nextflow run Biostatistics-Unit-HT/Flanders \
    -profile singularity \
    --coloc_h5ad_input /path/to/finemapping_output.h5ad \
    --run_joint_finemapping false \
    --coloc_id my_coloc_run \
    --outdir flanders_output \
    -w ./work \
    -resume
```

To include module identification and joint fine-mapping in this last example, remove `--run_joint_finemapping false`.

### Test dataset

```bash
nextflow run Biostatistics-Unit-HT/Flanders \
    -profile test,singularity \
    -w ./work \
    -resume
```

## Inputs

### Fine-mapping

Pass a tab-separated metadata table (.tsv) through `--summarystats_input` (see [example input table](example_data/example_data_finemapping_input_table.tsv)); the full column specification is documented in the [fine-mapping inputs wiki](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/Fine%E2%80%90mapping-inputs).

| Input | Required format |
|-------|-----------------|
| Summary statistics | TSV or CSV, optionally gzipped; paths and column mappings supplied in the metadata table |
| LD reference panel | A complete PLINK 1 `.bed/.bim/.fam` or PLINK 2 `.pgen/.pvar/.psam` file set |
| Metadata table | TSV with one row per input dataset and its study identifiers, column mappings, reference panel, and analysis settings |

Use absolute paths in the metadata table. In `bfile`, provide the reference-panel prefix **without a file extension**. PLINK 2 takes precedence if complete PLINK 1 and PLINK 2 sets exist for the same path and prefix.

Key settings are supplied **per row in the metadata table** (full column specification documented in the [fine-mapping inputs wiki](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/Fine%E2%80%90mapping-inputs)), rather than as pipeline command-line parameters:

| Column(s) | What to specify |
|-----------|-----------------|
| `input`, `study_id` | Summary-statistics path and study identifier |
| `chr_lab`, `pos_lab`, `rsid_lab`, `a1_lab`, `a0_lab`, `freq_lab`, `n_lab`, `effect_lab`, `se_lab`, `pvalue_lab` | Names of the corresponding summary-statistics columns; use `NA` for supported missing fields. |
| `bfile`, `grch`, `grch_bfile` | LD reference prefix and genome builds for summary statistics and reference panel (`37` or `38`). If `grch_bfile` is absent or empty, the pipeline assumes the LD reference uses the build specified in `grch`. Always provide `grch_bfile` when the reference and summary-statistics builds differ. |
| `process_bfile` | Whether to preprocess the LD reference. Preprocess the LD reference by standardizing variant IDs with alphabetically ordered alleles (**"chr"CHR:POS:EA:NEA** where EA is the first allele in alphabetical order), removing duplicate variants/positions, and optionally lifting GRCh37 coordinates to GRCh38. Set to `FALSE` only if the reference already meets these requirements. Processing is enabled when this column is omitted, and is forced for GRCh37 references when `--run_liftover true`. |
| `is_molQTL`, `key` | Whether the input contains multiple molecular phenotypes and the column identifying each phenotype |
| `type`, `sdY`, `s` | Trait type (`quant` or `cc`), quantitative-trait standard deviation, and case proportion for case-control traits. For quantitative traits, `sdY = 1` denotes a standardized phenotype; if `NA`, it will be estimated. Molecular-QTL inputs can use a per-phenotype `sdY` file with `phenotype_id` and `sdY` columns, covering every analysed phenotype with numeric, non-missing values. For case-control traits, provide `s` and set `sdY` to `NA`. |
| `maf` | Minor allele-frequency filter |
| `p_thresh1`, `p_thresh2`, `hole` | Locus-detection significance threshold, limit/grouping threshold, and maximum gap between grouped SNPs in base pairs |
| `cs_thresh` | Requested initial SuSiE credible-set coverage, e.g. `0.99`; see **Credible-set contents** below |



### Colocalization and joint fine-mapping

Supply fine-mapping AnnData through `--coloc_h5ad_input`. **Files must follow the [Flanders AnnData specification](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/csAnnData-specifications)**, an arbitrary `.h5ad` is not sufficient.

Multiple files can be selected using a **quoted glob**:

```bash
--coloc_h5ad_input '/path/to/finemapping_runs/*.h5ad'
```
When multiple fine-mapping AnnData files are supplied, their credible sets are concatenated into a single AnnData using the union of SNPs. The combined file is saved as `results/anndata/coloc_input_anndata.h5ad` and used for downstream analyses.

You can provide both `--summarystats_input` and `--coloc_h5ad_input` to combine newly generated fine-mapping results with existing results (inputs must use compatible genome builds).

Joint fine-mapping receives the complete pair-wise results table and the fine-mapping AnnData automatically from colocalization.

## Parameters


### Fine-mapping

| Parameter | Default | Effect |
|-----------|---------|--------|
| `--finemap_id` | `finemap_run` | Prefix for the fine-mapping AnnData filename |
| `--chromosomes` | `1-22` | Chromosomes to analyse; accepts lists/ranges such as `'1-6,8,10-22'`. X and Y are also accepted. |
| `--run_liftover` | `true` | Lift GRCh37 inputs to GRCh38 when needed. If disabled, summary statistics and LD references must already use matching builds. |
| `--large_locus_size` | `2000000` | Loci with size **greater than or equal to** this value in bp are excluded from fine-mapping and reported. |
| `--susie_max_iter` | `400` | Maximum SuSiE IBSS iterations |

Additional post-fit QC controls are `--susie_qc_cs_bf_thr` (`3`), `--susie_qc_pval_thr` (`1`), `--susie_qc_mean_r2_thr` (`0`), and `--susie_qc_min_r2_thr` (`0`). These filter credible sets by Bayes factor, top-SNP p-value, and purity criteria.

### Colocalization, modules, and joint fine-mapping

| Parameter | Default | Effect |
|-----------|---------|--------|
| `--run_colocalization` | `true` | Enable colocalization; supplying `--coloc_h5ad_input` also invokes this stage. |
| `--run_joint_finemapping` | `true` | Enable module identification, refinement and joint fine-mapping. **Requires colocalization.** |
| `--coloc_id` | `coloc_run` | Prefix for colocalization and final module output filenames |
| `--coloc_batch_size` | `2000` | Number of pair-wise tests per colocalization task |
| `--pph4_threshold` | `0.75` | Threshold for the filtered H4 table and module identification; also passed to the co-regulated module refinement procedure |
| `--pph3_threshold` | `0.75` | Threshold for the filtered H3 table; does not define module membership |
| `--module_refining_method` | `hyprcoloc_susie` | Refinement of partly connected modules: `hyprcoloc_susie` (default) or `coloc_percolation`. Fully connected modules bypass refinement; partly connected modules use the selected refinement method. |
| `--chunk_size` | `500` | Modules per joint fine-mapping task; positive integer |
| `--coloc_guide_table` | Unset | External CSV listing pairs to test; reuse the generated guide-table format. Overrides automatic guide generation and both exclusion settings below. |
| `--coloc_exclude_studies_table` | Unset | TSV with `study_id` and `phenotype_id`; exclude tests where either trait is listed. |
| `--coloc_skip_previous_studies` | `false` | When existing AnnData is supplied, exclude tests where both traits belong to those existing inputs; retain old–new tests. |

Exclusion settings and custom guide tables restrict the pair-wise evidence available for downstream module identification.

## Input transformations and output interpretation

### Variant handling and excluded loci

- In-sample LD is recommended. Only variants shared by the summary statistics and LD reference are used for fine-mapping.
- Harmonized variant IDs use `chrCHR:POS:EA:NEA`. Alleles are alphabetically ordered, with the first allele used as the effect allele; effect sizes and frequencies are adjusted accordingly. These are **not REF/ALT identifiers**.
- Loci overlapping the HLA interval are excluded: GRCh38 `chr6:28510120-33480577`, or GRCh37 `chr6:28477797-33448354`.
- Inspect `results/finemapping_exceptions/` for skipped loci and loci re-fitted with `L=1`. A completed pipeline does not imply that every input locus produced credible sets.

### Credible-set contents

**Fine-mapping outputs use expanded credible sets.** After SuSiE fitting and QC, each retained set contains the original credible-set SNPs **plus additional SNPs selected by an automatically determined threshold on that signal's SNP-level lABFs**. Original members are retained; this is not a fixed genomic-window extension.

`cs_thresh` controls the initial SuSiE coverage, not the final expanded membership. Recovery attempts can also reduce the requested coverage when no set is returned. Consequently, exported sets should not be interpreted as unmodified, fixed 99% SuSiE credible sets.

The expanded membership is written to the credible-set information tables and used to select lABFs stored in fine-mapping AnnData. It therefore also affects which pairs share SNPs and are considered for automatic colocalization. SNPs outside the stored membership have lABF imputed from the residual-lABF metadata during colocalization; sparse zeros in `X` should thus not be interpreted as measured lABFs for every SNP in the locus.

### AnnData outputs

| File | Observations (`obs`) | Variables (`var`) | Values |
|------|----------------------|-------------------|--------|
| Fine-mapping AnnData | Credible sets | SNPs | `X`: lABFs for stored expanded-set members; `beta` and `se` layers hold corresponding estimates |
| Module AnnData | Modules (`M_*`, also stored in `cs_name`) | SNPs | `X`: original lABFs for single-CS modules and joint lABFs for multi-CS modules |

Module metadata include `list_of_cs`, `list_of_study_id`, `list_of_phenotype_id`, `n_cs`, genomic coordinates, `most_likely_snp`, `a1`, `a0`, `list_of_most_likely_beta`, and `list_of_most_likely_se`. Effect-estimate lists describe contributing traits, not a single module-level beta or standard error. The master TSV retains credible-set-level rows and their module assignments.

## Outputs

Analysis paths below are relative to `<outdir>/results/`; `--outdir` defaults to `flanders_output`.

| Path | Contents |
|------|----------|
| `gwas_and_loci_tables/*_dataset_aligned.tsv.gz` | Harmonized, optionally lifted summary statistics |
| `gwas_and_loci_tables/*_loci.tsv` | Detected locus boundaries and sentinel-SNP information |
| `finemapping_exceptions/` | Skipped loci and SuSiE recovery reports; see the [exceptions wiki](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/Fine%E2%80%90mapping-exceptions) |
| `finemap/*_susie_finemap.rds` | Per-locus results; published with `--publish_susie true` |
| `coloc_info_tables/*_coloc_info_master_table.tsv` | Credible-set membership and summary metadata per study |
| `anndata/<finemap_id>_finemap_results.h5ad` | Fine-mapping AnnData |
| `anndata/coloc_input_anndata.h5ad` | Combined input AnnData; published when multiple AnnData inputs are merged |
| `coloc/coloc_guide_table.csv` | Automatically generated pairs to test; not regenerated when an external guide is supplied |
| `coloc/<coloc_id>_colocalization.table.all.tsv` | Complete results for the tested pairs |
| `coloc/<coloc_id>_colocalization.table.H4.tsv` | Results passing the H4 threshold |
| `coloc/<coloc_id>_colocalization.table.H3.tsv` | Results passing the H3 threshold |
| `master_table/<coloc_id>_master_table.tsv` | Credible-set-level module assignments and joint fine-mapping results |
| `master_table/<coloc_id>_modules_JF.h5ad` | Module-level AnnData |

Module outputs require joint fine-mapping to be enabled.

## Execution and resource settings

- Use `--coloc_batch_size` and `--chunk_size` to adjust task granularity. Smaller batches create more tasks; larger batches put more work in each task.
- Configure the executor and process resources in a separate config passed with `-c /path/to/site.config`. Inspect the trace and execution report before adjusting resources. Processes use the resource labels defined in `conf/base.config`.
- Automatic guide generation tests pairs sharing stored credible-set SNPs. Check the generated guide and result counts, especially after exclusions or when combining datasets.


## Convert Open Targets credible sets to AnnData

Flanders’ modular design allows precomputed fine-mapping results to be used directly for colocalization and joint fine-mapping. Use `bin/parse_credible_sets_refined.py` to convert Open Targets study and credible-set parquet files in the AnnData format required by the pipeline.

| Argument | Value |
|----------|-------|
| `--study` | Study parquet file |
| `--credset` | Credible-set parquet file |
| `--ref_panel` | Reference CSV with `chr`, `pos`, `a0`, and `a1`; used to add missing variants within the selected credible-set intervals. Allele columns can be in either order. |
| `--project` | Project identifier, e.g. `FINNGEN_R12`, `GTEx`, or `GCST` |
| `--type_sumstat` | `gwas` (default), `eqtl`, or `pqtl` |
| `--type_finemap` | `SuSie` (default) or `SuSiE-inf` |
| `--out` | Output prefix; writes `<out>.h5ad` |

```bash
python bin/parse_credible_sets_refined.py \
    --project FINNGEN_R12 \
    --type_sumstat gwas \
    --type_finemap SuSie \
    --study /path/to/study.parquet \
    --credset /path/to/credible_set.parquet \
    --ref_panel /path/to/ref_panel.csv \
    --out /path/to/ot_credsets
```

Pass the resulting file as `--coloc_h5ad_input /path/to/ot_credsets.h5ad` in the colocalization example above.

## Credits
Developed by the Biostatistics and Genome Analysis Units at [Human Technopole](https://humantechnopole.it/en/)<br>
-  [Arianna Landini](mailto:arianna.landini@fht.org)<br>
-  [Sodbo Sharapov](mailto:sodbo.sharapov@fht.org)<br>
-  [Edoardo Giacopuzzi](mailto:edoardo.giacopuzzi@fht.org)<br>
-  [Bruno Ariano](mailto:bruno.ariano@fht.org)<br>
-  [Nicola Pirastu](mailto:nicola.pirastu@fht.org)<br>
