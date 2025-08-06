# Flanders Pipeline

**Flanders : Finemapping coLocalization AND plEiotRopy Solver BETA**

The pipeline combines a set of tools to efficiently colocalise association signals across large sets of traits. The colocalization process is mainly composed of two steps: first the finemapping of associated loci and then colocalization itself. The pipeline optimize this process to allow for the analysis of large datasets.

## Main features

Describe here

## How to use

When you have Nextflow installed, you can run the pipeline directly from our repository using:

```bash
nextflow run Biostatistics-Unit-HT/Flanders [-r <version/commit>] \
   -profile [docker|singularity|conda] \
   --summarystats_input /path/to/input_table.tsv \
   --run_liftover T \
   --run_colocalization T \
   --outdir /path/to/output_directory \
   -w ./work \
```

Please read the [Quick Start Guide](quick-start.md) and the [input data section](inputs.md) for more details on how to prepare your input data and configure the pipeline.