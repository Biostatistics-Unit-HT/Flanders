#!/usr/bin/env python3
import bioalpha as bsc
import argparse
import pandas as pd
import re

def fix_ad_var(adata):
    """
    Standardize variant metadata in an AnnData object.

    Ensures:
      - 'snp' column exists in adata.var
      - Extracts chromosome ('chr') and position ('pos') from SNP identifiers
        in format 'chr{num}:{pos}:EA:RA'
    """
    var = adata.var.copy()

    # Add 'snp' column if missing
    if "snp" not in var.columns:
        var["snp"] = var.index.astype(str)

    # Extract chr and pos using regex
    var["chr"] = var["snp"].apply(lambda x: re.sub(r":.*", "", x))
    var["pos"] = var["snp"].apply(lambda x: int(re.sub(r".*:(\d+):.*", r"\1", x)))

    adata.var = var
    return adata


def main():
    parser = argparse.ArgumentParser(
        description="Concatenate multiple .h5ad files listed in a text file into one AnnData object."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to a text file listing .h5ad files to concatenate (one per line).",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        required=True,
        help="Path and file name for the output merged .h5ad file.",
    )
    args = parser.parse_args()

    # Read list of files
    with open(args.input, "r") as f:
        h5ad_files = [line.strip() for line in f if line.strip()]

    if not h5ad_files:
        raise ValueError("No valid file paths found in the input text file.")

    print(f"Loading {len(h5ad_files)} .h5ad files...")
    adata_list = [bsc.sc.read_h5ad(file) for file in h5ad_files]

    print("Concatenating...")
    adata = bsc.sc.concat(adata_list, join="outer", fill_value=0)

    print("Fixing variant metadata...")
    adata = fix_ad_var(adata)

    print(f"Writing output to: {args.output_file}")
    adata.write_h5ad(args.output_file)

    print("✅ Merge complete!")


if __name__ == "__main__":
    main()
