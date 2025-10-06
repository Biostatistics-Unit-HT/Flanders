#!/usr/bin/env python3
import bioalpha as bsc
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Concatenate multiple .h5ad files listed in a text file into one AnnData object."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to a text file listing .h5ad files to concatenate (one per line)."
    )
    parser.add_argument(
        "-o", "--output_file",
        required=True,
        help="Path and file name for the output merged .h5ad file."
    )
    args = parser.parse_args()

    # Read list of files
    with open(args.input_txt, "r") as f:
        h5ad_files = [line.strip() for line in f if line.strip()]

    if not h5ad_files:
        raise ValueError("No valid file paths found in the input text file.")

    print(f"Loading {len(h5ad_files)} .h5ad files...")
    adata_list = [bsc.sc.read_h5ad(file) for file in h5ad_files]

    print("Concatenating...")
    adata = bsc.sc.concat(adata_list, join='outer', fill_value=0)

    print(f"Writing output to: {args.output_h5ad}")
    adata.write_h5ad(args.output_h5ad)

    print("✅ Merge complete!")

if __name__ == "__main__":
    main()
