#!/usr/bin/env python

import anndata as an
import argparse

def main():
    parser = argparse.ArgumentParser(description="Extract study and phenotype IDs from AnnData object.")
    parser.add_argument("--input", required=True, help="Path to a file containing a list of anndata files.")
    parser.add_argument("--output", required=True, help="Path to the output file to save study and phenotype IDs.")
    args = parser.parse_args()

    with open(args.input, "r") as f:
        input_files = f.read().splitlines()

    dfs = []
    for f in input_files:
        print(f"Processing {f}")
        adata = an.read_h5ad(f
        out_df = adata.obs[["study_id", "phenotype_id"]].drop_duplicates()
        print(f"Number of extracted study-phenotype pairs: {len(out_df)}")
        dfs.append(out_df)
    
    print("Concatenating dataframes")
    out_df = pd.concat(dfs)
    out_df = out_df.drop_duplicates()

    print(f"Total number of extracted study-phenotype pairs: {len(out_df)}")
    out_df.to_csv(args.output, sep="\t", index=False)
    print(f"Study-phenotype pairs saved to {args.output}")

if __name__ == "__main__":
    main()
