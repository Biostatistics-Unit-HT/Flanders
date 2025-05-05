#!/usr/bin/env python

import pandas as pd
import argparse
import os
import sys

def fill_path(gwas_path, base_dir):
    if os.path.isabs(gwas_path):
        gwas_path
    else:
        gwas_path = base_dir + "/" + gwas_path
    return gwas_path

def read_gwas_file(gwas_path):
    if gwas_path.endswith('.gz'):
        return pd.read_csv(gwas_path, compression='gzip', sep='\t', nrows=1)
    else:
        return pd.read_csv(gwas_path, sep='\t', nrows=1)

def main():
    parser = argparse.ArgumentParser(description="Parsing columns for validation.")
    parser.add_argument("--launchdir", type=str, required=True, help="Launchdir for path")
    parser.add_argument("--liftover", action='store_true', help="Option for liftover")
    parser.add_argument("--table", type=str , required=True, help="Path of the list of inputs")
    args = parser.parse_args()

    table_files = pd.read_csv(args.table, sep = "\t")

    # Check that the builds are consistent
    if "grch_bfile" in table_files.columns:
        if args.liftover == False:
            if len(set(table_files[["grch_bfile"]]))>1:
                print("Error: multiple GRCh build detected in grch_bfile column but run_liftover is false.", file=sys.stderr)
                return sys.exit(1)
            if not (table_files['grch_bfile'] == table_files['grch']).all():
                print("Error: grch_bfile and grch columns do not match.", file=sys.stderr)
                return sys.exit(1)
    
    if args.liftover == False:
            if len(set(table_files[["grch"]]))>1:
                print("Error: multiple GRCh build detected in grch column but run_liftover is false.", file=sys.stderr)
                return sys.exit(1)
    
            
           
    for index, record in table_files.iterrows():
        gwas_path = fill_path(record["input"], args.launchdir)
        bfile_path = fill_path(record["bfile"], args.launchdir)
        # Check that the GWAS exist
        if not os.path.exists(gwas_path):
            print(f"The input GWAS {gwas_path} file was not found.", file=sys.stderr)
            return sys.exit(1)
        # Check that the bfiles exist
        if not all([os.path.exists(bfile_path + ".bim"), os.path.exists(bfile_path + ".bed"), os.path.exists(bfile_path + ".fam")]):
            print(f"One of the bfiles {bfile_path} was not found.", file=sys.stderr)
            return sys.exit(1)

        # Check that the desired columns are found in sumstats file
        gwas_sub = read_gwas_file(gwas_path)
        
        gwas_col = list(gwas_sub.columns)
        all_labels = ["pos_lab", "rsid_lab", "chr_lab", "a0_lab", "a1_lab", "effect_lab", "se_lab", "freq_lab", "n_lab"]
        if record['is_molQTL'] == "TRUE": all_labels.append('key')
        label_cols_in_df = [col for col in all_labels if col in table_files.columns]
        col_names_to_scan = [record[col] for col in label_cols_in_df if pd.notna(record[col])]

        print (f"Checking columns {col_names_to_scan} in file {gwas_path}.")
        if all(col in gwas_col for col in col_names_to_scan):
            continue
        else:
            print(f"Error: One of the columns defined in the input table are not present in file {gwas_path}.", file=sys.stderr)
            return sys.exit(1)

if __name__ == "__main__":
    main()
