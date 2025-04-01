import pandas as pd
import argparse
import os
import sys

def main():
    parser = argparse.ArgumentParser(description="Parsing columns for validation.")
    parser.add_argument("--launchdir", type=str, required=True, help="Launchdir for path")
    parser.add_argument("--liftover", action='store_true', help="Option for liftover")
    parser.add_argument("--table", type=str , required=True, help="Path of the list of inputs")


    args = parser.parse_args()
    def fill_path(gwas_path):
        if os.path.isabs(gwas_path):
            gwas_path
        else:
            gwas_path = args.launchdir + "/" + gwas_path
        return gwas_path
    
    def read_gwas_file(gwas_path):
        if gwas_path.endswith('.gz'):
            return pd.read_csv(gwas_path, compression='gzip', sep='\t', nrows=1)
        else:
            return pd.read_csv(gwas_path, sep='\t', nrows=1)

    
    table_files = pd.read_csv(fill_path(args.table), sep = "\t")

    # Check that the builds are consistent
    if args.liftover:
        if len(set(table_files[["grch"]]))>1:
            print("Error: multiple GRCh build detected in input file but run_liftover is false.", file=sys.stderr)
            return sys.exit(1)
           
    for index, record in table_files.iterrows():
        gwas_path = fill_path(record["input"])
        bfile_path = fill_path(record["bfile"])
        # Check that the GWAS exist
        if os.path.exists(gwas_path):
            # Check that the bfiles exist
            if os.path.exists(bfile_path + ".bim") and os.path.exists(bfile_path + ".bed") and os.path.exists(bfile_path + ".fam"):

                gwas_sub = read_gwas_file(gwas_path)
                # Check that the columns are found
                gwas_col = list(gwas_sub.columns)
                mandatory_columns = ["pos_lab", "rsid_lab", "chr_lab", "a0_lab", "a1_lab", "effect_lab", "se_lab"]
                if all(record[col] in gwas_col for col in mandatory_columns):
                    continue
                else:
                    print("Error: The columns from the GWAS are not contained in the input table.", file=sys.stderr)
                    return sys.exit(1)
            else:
                print(f"One of the bfiles {bfile_path} was not found.", file=sys.stderr)
                return sys.exit(1)
        else:
            print(f"The input GWAS {gwas_path} file was not found.", file=sys.stderr)
            return sys.exit(1)
    return 0

if __name__ == "__main__":
    main()
